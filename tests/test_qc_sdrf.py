"""Standards-pending QC-to-SDRF proposal and existing-refinement regression tests."""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np

from prideqc.metrics import QCMetricCalculator, RunSummary
from prideqc.models import AnalysisResult, Metric, RunMetadata, Spectrum
from prideqc.pipeline import FileOutcome, Workflow, WorkflowOptions
from prideqc.qc_sdrf import (
    QC_SDRF_PROPOSALS,
    collect_qc_sdrf_proposals,
    metric_annotations,
    write_qc_sdrf_proposals,
)
from prideqc.sdrf import SDRFDocument
from prideqc.validation import ValidationReport


def _result(
    name: str,
    ms1: object = 0,
    ms2: object = 0,
    duration: object = 0.0,
) -> AnalysisResult:
    return AnalysisResult(
        input_path=Path(name),
        metadata=RunMetadata(),
        metrics=[
            Metric("NumberOfSpectra_MS1", ms1),
            Metric("NumberOfSpectra_MS2", ms2),
            Metric("ChromatographyDuration", duration),
        ],
        annotations=[],
        warnings=[],
        engine_version="test",
        elapsed_seconds=1.0,
    )


class QCSDRFProposalTests(unittest.TestCase):
    def test_calculator_retains_empty_scans_and_all_ms_level_rt_span(self) -> None:
        run = RunSummary()
        for level, rt in [(1, 30.0), (2, 6.0), (2, 126.0), (3, float("nan"))]:
            run.consume_spectrum(
                Spectrum(
                    ms_level=level,
                    rt=rt,
                    mz=np.array([], dtype=float),
                    intensity=np.array([], dtype=float),
                )
            )
        metrics = {m.key: m.value for m in QCMetricCalculator().calculate(run)}
        self.assertEqual(metrics["NumberOfSpectra_MS1"], 1)
        self.assertEqual(metrics["NumberOfSpectra_MS2"], 2)
        self.assertEqual(metrics["ChromatographyDuration"], 120.0)

    def test_calculator_missing_rt_and_single_scan_zero_span(self) -> None:
        run = RunSummary()
        run.consume_spectrum(
            Spectrum(
                ms_level=2,
                rt=float("inf"),
                mz=np.array([], dtype=float),
                intensity=np.array([], dtype=float),
            )
        )
        values = {m.key: m.value for m in QCMetricCalculator().calculate(run)}
        self.assertEqual(values["NumberOfSpectra_MS1"], 0)
        self.assertEqual(values["NumberOfSpectra_MS2"], 1)
        self.assertIsNone(values["ChromatographyDuration"])
        run.consume_spectrum(
            Spectrum(
                ms_level=1,
                rt=25.0,
                mz=np.array([], dtype=float),
                intensity=np.array([], dtype=float),
            )
        )
        values = {m.key: m.value for m in QCMetricCalculator().calculate(run)}
        self.assertEqual(values["ChromatographyDuration"], 0.0)

    def test_ontology_mapping_and_zero_counts_duration(self) -> None:
        annotations = metric_annotations(_result("zero.raw"))
        self.assertEqual(
            [item.value["psi_ms_accession"] for item in annotations],
            ["MS:4000059", "MS:4000060", "MS:4000067"],
        )
        self.assertEqual([item.sdrf_value for item in annotations], ["0", "0", "0 s"])
        self.assertEqual(
            [item.value["unit_accession"] for item in annotations],
            [
                "UO:0000189",
                "UO:0000189",
                "UO:0000010",
            ],
        )
        self.assertEqual(len(QC_SDRF_PROPOSALS), 3)

    def test_missing_invalid_nonfinite_never_coerced_to_zero(self) -> None:
        for invalid in [None, float("nan"), float("inf"), float("-inf"), -1, True, "2"]:
            with self.subTest(invalid=invalid):
                result = _result("test.raw", ms1=invalid, ms2=invalid, duration=invalid)
                self.assertFalse(metric_annotations(result))
                rows = collect_qc_sdrf_proposals([result])
                self.assertTrue(all(not row.value for row in rows))

    def test_float_count_rejected_but_float_duration_accepted(self) -> None:
        annotations = metric_annotations(_result("test.raw", 1.0, 3, 12.5))
        self.assertEqual([a.sdrf_value for a in annotations], ["3", "12.5 s"])

    def test_duplicate_metrics_fail_closed(self) -> None:
        result = _result("run.raw", 5, 10, 50.0)
        result.metrics.append(Metric("NumberOfSpectra_MS1", 999))
        self.assertEqual(collect_qc_sdrf_proposals([result])[0].status, "ambiguous_metric")
        self.assertEqual(len(metric_annotations(result)), 2)

    def test_repeated_rows_and_distinct_runs_not_averaged(self) -> None:
        sdrf = SDRFDocument(
            ["comment[data file]", "factor value[treatment]"],
            [["a.raw", "X"], ["a.raw", "Y"], ["b.raw", "Z"]],
        )
        before = [row.copy() for row in sdrf.rows]
        records = collect_qc_sdrf_proposals(
            [
                _result("a.raw", 0, 8, 100.25),
                _result("b.raw", 9, 0, 300.5),
            ],
            sdrf,
        )
        self.assertEqual(len(records), 9)
        self.assertEqual([r.value for r in records[::3]], ["0", "0", "9"])
        self.assertEqual([r.value for r in records[2::3]], ["100.25 s", "100.25 s", "300.5 s"])
        self.assertEqual(sdrf.rows, before)
        self.assertEqual(sdrf.columns, ["comment[data file]", "factor value[treatment]"])

    def test_alias_exact_case_and_unmatched(self) -> None:
        sdrf = SDRFDocument(["comment[data file]"], [["ALIAS.raw"], ["Run.raw"], ["run.raw"]])
        results = [_result("Run.raw", 4, 5, 6.0)]
        unmapped = collect_qc_sdrf_proposals(results, sdrf)
        self.assertEqual(
            [r.status for r in unmapped[::3]], ["unmatched", "proposed_unregistered", "unmatched"]
        )
        mapped = collect_qc_sdrf_proposals(results, sdrf, {"ALIAS.raw": "Run.raw"})
        self.assertEqual(
            [r.status for r in mapped[::3]],
            ["proposed_unregistered", "proposed_unregistered", "unmatched"],
        )
        with self.assertRaisesRegex(ValueError, "not analyzed"):
            collect_qc_sdrf_proposals(results, sdrf, {"alias.raw": "absent.raw"})

    def test_ambiguous_basename_mapping_raises(self) -> None:
        sdrf = SDRFDocument(["comment[data file]"], [["same.raw"]])
        with self.assertRaisesRegex(ValueError, "Ambiguous data-file mapping"):
            collect_qc_sdrf_proposals([_result("/one/same.raw"), _result("/two/same.raw")], sdrf)

    def test_existing_conflict_duplicate_column_and_missing(self) -> None:
        sdrf = SDRFDocument(
            ["comment[data file]", "comment[number of MS1 spectra]"],
            [["one.raw", "4"], ["two.raw", "not available"], ["three.raw", "0"]],
        )
        original = [r.copy() for r in sdrf.rows]
        records = collect_qc_sdrf_proposals(
            [_result("one.raw", 9), _result("two.raw", 0), _result("three.raw", 0)], sdrf
        )
        self.assertEqual(
            [r.status for r in records[::3]],
            [
                "existing_conflict_unregistered",
                "existing_missing_unregistered",
                "existing_match_unregistered",
            ],
        )
        self.assertEqual(sdrf.rows, original)
        repeated = SDRFDocument(
            [
                "comment[data file]",
                "comment[number of MS1 spectra]",
                "comment[number of MS1 spectra]",
            ],
            [["one.raw", "4", "5"]],
        )
        self.assertEqual(
            collect_qc_sdrf_proposals([_result("one.raw")], repeated)[0].status, "ambiguous_column"
        )

    def test_unavailable_duration_preserved_and_json_sidecar_idempotent(self) -> None:
        sdrf = SDRFDocument(["comment[data file]"], [["one.raw"]])
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            records = write_qc_sdrf_proposals(root, [_result("one.raw", 0, 2, None)], sdrf)
            self.assertEqual(records[2].status, "unavailable_metric")
            self.assertEqual(records[2].value, "")
            payload = json.loads((root / "sdrf-qc-proposals.json").read_text())
            self.assertEqual(payload["status"], "standards_pending_not_written_to_sdrf")
            self.assertEqual(payload["mapping"][0]["accession"], "MS:4000059")
            old_tsv = (root / "sdrf-qc-proposals.tsv").read_bytes()
            write_qc_sdrf_proposals(root, [_result("one.raw", 0, 2, None)], sdrf)
            self.assertEqual((root / "sdrf-qc-proposals.tsv").read_bytes(), old_tsv)
            self.assertEqual(sdrf.columns, ["comment[data file]"])

    def test_analyze_cohort_flow_writes_sidecar_without_sdrf_mutation(self) -> None:
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            source = root / "measured.raw"
            source.touch()
            original = root / "study.sdrf.tsv"
            original.write_text("comment[data file]\tfactor value[condition]\nmeasured.raw\tA\n")
            result = _result("measured.raw", 0, 26, 3.5)
            validator = MagicMock()
            validator.validate.side_effect = [
                ValidationReport(str(original), "ms-proteomics", False, (), "mock"),
                ValidationReport("refined", "ms-proteomics", False, (), "mock"),
            ]
            output = root / "analysis"
            with patch(
                "prideqc.pipeline._analyze_file",
                return_value=FileOutcome(source, result),
            ):
                manifest = Workflow(
                    WorkflowOptions(refine_sdrf_qc=True),
                    validator=validator,
                ).run([source], output, sdrf=original)
            self.assertTrue(manifest["success"])
            self.assertEqual(
                manifest["cohort_refinement"]["run_qc_proposals"],
                "sdrf-qc-proposals.tsv",
            )
            self.assertTrue((output / "sdrf-qc-proposals.json").exists())
            refined = SDRFDocument.read(output / "refined.sdrf.tsv")
            self.assertEqual(refined.indices("comment[number of MS1 spectra]"), [])
            self.assertEqual(refined.indices("comment[chromatography duration]"), [])

    def test_existing_summaries_workflow_does_not_insert_unregistered_fields(self) -> None:
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            source = root / "summaries" / "PXD123456"
            source.mkdir(parents=True)
            result = _result("one.raw", 8, 0, 70.0)
            result.project_accession = "PXD123456"
            (source / "one.raw.summary.json").write_text(json.dumps(result.to_dict()))
            original = root / "study.sdrf.tsv"
            original.write_text("comment[data file]\tfactor value[treatment]\none.raw\tcontrol\n")
            validator = MagicMock()
            validator.validate.side_effect = [
                ValidationReport(str(original), "ms-proteomics", False, (), "mock"),
                ValidationReport("refined", "ms-proteomics", False, (), "mock"),
            ]
            output = root / "output"
            manifest = Workflow(
                WorkflowOptions(refine_sdrf_qc=True), validator=validator
            ).refine_existing_sdrf(
                source,
                output,
                sdrf=original,
            )
            self.assertTrue(manifest["success"])
            self.assertEqual(manifest["run_qc_proposals"], "sdrf-qc-proposals.tsv")
            refined = SDRFDocument.read(output / "refined.sdrf.tsv")
            self.assertFalse(any("number of MS" in c.lower() for c in refined.columns))
            self.assertNotIn("comment[chromatography duration]", refined.columns)
            proposals = collect_qc_sdrf_proposals([result], SDRFDocument.read(original))
            self.assertEqual([r.value for r in proposals], ["8", "0", "70 s"])
            self.assertEqual(
                original.read_text(),
                "comment[data file]\tfactor value[treatment]\none.raw\tcontrol\n",
            )


if __name__ == "__main__":
    unittest.main()
