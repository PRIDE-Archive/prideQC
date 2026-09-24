"""Cohort-level reanalysis and SDRF refinement contracts."""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import MagicMock, patch

from prideqc.cohort import (
    COHORT_FRAGMENT_FIELD,
    COHORT_MODIFICATION_FIELD,
    COHORT_PRECURSOR_FIELD,
    COHORT_PUTATIVE_MODIFICATION_FIELD,
    PUTATIVE_MODIFICATION_COLUMN,
    SemanticEvidence,
    read_semantic_evidence,
    synthesize_cohort,
)
from prideqc.models import AnalysisResult, Annotation, CVTerm, EvidenceKind, Metric, RunMetadata
from prideqc.pipeline import FileOutcome, Workflow, WorkflowOptions
from prideqc.sdrf import SDRFDocument
from prideqc.validation import ValidationReport


def _result(
    name: str,
    *,
    precursor: float = 8.0,
    fragment: float = 20.0,
    duration: float = 5400.0,
    mass_shift: float | None = None,
    candidates: list[tuple[str, str]] | None = None,
) -> AnalysisResult:
    annotations = [
        Annotation(
            "suggested_precursor_search_tolerance_ppm",
            {"unit": "ppm", "suggested_tolerance": precursor, "confidence": "high"},
            EvidenceKind.INFERRED,
            "test precursor precision",
            support=1000,
            total=1200,
        ),
        Annotation(
            "suggested_fragment_search_tolerance_ppm",
            {
                "unit": "ppm",
                "suggested_tolerance": fragment,
                "confidence": "high",
                "resolution_regime": "high-resolution",
            },
            EvidenceKind.INFERRED,
            "test fragment precision",
            support=10000,
            total=12000,
        ),
    ]
    if mass_shift is not None:
        candidate_rows = [
            {
                "unimod_accession": accession,
                "name": label,
                "candidate_category": "biological-ptm",
            }
            for accession, label in (candidates or [("UniMod:21", "Phosphorylation")])
        ]
        annotations.append(
            Annotation(
                "putative_modification_mass_shifts",
                [
                    {
                        "delta_mass_da": mass_shift,
                        "pair_support": 120,
                        "unique_spectrum_support": 80,
                        "classification": "putative-ptm",
                        "confidence": "high-support",
                        "unimod_candidates": candidate_rows,
                        "diagnostic_unimod_candidates": candidate_rows,
                    }
                ],
                EvidenceKind.INFERRED,
                "test mass shift",
            )
        )
    return AnalysisResult(
        input_path=Path(name),
        metadata=RunMetadata(instruments=[CVTerm("MS:1002416", "Orbitrap Fusion")]),
        metrics=[
            Metric("ChromatographyDuration", duration),
            Metric("IsolationWidth_MS2_Median", 1.6),
            Metric("NumberOfSpectra_MS1", 4000),
            Metric("NumberOfSpectra_MS2", 40000),
        ],
        annotations=annotations,
        warnings=[],
        engine_version="test",
        elapsed_seconds=1.0,
    )


class CohortRefinementTests(unittest.TestCase):
    def test_cohort_tolerances_use_one_conservative_setting_per_group(self) -> None:
        results = [
            _result("a.raw", precursor=7.1, fragment=20.1),
            _result("b.raw", precursor=8.2, fragment=21.2),
            _result("c.raw", precursor=9.4, fragment=24.1),
            _result("d.raw", precursor=8.8, fragment=22.5),
        ]
        synthesis = synthesize_cohort(results)
        self.assertEqual(len(synthesis.groups), 1)
        for result in results:
            precursor = next(
                item for item in result.annotations if item.field == COHORT_PRECURSOR_FIELD
            )
            fragment = next(
                item for item in result.annotations if item.field == COHORT_FRAGMENT_FIELD
            )
            self.assertEqual(precursor.sdrf_value, "10 ppm")
            self.assertEqual(fragment.sdrf_value, "25 ppm")
            self.assertEqual(precursor.value["common_max"], 9.4)
            self.assertEqual(fragment.value["common_max"], 24.1)
            self.assertEqual((precursor.support, precursor.total), (4, 4))
            self.assertEqual(precursor.value["source_measurement_support"], 4000)

    def test_chromatography_duration_can_separate_experiment_groups(self) -> None:
        results = [
            *[_result(f"short-{index}.raw", duration=5400 + index) for index in range(4)],
            *[_result(f"long-{index}.raw", duration=7200 + index) for index in range(4)],
        ]
        synthesis = synthesize_cohort(results)
        self.assertEqual(len(synthesis.groups), 2)
        short_groups = {synthesis.assignments[f"short-{index}.raw"] for index in range(4)}
        long_groups = {synthesis.assignments[f"long-{index}.raw"] for index in range(4)}
        self.assertEqual(len(short_groups), 1)
        self.assertEqual(len(long_groups), 1)
        self.assertNotEqual(short_groups, long_groups)

    def test_only_semantic_supported_recurrent_ptm_family_is_sdrf_eligible(self) -> None:
        results = [
            _result(f"run-{index}.raw", mass_shift=79.9663 + index * 1e-5)
            for index in range(10)
        ]
        evidence = {
            ("PXDTEST", "UNIMOD:21"): SemanticEvidence(
                "supported", ("publication",), ("phosphopeptide enrichment",)
            )
        }
        synthesis = synthesize_cohort(
            results, semantic_evidence=evidence, project_accession="PXDTEST"
        )
        self.assertEqual(len(synthesis.ptm_families), 1)
        self.assertEqual(synthesis.ptm_families[0]["unimod_accession"], "UniMod:21")
        self.assertEqual(synthesis.ptm_families[0]["semantic_evidence_status"], "supported")
        self.assertGreaterEqual(
            synthesis.ptm_families[0]["raw_prevalence_probability"], 0.99
        )
        for result in results:
            modification = next(
                item for item in result.annotations if item.field == COHORT_MODIFICATION_FIELD
            )
            self.assertEqual(
                modification.sdrf_value, "NT=Phosphorylation;AC=UniMod:21"
            )
            putative = next(
                item
                for item in result.annotations
                if item.field == COHORT_PUTATIVE_MODIFICATION_FIELD
            )
            self.assertIn(
                "prideqc putative modification: NT=Phosphorylation;AC=UniMod:21",
                putative.sdrf_value or "",
            )

    def test_mass_only_recurrent_ptm_family_is_review_only(self) -> None:
        results = [
            _result(f"run-{index}.raw", mass_shift=14.0155 + index * 1e-5,
                    candidates=[("UniMod:34", "Methylation")])
            for index in range(10)
        ]
        synthesis = synthesize_cohort(results, project_accession="PXD000612")
        self.assertEqual(synthesis.ptm_families, [])
        self.assertEqual(len(synthesis.ptm_review_families), 1)
        review = synthesis.ptm_review_families[0]
        self.assertEqual(review["unimod_accession"], "UniMod:34")
        self.assertEqual(review["semantic_evidence_status"], "not-evaluated")
        self.assertEqual(review["sdrf_status"], "hold-semantic-not-supported")
        self.assertFalse(
            any(
                item.field == COHORT_MODIFICATION_FIELD
                for result in results
                for item in result.annotations
            )
        )
        for result in results:
            putative = [
                item
                for item in result.annotations
                if item.field == COHORT_PUTATIVE_MODIFICATION_FIELD
            ]
            self.assertEqual(len(putative), 1)
            self.assertEqual(putative[0].sdrf_column, PUTATIVE_MODIFICATION_COLUMN)
            self.assertIn(
                "prideqc putative modification: NT=Methylation;AC=UniMod:34",
                putative[0].sdrf_value or "",
            )
            self.assertIn("DM=14.015", putative[0].sdrf_value or "")
            self.assertIn("STATUS=hold-semantic-not-supported", putative[0].sdrf_value or "")

    def test_mass_ambiguous_recurrent_family_is_not_written_as_modification(self) -> None:
        results = [
            _result(
                f"run-{index}.raw",
                mass_shift=79.9663,
                candidates=[
                    ("UniMod:21", "Phosphorylation"),
                    ("UniMod:99913", "Phosphorylation Decoy"),
                ],
            )
            for index in range(10)
        ]
        synthesis = synthesize_cohort(results)
        self.assertEqual(synthesis.ptm_families, [])
        self.assertFalse(
            any(
                item.field == COHORT_MODIFICATION_FIELD
                for result in results
                for item in result.annotations
            )
        )

    def test_large_chromatography_duration_gap_splits_before_multivariate_clustering(self) -> None:
        results = [
            *[_result(f"short-{index}.raw", duration=7200 + index) for index in range(4)],
            *[_result(f"long-{index}.raw", duration=15900 + index) for index in range(6)],
        ]
        synthesis = synthesize_cohort(results)
        self.assertEqual(len(synthesis.groups), 2)
        short_groups = {synthesis.assignments[f"short-{index}.raw"] for index in range(4)}
        long_groups = {synthesis.assignments[f"long-{index}.raw"] for index in range(6)}
        self.assertEqual(len(short_groups), 1)
        self.assertEqual(len(long_groups), 1)
        self.assertNotEqual(short_groups, long_groups)
        self.assertTrue(
            any("univariate chromatography duration split" in item for item in synthesis.evidence)
        )

    def test_study_evidence_reader_uses_v4_schema(self) -> None:
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "evidence.tsv"
            path.write_text(
                "pxd_accession\tunimod_accession\tevidence_status\tevidence_source\tevidence_note\n"
                "PXD000612\tUniMod:21\tsupported\tpublication\tphosphorylation study\n",
                encoding="utf-8",
            )
            evidence = read_semantic_evidence(path)
        item = evidence[("PXD000612", "UNIMOD:21")]
        self.assertEqual(item.status, "supported")
        self.assertEqual(item.sources, ("publication",))

    def test_new_sdrf_columns_are_inserted_before_factor_columns(self) -> None:
        result = _result("run.raw")
        result.annotations.extend(
            [
                Annotation(
                    "collision_energy_ms2",
                    "25 eV",
                    EvidenceKind.OBSERVED,
                    "mzML activation energy",
                    sdrf_column="comment[collision energy]",
                    sdrf_value="25 eV",
                ),
                Annotation(
                    COHORT_MODIFICATION_FIELD,
                    {"experiment_group": "Experiment group 1"},
                    EvidenceKind.INFERRED,
                    "cohort PTM",
                    sdrf_column="comment[modification parameters]",
                    sdrf_value="NT=Phosphorylation;AC=UniMod:21",
                ),
            ]
        )
        document = SDRFDocument(
            [
                "comment[data file]",
                "comment[modification parameters]",
                "factor value[condition]",
            ],
            [["run.raw", "NT=Oxidation;AC=UniMod:35", "control"]],
        )
        document.annotate(
            [result],
            include_inferred_fields=frozenset({COHORT_MODIFICATION_FIELD}),
            append_columns=frozenset({"comment[modification parameters]"}),
        )
        self.assertEqual(document.columns[-1], "factor value[condition]")
        self.assertEqual(document.rows[0][-1], "control")
        self.assertIn("comment[collision energy]", document.columns[:-1])
        self.assertEqual(document.columns[-2], "comment[modification parameters]")

    def test_sdrf_refinement_replaces_tolerance_but_appends_new_modification(self) -> None:
        result = _result("run.raw")
        result.annotations.extend(
            [
                Annotation(
                    COHORT_PRECURSOR_FIELD,
                    {"experiment_group": "Experiment group 1"},
                    EvidenceKind.INFERRED,
                    "cohort tolerance",
                    sdrf_column="comment[precursor mass tolerance]",
                    sdrf_value="10 ppm",
                ),
                Annotation(
                    COHORT_MODIFICATION_FIELD,
                    {"unimod_accession": "UniMod:21"},
                    EvidenceKind.INFERRED,
                    "cohort PTM",
                    sdrf_column="comment[modification parameters]",
                    sdrf_value="NT=Phosphorylation;AC=UniMod:21",
                ),
            ]
        )
        document = SDRFDocument(
            [
                "comment[data file]",
                "comment[precursor mass tolerance]",
                "comment[modification parameters]",
            ],
            [["run.raw", "20 ppm", "NT=Oxidation;AC=UniMod:35"]],
        )
        changes = document.annotate(
            [result],
            include_inferred_fields=frozenset(
                {COHORT_PRECURSOR_FIELD, COHORT_MODIFICATION_FIELD}
            ),
            overwrite_fields=frozenset({COHORT_PRECURSOR_FIELD}),
            append_columns=frozenset({"comment[modification parameters]"}),
        )
        self.assertEqual(document.rows[0][1], "10 ppm")
        self.assertEqual(document.rows[0][2], "NT=Oxidation;AC=UniMod:35")
        self.assertEqual(document.rows[0][3], "NT=Phosphorylation;AC=UniMod:21")
        self.assertEqual(
            [change.status for change in changes[-2:]], ["replaced", "appended"]
        )

    def test_workflow_writes_original_refined_change_table_and_human_log(self) -> None:
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            source = root / "run.raw"
            source.write_bytes(b"raw")
            sdrf = root / "input.sdrf.tsv"
            sdrf.write_text(
                "comment[data file]\tcomment[precursor mass tolerance]"
                "\tcomment[modification parameters]\n"
                "run.raw\t20 ppm\tNT=Oxidation;AC=UniMod:35\n",
                encoding="utf-8",
            )
            result = _result(
                "run.raw", precursor=9.4, fragment=20.2, mass_shift=79.9663
            )
            validator = MagicMock()
            validator.validate.side_effect = [
                ValidationReport(str(sdrf), "ms-proteomics", False, (), "test"),
                ValidationReport("refined", "ms-proteomics", False, (), "test"),
            ]
            output = root / "qc"
            with patch(
                "prideqc.pipeline._analyze_file",
                return_value=FileOutcome(source, result),
            ):
                manifest = Workflow(
                    WorkflowOptions(refine_sdrf_qc=True), validator=validator
                ).run([source], output, sdrf=sdrf)

            self.assertTrue(manifest["success"])
            self.assertEqual((output / "original.sdrf.tsv").read_bytes(), sdrf.read_bytes())
            refined = SDRFDocument.read(output / "refined.sdrf.tsv")
            self.assertEqual(refined.rows[0][1], "10 ppm")
            self.assertIn("NT=Oxidation;AC=UniMod:35", refined.rows[0])
            self.assertNotIn("NT=Phosphorylation;AC=UniMod:21", [
                refined.rows[0][index]
                for index in refined.indices("comment[modification parameters]")
            ])
            changes = (output / "sdrf-changes.tsv").read_text(encoding="utf-8")
            self.assertIn("annotation_field", changes)
            self.assertIn("cohort_precursor_search_tolerance", changes)
            log = (output / "sdrf-refinement.log.txt").read_text(encoding="utf-8")
            self.assertIn("applied_changes=", log)
            self.assertIn("comment[precursor mass tolerance]", log)
            self.assertIn("sdrf_eligible_ptm_families=0", log)
            self.assertTrue((output / "cohort-refinement.json").exists())
            self.assertTrue((output / "llm-refinement-packet.json").exists())
            self.assertTrue((output / "llm-adjudication-request.json").exists())
            packet = json.loads(
                (output / "llm-refinement-packet.json").read_text(encoding="utf-8")
            )
            self.assertEqual(packet["packet_scope"], "accession-sdrf")
            self.assertEqual(len(packet["runs"]), 1)
            request = json.loads(
                (output / "llm-adjudication-request.json").read_text(encoding="utf-8")
            )
            self.assertEqual(request["request_scope"], "accession-sdrf")
            self.assertEqual(
                manifest["cohort_refinement"]["llm_adjudication_request"],
                "llm-adjudication-request.json",
            )

    def test_existing_results_recovers_accession_and_holds_mass_only_ptm(self) -> None:
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            results_root = root / "results" / "PXD000612"
            results_root.mkdir(parents=True)
            for index in range(10):
                result = _result(
                    f"run-{index}.raw",
                    precursor=7.0 + index * 0.05,
                    fragment=12.0 + index * 0.05,
                    mass_shift=14.0155 + index * 1e-5,
                    candidates=[("UniMod:34", "Methylation")],
                )
                task = results_root / f"task-{index}"
                task.mkdir()
                (task / f"run-{index}.raw.summary.json").write_text(
                    json.dumps(result.to_dict()), encoding="utf-8"
                )
            sdrf = root / "PXD000612.holdout10.sdrf.tsv"
            sdrf.write_text(
                "comment[data file]\tcomment[precursor mass tolerance]"
                "\tfactor value[condition]\n"
                + "".join(
                    f"run-{index}.raw\t20 ppm\tcontrol\n" for index in range(10)
                ),
                encoding="utf-8",
            )
            validator = MagicMock()
            validator.validate.side_effect = [
                ValidationReport(str(sdrf), "ms-proteomics", False, (), "test"),
                ValidationReport("refined", "ms-proteomics", False, (), "test"),
            ]
            output = root / "refinement"
            manifest = Workflow(
                WorkflowOptions(refine_sdrf_qc=True), validator=validator
            ).refine_existing_sdrf(results_root, output, sdrf=sdrf)

            self.assertTrue(manifest["success"])
            self.assertEqual(manifest["project_accession"], "PXD000612")
            self.assertEqual(manifest["sdrf_eligible_ptm_families"], 0)
            self.assertEqual(manifest["ptm_review_families"], 1)
            self.assertEqual(manifest["llm_refinement_packet"], "llm-refinement-packet.json")
            self.assertEqual(
                manifest["llm_adjudication_request"], "llm-adjudication-request.json"
            )
            self.assertTrue((output / "llm-adjudication-request.json").exists())
            packet = json.loads(
                (output / "llm-refinement-packet.json").read_text(encoding="utf-8")
            )
            self.assertEqual(packet["project_accession"], "PXD000612")
            self.assertEqual(len(packet["runs"]), 10)
            refined = SDRFDocument.read(output / "refined.sdrf.tsv")
            self.assertEqual(refined.columns[-1], "factor value[condition]")
            putative_indices = refined.indices(PUTATIVE_MODIFICATION_COLUMN)
            self.assertEqual(len(putative_indices), 1)
            self.assertTrue(
                all(
                    "prideqc putative modification: NT=Methylation;AC=UniMod:34"
                    in row[putative_indices[0]]
                    for row in refined.rows
                )
            )
            standard_modification_indices = refined.indices("comment[modification parameters]")
            self.assertFalse(
                any(
                    "UniMod:34" in row[index]
                    for row in refined.rows
                    for index in standard_modification_indices
                )
            )
            self.assertIn(
                "STATUS=hold-semantic-not-supported", refined.rows[0][putative_indices[0]]
            )
            log = (output / "sdrf-refinement.log.txt").read_text(encoding="utf-8")
            self.assertIn("unimod=UniMod:34", log)
            self.assertIn("sdrf_status=hold-semantic-not-supported", log)


    def test_existing_results_support_no_original_sdrf_adjudication_mode(self) -> None:
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            results_root = root / "results" / "PXD041271"
            results_root.mkdir(parents=True)
            for index in range(4):
                result = _result(
                    f"run-{index}.raw",
                    precursor=8.0 + index * 0.1,
                    fragment=20.0 + index * 0.1,
                    mass_shift=14.0155 + index * 1e-5,
                    candidates=[("UniMod:34", "Methylation")],
                )
                result.project_accession = "PXD041271"
                task = results_root / f"task-{index}"
                task.mkdir()
                (task / f"run-{index}.raw.summary.json").write_text(
                    json.dumps(result.to_dict()), encoding="utf-8"
                )
            validator = MagicMock()
            output = root / "manifest-only"
            manifest = Workflow(
                WorkflowOptions(refine_sdrf_qc=True), validator=validator
            ).refine_existing_sdrf(results_root, output, sdrf=None)

            self.assertTrue(manifest["success"])
            self.assertEqual(manifest["input_mode"], "no-original-sdrf")
            self.assertEqual(
                manifest["mode"], "existing-results-no-original-sdrf-adjudication"
            )
            self.assertEqual(manifest["summary_count"], 4)
            self.assertEqual(manifest["project_accession"], "PXD041271")
            self.assertFalse(manifest["sdrf_writeback_supported"])
            self.assertIsNone(manifest["original_sdrf"])
            self.assertIsNone(manifest["refined_sdrf"])
            self.assertFalse((output / "original.sdrf.tsv").exists())
            self.assertFalse((output / "refined.sdrf.tsv").exists())
            self.assertFalse((output / "sdrf-validation.json").exists())
            self.assertTrue((output / "cohort-refinement.json").exists())
            self.assertTrue((output / "llm-refinement-packet.json").exists())
            self.assertTrue((output / "llm-adjudication-request.json").exists())
            packet = json.loads(
                (output / "llm-refinement-packet.json").read_text(encoding="utf-8")
            )
            request = json.loads(
                (output / "llm-adjudication-request.json").read_text(encoding="utf-8")
            )
            self.assertEqual(packet["provenance"]["input_mode"], "no-original-sdrf")
            self.assertEqual(request["input_mode"], "no-original-sdrf")
            self.assertFalse(packet["sdrf"]["available"])
            self.assertTrue(request["decisions"])
            self.assertTrue(all(item["target_rows"] == [] for item in request["decisions"]))
            validator.validate.assert_not_called()

    def test_no_original_sdrf_requires_unambiguous_project_accession(self) -> None:
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            results_root = root / "results"
            results_root.mkdir()
            result = _result("run.raw", precursor=8.0, fragment=20.0)
            (results_root / "run.raw.summary.json").write_text(
                json.dumps(result.to_dict()), encoding="utf-8"
            )
            with self.assertRaisesRegex(ValueError, "requires one unambiguous"):
                Workflow(WorkflowOptions(refine_sdrf_qc=True)).refine_existing_sdrf(
                    results_root, root / "out", sdrf=None
                )

    def test_existing_results_refinement_combines_file_array_summaries(self) -> None:
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            results_root = root / "array-results"
            results_root.mkdir()
            for index in range(10):
                result = _result(
                    f"run-{index}.raw",
                    precursor=8.0 + index * 0.1,
                    fragment=20.0 + index * 0.1,
                    mass_shift=79.9663 + index * 1e-5,
                )
                result.project_accession = "PXD041271"
                task = results_root / f"task-{index}"
                task.mkdir()
                (task / f"run-{index}.raw.summary.json").write_text(
                    json.dumps(result.to_dict()), encoding="utf-8"
                )
            sdrf = root / "input.sdrf.tsv"
            sdrf.write_text(
                "comment[data file]\tcomment[precursor mass tolerance]\n"
                + "".join(f"run-{index}.raw\t20 ppm\n" for index in range(10)),
                encoding="utf-8",
            )
            validator = MagicMock()
            validator.validate.side_effect = [
                ValidationReport(str(sdrf), "ms-proteomics", False, (), "test"),
                ValidationReport("refined", "ms-proteomics", False, (), "test"),
            ]
            evidence = root / "study-evidence.tsv"
            evidence.write_text(
                "pxd_accession\tunimod_accession\tevidence_status\n"
                "PXD041271\tUniMod:21\tsupported\n",
                encoding="utf-8",
            )
            output = results_root / "sdrf-refinement"
            manifest = Workflow(
                WorkflowOptions(
                    refine_sdrf_qc=True, ptm_study_evidence=str(evidence)
                ),
                validator=validator,
            ).refine_existing_sdrf(results_root, output, sdrf=sdrf)

            self.assertTrue(manifest["success"])
            self.assertEqual(manifest["summary_count"], 10)
            self.assertEqual(manifest["project_accession"], "PXD041271")
            self.assertEqual(manifest["sdrf_eligible_ptm_families"], 1)
            refined = SDRFDocument.read(output / "refined.sdrf.tsv")
            self.assertTrue(all(row[1] == "9 ppm" for row in refined.rows))
            self.assertTrue(
                all("NT=Phosphorylation;AC=UniMod:21" in row for row in refined.rows)
            )


if __name__ == "__main__":
    unittest.main()
