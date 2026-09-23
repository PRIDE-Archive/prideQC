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

    def test_only_strict_unambiguous_recurrent_ptm_family_is_sdrf_eligible(self) -> None:
        results = [
            _result(f"run-{index}.raw", mass_shift=79.9663 + index * 1e-5)
            for index in range(10)
        ]
        synthesis = synthesize_cohort(results)
        self.assertEqual(len(synthesis.ptm_families), 1)
        self.assertEqual(synthesis.ptm_families[0]["unimod_accession"], "UniMod:21")
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
            self.assertNotIn("NT=Phosphorylation;AC=UniMod:21", refined.rows[0])
            changes = (output / "sdrf-changes.tsv").read_text(encoding="utf-8")
            self.assertIn("annotation_field", changes)
            self.assertIn("cohort_precursor_search_tolerance", changes)
            log = (output / "sdrf-refinement.log.txt").read_text(encoding="utf-8")
            self.assertIn("applied_changes=", log)
            self.assertIn("comment[precursor mass tolerance]", log)
            self.assertIn("sdrf_eligible_ptm_families=0", log)
            self.assertTrue((output / "cohort-refinement.json").exists())

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
            output = results_root / "sdrf-refinement"
            manifest = Workflow(
                WorkflowOptions(refine_sdrf_qc=True), validator=validator
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
