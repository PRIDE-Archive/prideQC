from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from prideqc.refinement_adjudication import (
    REQUEST_SCHEMA_VERSION,
    REQUEST_SCOPE,
    RESPONSE_SCHEMA_VERSION,
)
from prideqc.refinement_packet import PACKET_SCHEMA_VERSION
from prideqc.sdrf import SDRFDocument
from prideqc.submission import apply_adjudicated_sdrf
from prideqc.validation import ValidationReport


class _Validator:
    def __init__(self, issues: list[tuple[str, ...]]) -> None:
        self.issues = list(issues)
        self.paths: list[Path] = []

    def validate(self, path: Path) -> ValidationReport:
        self.paths.append(path)
        findings = self.issues.pop(0) if self.issues else ()
        return ValidationReport(
            str(path),
            "ms-proteomics",
            False,
            findings,
            "test",
        )


def _request(decisions: list[dict[str, object]], *, mode: str = "sdrf-backed") -> dict[str, object]:
    source_hash = "a" * 64
    return {
        "schema_version": REQUEST_SCHEMA_VERSION,
        "request_scope": REQUEST_SCOPE,
        "request_id": f"sha256:{source_hash}",
        "project_accession": "PXD123456",
        "input_mode": mode,
        "source_packet": {
            "schema_version": PACKET_SCHEMA_VERSION,
            "sha256": source_hash,
            "prideqc_version": "0.2.0",
        },
        "sdrf": {
            "source_name": "PXD123456.benchmark.sdrf.tsv",
            "sha256": "b" * 64,
        },
        "contract": {},
        "decisions": decisions,
    }


def _decision(
    decision_id: str,
    *,
    target_field: str,
    target_run: str,
    data_file: str,
    values: list[str],
    candidate: str,
    write_semantics: str = "fill_or_replace_canonical_value",
    decision_type: str = "mass_tolerance",
    column_count: int = 1,
) -> dict[str, object]:
    return {
        "decision_id": decision_id,
        "decision_type": decision_type,
        "experiment_group": "Experiment group 1",
        "target_field": target_field,
        "write_semantics": write_semantics,
        # Deliberately wrong/benchmark-local. Full application must ignore this.
        "target_rows": [999],
        "target_runs": [target_run],
        "original": {
            "column_present": column_count > 0,
            "column_count": column_count,
            "rows": [
                {
                    "row": 999,
                    "data_file": data_file,
                    "values": values,
                }
            ],
        },
        "candidate_values": [candidate],
        "evidence": {},
        "evidence_semantics": {},
        "local_context": {"ptm_families": []},
    }


def _response(request: dict[str, object], rows: list[dict[str, object]]) -> dict[str, object]:
    return {
        "schema_version": RESPONSE_SCHEMA_VERSION,
        "request_id": request["request_id"],
        "project_accession": request["project_accession"],
        "input_mode": request["input_mode"],
        "decisions": rows,
    }


def _write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")


class SubmissionApplicationTests(unittest.TestCase):
    def test_applies_only_accepted_candidates_to_full_sdrf_and_stamps_provenance(self) -> None:
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            sdrf = root / "PXD123456.sdrf.tsv"
            sdrf.write_text(
                "source name\tassay name\ttechnology type\tcomment[data file]\t"
                "comment[precursor mass tolerance]\tcomment[fragment mass tolerance]\t"
                "comment[modification parameters]\tfactor value[condition]\n"
                "s1\ta1\tproteomic profiling by mass spectrometry\trun-a.raw\t20 ppm\t0.02 Da\t"
                "NT=Oxidation;AC=UniMod:35\tA\n"
                "s2\ta2\tproteomic profiling by mass spectrometry\trun-b.raw\t20 ppm\t0.02 Da\t"
                "NT=Oxidation;AC=UniMod:35\tB\n"
                "s3\ta3\tproteomic profiling by mass spectrometry\trun-c.raw\t20 ppm\t0.02 Da\t"
                "NT=Oxidation;AC=UniMod:35\tC\n",
                encoding="utf-8",
            )
            request = _request(
                [
                    _decision(
                        "Experiment group 1:precursor-mass-tolerance",
                        target_field="comment[precursor mass tolerance]",
                        target_run="run-a.raw",
                        data_file="run-a.raw",
                        values=["20 ppm"],
                        candidate="10 ppm",
                    ),
                    _decision(
                        "Experiment group 1:fragment-mass-tolerance",
                        target_field="comment[fragment mass tolerance]",
                        target_run="run-a.raw",
                        data_file="run-a.raw",
                        values=["0.02 Da"],
                        candidate="11 ppm",
                    ),
                    _decision(
                        "Experiment group 1:modification-family:14.015500",
                        target_field="comment[modification parameters]",
                        target_run="run-b.raw",
                        data_file="run-b.raw",
                        values=["NT=Oxidation;AC=UniMod:35"],
                        candidate="NT=Methylation;AC=UniMod:34",
                        write_semantics=(
                            "append_one_allowed_canonical_value_if_accepted_and_missing"
                        ),
                        decision_type="modification",
                    ),
                ]
            )
            response = _response(
                request,
                [
                    {
                        "decision_id": "Experiment group 1:fragment-mass-tolerance",
                        "decision": "reject",
                        "selected_value": None,
                        "reason": "Preserve the original.",
                    },
                    {
                        "decision_id": "Experiment group 1:modification-family:14.015500",
                        "decision": "accept",
                        "selected_value": "NT=Methylation;AC=UniMod:34",
                        "reason": "Use the supplied candidate.",
                    },
                    {
                        "decision_id": "Experiment group 1:precursor-mass-tolerance",
                        "decision": "accept",
                        "selected_value": "10 ppm",
                        "reason": "Use the supplied candidate.",
                    },
                ],
            )
            request_path = root / "request.json"
            decisions_path = root / "decisions.json"
            _write_json(request_path, request)
            _write_json(decisions_path, response)

            output = root / "out"
            manifest = apply_adjudicated_sdrf(
                sdrf=sdrf,
                request_path=request_path,
                decisions_path=decisions_path,
                output_directory=output,
                validator=_Validator([(), (), ()]),
            )

            self.assertTrue(manifest["submission_ready"])
            self.assertEqual(manifest["changed_cell_count"], 2)
            self.assertEqual(manifest["annotation_tool"], "prideQC v0.2.0")
            final = SDRFDocument.read(output / "PXD123456.sdrf.tsv")
            precursor = final.indices("comment[precursor mass tolerance]")[0]
            fragment = final.indices("comment[fragment mass tolerance]")[0]
            mods = final.indices("comment[modification parameters]")
            tool = final.indices("comment[sdrf annotation tool]")[0]
            self.assertEqual(final.rows[0][precursor], "10 ppm")
            self.assertEqual(final.rows[0][fragment], "0.02 Da")
            self.assertEqual(final.rows[2][precursor], "20 ppm")
            self.assertIn("NT=Methylation;AC=UniMod:34", [final.rows[1][i] for i in mods])
            self.assertNotIn("NT=Methylation;AC=UniMod:34", [final.rows[2][i] for i in mods])
            self.assertTrue(all(row[tool] == "prideQC v0.2.0" for row in final.rows))

    def test_refuses_full_sdrf_original_value_drift(self) -> None:
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            sdrf = root / "PXD123456.sdrf.tsv"
            sdrf.write_text(
                "comment[data file]\tcomment[precursor mass tolerance]\n"
                "run.raw\t30 ppm\n",
                encoding="utf-8",
            )
            request = _request(
                [
                    _decision(
                        "Experiment group 1:precursor-mass-tolerance",
                        target_field="comment[precursor mass tolerance]",
                        target_run="run.raw",
                        data_file="run.raw",
                        values=["20 ppm"],
                        candidate="10 ppm",
                    )
                ]
            )
            response = _response(
                request,
                [
                    {
                        "decision_id": "Experiment group 1:precursor-mass-tolerance",
                        "decision": "accept",
                        "selected_value": "10 ppm",
                        "reason": "Use candidate.",
                    }
                ],
            )
            request_path = root / "request.json"
            decisions_path = root / "decisions.json"
            _write_json(request_path, request)
            _write_json(decisions_path, response)

            with self.assertRaisesRegex(ValueError, "original SDRF values drifted"):
                apply_adjudicated_sdrf(
                    sdrf=sdrf,
                    request_path=request_path,
                    decisions_path=decisions_path,
                    output_directory=root / "out",
                    validator=_Validator([()]),
                )

    def test_exact_file_map_rebinds_converted_sdrf_name_to_raw_target(self) -> None:
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            sdrf = root / "PXD123456.sdrf.tsv"
            sdrf.write_text(
                "comment[data file]\tcomment[precursor mass tolerance]\n"
                "run.mzML\t20 ppm\n",
                encoding="utf-8",
            )
            request = _request(
                [
                    _decision(
                        "Experiment group 1:precursor-mass-tolerance",
                        target_field="comment[precursor mass tolerance]",
                        target_run="run.raw",
                        data_file="run.mzML",
                        values=["20 ppm"],
                        candidate="10 ppm",
                    )
                ]
            )
            response = _response(
                request,
                [
                    {
                        "decision_id": "Experiment group 1:precursor-mass-tolerance",
                        "decision": "accept",
                        "selected_value": "10 ppm",
                        "reason": "Use candidate.",
                    }
                ],
            )
            request_path = root / "request.json"
            decisions_path = root / "decisions.json"
            mapping = root / "file-map.json"
            _write_json(request_path, request)
            _write_json(decisions_path, response)
            _write_json(mapping, {"run.mzML": "run.raw"})

            manifest = apply_adjudicated_sdrf(
                sdrf=sdrf,
                request_path=request_path,
                decisions_path=decisions_path,
                output_directory=root / "out",
                file_map=mapping,
                validator=_Validator([(), (), ()]),
            )
            self.assertTrue(manifest["submission_ready"])

    def test_invalid_candidate_output_is_not_stamped_or_submission_ready(self) -> None:
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            sdrf = root / "PXD123456.sdrf.tsv"
            sdrf.write_text(
                "comment[data file]\tcomment[precursor mass tolerance]\n"
                "run.raw\t20 ppm\n",
                encoding="utf-8",
            )
            request = _request(
                [
                    _decision(
                        "Experiment group 1:precursor-mass-tolerance",
                        target_field="comment[precursor mass tolerance]",
                        target_run="run.raw",
                        data_file="run.raw",
                        values=["20 ppm"],
                        candidate="10 ppm",
                    )
                ]
            )
            response = _response(
                request,
                [
                    {
                        "decision_id": "Experiment group 1:precursor-mass-tolerance",
                        "decision": "accept",
                        "selected_value": "10 ppm",
                        "reason": "Use candidate.",
                    }
                ],
            )
            request_path = root / "request.json"
            decisions_path = root / "decisions.json"
            _write_json(request_path, request)
            _write_json(decisions_path, response)

            output = root / "out"
            manifest = apply_adjudicated_sdrf(
                sdrf=sdrf,
                request_path=request_path,
                decisions_path=decisions_path,
                output_directory=output,
                validator=_Validator([(), ("invalid candidate",)]),
            )
            self.assertFalse(manifest["submission_ready"])
            self.assertIsNone(manifest["annotation_tool"])
            self.assertIsNone(manifest["submission_sdrf"])
            candidate = SDRFDocument.read(output / "candidate.sdrf.tsv")
            self.assertEqual(candidate.indices("comment[sdrf annotation tool]"), [])

    def test_explicit_technology_type_normalization_can_make_structural_update(self) -> None:
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            sdrf = root / "PXD123456.sdrf.tsv"
            sdrf.write_text(
                "comment[data file]\ttechnology type\tassay name\n"
                "run.raw\tproteomic profiling by mass spectrometry\ta1\n",
                encoding="utf-8",
            )
            request = _request([])
            response = _response(request, [])
            request_path = root / "request.json"
            decisions_path = root / "decisions.json"
            _write_json(request_path, request)
            _write_json(decisions_path, response)

            output = root / "out"
            manifest = apply_adjudicated_sdrf(
                sdrf=sdrf,
                request_path=request_path,
                decisions_path=decisions_path,
                output_directory=output,
                normalize_technology_type_order=True,
                validator=_Validator(
                    [
                        ("technology type must follow assay name",),
                        (),
                        (),
                    ]
                ),
            )
            self.assertTrue(manifest["submission_ready"])
            self.assertEqual(len(manifest["structural_changes"]), 1)
            self.assertEqual(manifest["annotation_tool"], "prideQC v0.2.0")
            final = SDRFDocument.read(output / "PXD123456.sdrf.tsv")
            assay = final.indices("assay name")[0]
            technology = final.indices("technology type")[0]
            self.assertEqual(technology, assay + 1)

    def test_no_original_sdrf_adjudication_cannot_be_applied(self) -> None:
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            sdrf = root / "PXD123456.sdrf.tsv"
            sdrf.write_text("comment[data file]\nrun.raw\n", encoding="utf-8")
            request = _request([], mode="no-original-sdrf")
            response = _response(request, [])
            request_path = root / "request.json"
            decisions_path = root / "decisions.json"
            _write_json(request_path, request)
            _write_json(decisions_path, response)

            with self.assertRaisesRegex(ValueError, "requires an sdrf-backed"):
                apply_adjudicated_sdrf(
                    sdrf=sdrf,
                    request_path=request_path,
                    decisions_path=decisions_path,
                    output_directory=root / "out",
                    validator=_Validator([()]),
                )


if __name__ == "__main__":
    unittest.main()
