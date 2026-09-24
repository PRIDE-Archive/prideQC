"""Deterministic accession-level evidence packet contracts."""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import jsonschema

from prideqc import __version__
from prideqc.cohort import synthesize_cohort
from prideqc.refinement_adjudication import (
    build_llm_adjudication_request,
    validate_llm_adjudication_request,
)
from prideqc.refinement_packet import (
    PACKET_SCHEMA_VERSION,
    build_llm_refinement_packet,
    validate_llm_refinement_packet,
)
from prideqc.sdrf import SDRFDocument
from tests.test_cohort_refinement import _result


class RefinementPacketTests(unittest.TestCase):
    def _fixture(self, root: Path):
        results = [
            *[
                _result(
                    f"short-{index}.raw",
                    precursor=7.0 + index * 0.05,
                    fragment=10.0 + index * 0.05,
                    duration=7200 + index,
                    mass_shift=14.0155 + index * 1e-5,
                    candidates=[("UniMod:34", "Methylation")],
                )
                for index in range(4)
            ],
            *[
                _result(
                    f"long-{index}.raw",
                    precursor=7.4 + index * 0.05,
                    fragment=12.0 + index * 0.05,
                    duration=15900 + index,
                    mass_shift=14.0155 + index * 1e-5,
                    candidates=[("UniMod:34", "Methylation")],
                )
                for index in range(6)
            ],
        ]
        for result in results:
            result.project_accession = "PXD000612"
        lines = [
            "comment[data file]\tcomment[precursor mass tolerance]"
            "\tcomment[fragment mass tolerance]\tcomment[modification parameters]"
            "\tcomment[modification parameters]\tfactor value[condition]\n"
        ]
        for result in results:
            lines.append(
                f"{result.input_path.name}\tnot available\tnot available"
                "\tNT=Oxidation;AC=UniMod:35\tNT=Carbamidomethyl;AC=UniMod:4\tcontrol\n"
            )
        sdrf_path = root / "PXD000612.sdrf.tsv"
        sdrf_path.write_text("".join(lines), encoding="utf-8")
        document = SDRFDocument.read(sdrf_path)
        synthesis = synthesize_cohort(results, project_accession="PXD000612")
        artifacts = {
            result.input_path.name: {
                "summary_json": f"task/{result.input_path.name}.summary.json",
                "mzqc": f"task/{result.input_path.name}.mzQC",
            }
            for result in results
        }
        packet = build_llm_refinement_packet(
            document,
            results,
            synthesis,
            project_accession="PXD000612",
            sdrf_path=sdrf_path,
            prideqc_version=__version__,
            source_artifacts=artifacts,
        )
        return packet, results, synthesis, sdrf_path

    def test_packet_is_one_accession_scope_with_all_runs_groups_and_tolerances(self) -> None:
        with tempfile.TemporaryDirectory() as folder:
            packet, results, synthesis, _ = self._fixture(Path(folder))

        validate_llm_refinement_packet(packet)
        self.assertEqual(packet["schema_version"], PACKET_SCHEMA_VERSION)
        self.assertEqual(packet["packet_scope"], "accession-sdrf")
        self.assertEqual(packet["project_accession"], "PXD000612")
        self.assertEqual(len(packet["runs"]), len(results))
        self.assertEqual(len(packet["experiment_groups"]), len(synthesis.groups))
        self.assertEqual(len(packet["experiment_groups"]), 2)

        tolerance_decisions = [
            item
            for item in packet["decision_candidates"]
            if item["decision_type"] == "mass_tolerance"
        ]
        self.assertEqual(len(tolerance_decisions), 4)
        fields = {item["target_field"] for item in tolerance_decisions}
        self.assertEqual(
            fields,
            {
                "comment[precursor mass tolerance]",
                "comment[fragment mass tolerance]",
            },
        )
        self.assertTrue(
            all(item["evidence"]["per_run_estimates"] for item in tolerance_decisions)
        )
        self.assertTrue(
            all(
                run["tolerance_evidence"]["precursor"]["status"] == "available"
                and run["tolerance_evidence"]["fragment"]["status"] == "available"
                for run in packet["runs"]
            )
        )
        self.assertEqual(len(packet["ptm_context"]), 2)
        self.assertTrue(all(not family["actionable"] for family in packet["ptm_context"]))
        self.assertTrue(
            all(
                group["actionable_ptm_review_family_count"] == 1
                for group in packet["experiment_groups"]
            )
        )

    def test_packet_contains_ptm_candidate_evidence_and_probability_guard(self) -> None:
        with tempfile.TemporaryDirectory() as folder:
            packet, _, _, _ = self._fixture(Path(folder))

        ptms = [
            item
            for item in packet["decision_candidates"]
            if item["decision_type"] == "modification"
        ]
        self.assertEqual(len(ptms), 2)
        self.assertEqual(len(packet["ptm_context"]), 2)
        for item in ptms:
            self.assertEqual(item["evidence"]["review_gate"], "strict-cohort-ptm-review-family")
            self.assertEqual(
                item["candidate_values"], ["NT=Methylation;AC=UniMod:34"]
            )
            self.assertEqual(
                item["evidence"]["candidate_options"][0]["accession"], "UniMod:34"
            )
            self.assertGreaterEqual(item["evidence"]["supporting_runs"], 3)
            self.assertTrue(item["evidence"]["per_run_support"])
            self.assertFalse(
                item["evidence_semantics"][
                    "recurrent_family_probability_is_identity_probability"
                ]
            )


    def test_packet_keeps_mass_ambiguous_family_as_non_actionable_context(self) -> None:
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            results = [
                _result(
                    f"run-{index}.raw",
                    mass_shift=79.9663 + index * 1e-5,
                    candidates=[
                        ("UniMod:21", "Phosphorylation"),
                        ("UniMod:99913", "Alternative phosphate-like candidate"),
                    ],
                )
                for index in range(10)
            ]
            sdrf_path = root / "PXDTEST.sdrf.tsv"
            sdrf_path.write_text(
                "comment[data file]\tcomment[modification parameters]\n"
                + "".join(
                    f"run-{index}.raw\tNT=Oxidation;AC=UniMod:35\n"
                    for index in range(10)
                ),
                encoding="utf-8",
            )
            document = SDRFDocument.read(sdrf_path)
            synthesis = synthesize_cohort(results, project_accession="PXDTEST")
            self.assertEqual(synthesis.ptm_review_families, [])
            packet = build_llm_refinement_packet(
                document,
                results,
                synthesis,
                project_accession="PXDTEST",
                sdrf_path=sdrf_path,
                prideqc_version=__version__,
            )

        ptm_decisions = [
            item
            for item in packet["decision_candidates"]
            if item["decision_type"] == "modification"
        ]
        self.assertEqual(ptm_decisions, [])
        self.assertEqual(len(packet["ptm_context"]), 1)
        context = packet["ptm_context"][0]
        self.assertFalse(context["actionable"])
        self.assertTrue(context["mass_identity_ambiguous"])
        self.assertEqual(len(context["candidate_values"]), 2)
        self.assertEqual(
            {option["accession"] for option in context["candidate_options"]},
            {"UniMod:21", "UniMod:99913"},
        )

    def test_broad_ptm_context_does_not_expand_actionable_decisions(self) -> None:
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            results = [
                _result(
                    f"run-{index}.raw",
                    mass_shift=14.0155 + index * 1e-5,
                    candidates=[("UniMod:34", "Methylation")],
                )
                for index in range(10)
            ]
            extra_families = [
                (79.9663, [("UniMod:21", "Phosphorylation"), ("UniMod:99913", "Alt")]),
                (
                    57.0215,
                    [
                        ("UniMod:4", "Iodoacetamide derivative"),
                        ("UniMod:1263", "Addition of Glycine"),
                    ],
                ),
                (42.0103, [("UniMod:1", "Acetylation"), ("UniMod:52", "Guanidination")]),
            ]
            for index, result in enumerate(results):
                annotation = next(
                    item
                    for item in result.annotations
                    if item.field == "putative_modification_mass_shifts"
                )
                assert isinstance(annotation.value, list)
                for mass, candidates in extra_families:
                    rows = [
                        {
                            "unimod_accession": accession,
                            "name": name,
                            "candidate_category": "biological-ptm",
                        }
                        for accession, name in candidates
                    ]
                    annotation.value.append(
                        {
                            "delta_mass_da": mass + index * 1e-5,
                            "pair_support": 120,
                            "unique_spectrum_support": 80,
                            "classification": "putative-ptm",
                            "confidence": "high-support",
                            "unimod_candidates": rows,
                            "diagnostic_unimod_candidates": rows,
                        }
                    )

            sdrf_path = root / "PXDTEST.sdrf.tsv"
            sdrf_path.write_text(
                "comment[data file]\tcomment[modification parameters]\n"
                + "".join(
                    f"run-{index}.raw\tNT=Oxidation;AC=UniMod:35\n"
                    for index in range(10)
                ),
                encoding="utf-8",
            )
            document = SDRFDocument.read(sdrf_path)
            synthesis = synthesize_cohort(results, project_accession="PXDTEST")
            self.assertEqual(len(synthesis.ptm_review_families), 1)
            packet = build_llm_refinement_packet(
                document,
                results,
                synthesis,
                project_accession="PXDTEST",
                sdrf_path=sdrf_path,
                prideqc_version=__version__,
            )

        ptm_decisions = [
            item
            for item in packet["decision_candidates"]
            if item["decision_type"] == "modification"
        ]
        self.assertEqual(len(packet["ptm_context"]), 4)
        self.assertEqual(len(ptm_decisions), 1)
        self.assertEqual(
            ptm_decisions[0]["candidate_values"],
            ["NT=Methylation;AC=UniMod:34"],
        )

    def test_packet_preserves_repeated_original_modification_context(self) -> None:
        with tempfile.TemporaryDirectory() as folder:
            packet, _, _, _ = self._fixture(Path(folder))

        first = packet["runs"][0]["original_values"]["comment[modification parameters]"]
        self.assertEqual(first["column_count"], 2)
        self.assertEqual(
            first["rows"][0]["values"],
            [
                "NT=Oxidation;AC=UniMod:35",
                "NT=Carbamidomethyl;AC=UniMod:4",
            ],
        )


    def test_packet_supports_no_original_sdrf_mode(self) -> None:
        with tempfile.TemporaryDirectory():
            results = [
                _result(
                    f"run-{index}.raw",
                    mass_shift=14.0155 + index * 1e-5,
                    candidates=[("UniMod:34", "Methylation")],
                )
                for index in range(4)
            ]
            synthesis = synthesize_cohort(results, project_accession="PXDTEST")
            packet = build_llm_refinement_packet(
                None,
                results,
                synthesis,
                project_accession="PXDTEST",
                sdrf_path=None,
                prideqc_version=__version__,
            )
            validate_llm_refinement_packet(packet)
            packet_schema = json.loads(
                Path("src/prideqc/data/llm-refinement-packet-v1.schema.json").read_text(
                    encoding="utf-8"
                )
            )
            jsonschema.validate(packet, packet_schema)
            request = build_llm_adjudication_request(packet)
            validate_llm_adjudication_request(request)
            request_schema = json.loads(
                Path("src/prideqc/data/llm-adjudication-request-v1.schema.json").read_text(
                    encoding="utf-8"
                )
            )
            jsonschema.validate(request, request_schema)

        self.assertFalse(packet["sdrf"]["available"])
        self.assertIsNone(packet["sdrf"]["source_name"])
        self.assertIsNone(packet["sdrf"]["sha256"])
        self.assertEqual(packet["sdrf"]["row_count"], 0)
        self.assertFalse(packet["limitations"]["original_sdrf_available"])
        self.assertFalse(packet["limitations"]["sdrf_writeback_supported"])
        self.assertTrue(all(run["sdrf_rows"] == [] for run in packet["runs"]))
        self.assertTrue(all(run["sdrf_data_files"] == [] for run in packet["runs"]))
        self.assertTrue(
            all(
                context["status"] == "unavailable"
                for run in packet["runs"]
                for context in run["original_values"].values()
            )
        )
        self.assertTrue(all(item["target_rows"] == [] for item in request["decisions"]))
        self.assertTrue(
            all(
                item["write_semantics"] == "evidence_only_no_original_sdrf"
                for item in request["decisions"]
            )
        )

    def test_packet_validates_against_packaged_json_schema_and_is_deterministic(self) -> None:
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            packet, results, synthesis, sdrf_path = self._fixture(root)
            document = SDRFDocument.read(sdrf_path)
            again = build_llm_refinement_packet(
                document,
                results,
                synthesis,
                project_accession="PXD000612",
                sdrf_path=sdrf_path,
                prideqc_version=__version__,
                source_artifacts={
                    result.input_path.name: {
                        "summary_json": f"task/{result.input_path.name}.summary.json",
                        "mzqc": f"task/{result.input_path.name}.mzQC",
                    }
                    for result in results
                },
            )

        schema = json.loads(
            Path("src/prideqc/data/llm-refinement-packet-v1.schema.json").read_text(
                encoding="utf-8"
            )
        )
        jsonschema.validate(packet, schema)
        self.assertEqual(packet, again)
        self.assertEqual(
            json.dumps(packet, indent=2, ensure_ascii=False),
            json.dumps(again, indent=2, ensure_ascii=False),
        )


if __name__ == "__main__":
    unittest.main()
