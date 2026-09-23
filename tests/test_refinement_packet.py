"""Deterministic accession-level evidence packet contracts."""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import jsonschema

from prideqc import __version__
from prideqc.cohort import synthesize_cohort
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

    def test_packet_contains_ptm_candidate_evidence_and_probability_guard(self) -> None:
        with tempfile.TemporaryDirectory() as folder:
            packet, _, _, _ = self._fixture(Path(folder))

        ptms = [
            item
            for item in packet["decision_candidates"]
            if item["decision_type"] == "modification"
        ]
        self.assertEqual(len(ptms), 2)
        for item in ptms:
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


    def test_packet_retains_mass_ambiguous_candidate_set_for_adjudication(self) -> None:
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

        ptms = [
            item
            for item in packet["decision_candidates"]
            if item["decision_type"] == "modification"
        ]
        self.assertEqual(len(ptms), 1)
        self.assertEqual(len(ptms[0]["candidate_values"]), 2)
        self.assertTrue(ptms[0]["evidence"]["mass_identity_ambiguous"])
        self.assertEqual(
            {option["accession"] for option in ptms[0]["evidence"]["candidate_options"]},
            {"UniMod:21", "UniMod:99913"},
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
