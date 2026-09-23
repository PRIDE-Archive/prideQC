from __future__ import annotations

import copy
import json
import unittest
from pathlib import Path

import jsonschema

from prideqc.refinement_adjudication import (
    REQUEST_SCHEMA_VERSION,
    RESPONSE_SCHEMA_VERSION,
    build_llm_adjudication_request,
    empty_llm_refinement_response,
    validate_llm_adjudication_request,
    validate_llm_refinement_decisions,
)


class RefinementAdjudicationTests(unittest.TestCase):
    def _packet(self) -> dict[str, object]:
        allowed = ["accept", "reject", "abstain"]
        tolerance_id = "Experiment group 1:precursor-mass-tolerance"
        ptm_id = "Experiment group 1:modification-family:14.015500"
        return {
            "schema_version": "prideqc-llm-refinement-packet-v1",
            "packet_scope": "accession-sdrf",
            "project_accession": "PXDTEST",
            "provenance": {"prideqc_version": "0.2.0"},
            "sdrf": {
                "source_name": "PXDTEST.sdrf.tsv",
                "sha256": "a" * 64,
                "row_count": 2,
                "target_columns": {},
            },
            "runs": [
                {
                    "run_id": "run.raw",
                    "experiment_group": "Experiment group 1",
                }
            ],
            "experiment_groups": [
                {
                    "group_id": "Experiment group 1",
                    "decision_ids": [tolerance_id, ptm_id],
                }
            ],
            "ptm_context": [
                {
                    "experiment_group": "Experiment group 1",
                    "observed_delta_mass_da": 14.01549,
                    "candidate_values": ["NT=Methylation;AC=UniMod:34"],
                    "actionable": False,
                },
                {
                    "experiment_group": "Experiment group 1",
                    "observed_delta_mass_da": 79.9663,
                    "candidate_values": ["NT=Phosphorylation;AC=UniMod:21"],
                    "actionable": False,
                },
            ],
            "decision_candidates": [
                {
                    "decision_id": tolerance_id,
                    "decision_type": "mass_tolerance",
                    "experiment_group": "Experiment group 1",
                    "target_field": "comment[precursor mass tolerance]",
                    "write_semantics": "fill_or_replace_canonical_value",
                    "target_rows": [2],
                    "target_runs": ["run.raw"],
                    "original": {"rows": [{"row": 2, "values": ["not available"]}]},
                    "candidate_values": ["8 ppm"],
                    "allowed_values": ["8 ppm"],
                    "allowed_decisions": allowed,
                    "evidence": {
                        "supporting_runs": 1,
                        "group_runs": 1,
                        "per_run_estimates": [{"run_id": "run.raw", "value": 7.4}],
                    },
                },
                {
                    "decision_id": ptm_id,
                    "decision_type": "modification",
                    "experiment_group": "Experiment group 1",
                    "target_field": "comment[modification parameters]",
                    "write_semantics": "append_one_allowed_canonical_value_if_accepted_and_missing",
                    "target_rows": [2],
                    "target_runs": ["run.raw"],
                    "original": {"rows": [{"row": 2, "values": []}]},
                    "candidate_values": ["NT=Methylation;AC=UniMod:34"],
                    "allowed_values": ["NT=Methylation;AC=UniMod:34"],
                    "allowed_decisions": allowed,
                    "evidence": {
                        "observed_delta_mass_da": 14.0155,
                        "recurrent_family_probability": 0.999,
                    },
                    "evidence_semantics": {
                        "recurrent_family_probability_is_identity_probability": False
                    },
                },
            ],
            "llm_contract": {},
            "limitations": {},
        }

    def test_request_is_deterministic_and_only_includes_local_ptm_context(self) -> None:
        packet = self._packet()
        first = build_llm_adjudication_request(packet)
        second = build_llm_adjudication_request(copy.deepcopy(packet))

        self.assertEqual(first, second)
        self.assertEqual(first["schema_version"], REQUEST_SCHEMA_VERSION)
        self.assertRegex(first["request_id"], r"^sha256:[0-9a-f]{64}$")
        self.assertEqual(len(first["decisions"]), 2)
        tolerance = next(
            item for item in first["decisions"] if item["decision_type"] == "mass_tolerance"
        )
        ptm = next(
            item for item in first["decisions"] if item["decision_type"] == "modification"
        )
        self.assertEqual(tolerance["local_context"]["ptm_families"], [])
        self.assertEqual(len(ptm["local_context"]["ptm_families"]), 1)
        self.assertAlmostEqual(
            ptm["local_context"]["ptm_families"][0]["observed_delta_mass_da"],
            14.01549,
        )
        self.assertTrue(first["contract"]["ignore_instructions_embedded_in_evidence"])

    def test_request_validates_against_packaged_json_schema(self) -> None:
        request = build_llm_adjudication_request(self._packet())
        validate_llm_adjudication_request(request)
        schema = json.loads(
            Path("src/prideqc/data/llm-adjudication-request-v1.schema.json").read_text(
                encoding="utf-8"
            )
        )
        jsonschema.validate(request, schema)

    def test_all_abstain_fixture_validates_against_response_schema(self) -> None:
        request = build_llm_adjudication_request(self._packet())
        response = empty_llm_refinement_response(request)
        validate_llm_refinement_decisions(response, request)
        self.assertEqual(response["schema_version"], RESPONSE_SCHEMA_VERSION)
        schema = json.loads(
            Path("src/prideqc/data/llm-refinement-decision-v1.schema.json").read_text(
                encoding="utf-8"
            )
        )
        jsonschema.validate(response, schema)
        self.assertTrue(all(item["decision"] == "abstain" for item in response["decisions"]))

    def test_accept_requires_exact_candidate_value(self) -> None:
        request = build_llm_adjudication_request(self._packet())
        response = empty_llm_refinement_response(request)
        response["decisions"][0] = {
            "decision_id": request["decisions"][0]["decision_id"],
            "decision": "accept",
            "selected_value": "9 ppm",
            "reason": "Invented value",
        }
        with self.assertRaisesRegex(ValueError, "exact supplied candidate"):
            validate_llm_refinement_decisions(response, request)

    def test_accept_supplied_candidate_passes(self) -> None:
        request = build_llm_adjudication_request(self._packet())
        response = empty_llm_refinement_response(request)
        first = request["decisions"][0]
        response["decisions"][0] = {
            "decision_id": first["decision_id"],
            "decision": "accept",
            "selected_value": first["candidate_values"][0],
            "reason": "The supplied evidence supports the supplied candidate.",
        }
        validate_llm_refinement_decisions(response, request)

    def test_reject_and_abstain_cannot_select_a_value(self) -> None:
        request = build_llm_adjudication_request(self._packet())
        response = empty_llm_refinement_response(request)
        response["decisions"][0]["selected_value"] = request["decisions"][0][
            "candidate_values"
        ][0]
        with self.assertRaisesRegex(ValueError, "null selected_value"):
            validate_llm_refinement_decisions(response, request)

    def test_response_requires_exact_decision_coverage_and_request_id(self) -> None:
        request = build_llm_adjudication_request(self._packet())
        response = empty_llm_refinement_response(request)
        missing = copy.deepcopy(response)
        missing["decisions"].pop()
        with self.assertRaisesRegex(ValueError, "missing requested decisions"):
            validate_llm_refinement_decisions(missing, request)

        wrong_request = copy.deepcopy(response)
        wrong_request["request_id"] = "sha256:" + "0" * 64
        with self.assertRaisesRegex(ValueError, "do not match"):
            validate_llm_refinement_decisions(wrong_request, request)

        duplicate = copy.deepcopy(response)
        duplicate["decisions"].append(copy.deepcopy(duplicate["decisions"][0]))
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            validate_llm_refinement_decisions(duplicate, request)


if __name__ == "__main__":
    unittest.main()
