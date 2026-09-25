from __future__ import annotations

import json
import unittest
from pathlib import Path
from unittest.mock import patch

from prideqc.refinement_model import (
    MODEL_INPUT_PROJECTION_VERSION,
    SYSTEM_PROMPT_VERSION,
    LocalLlamaCppAdapter,
    build_llama_chat_payload,
)


class RefinementModelScalingTests(unittest.TestCase):
    @staticmethod
    def _request(decision: dict[str, object]) -> dict[str, object]:
        source_hash = "a" * 64
        return {
            "schema_version": "prideqc-llm-adjudication-request-v1",
            "request_scope": "accession-sdrf",
            "request_id": f"sha256:{source_hash}",
            "project_accession": "PXDTEST",
            "input_mode": "no-original-sdrf",
            "source_packet": {
                "schema_version": "prideqc-llm-refinement-packet-v1",
                "sha256": source_hash,
                "prideqc_version": "0.2.0",
            },
            "sdrf": {"source_name": None, "sha256": None},
            "contract": {},
            "decisions": [decision],
        }

    @staticmethod
    def _tolerance_decision(run_count: int) -> dict[str, object]:
        per_run = []
        for index in range(run_count):
            per_run.append(
                {
                    "run_id": f"run-{index:03d}.raw",
                    "status": "available",
                    "source_field": "suggested_precursor_search_tolerance_ppm",
                    "value": 8.0 + (index % 9) * 0.25,
                    "unit": "ppm",
                    "confidence": "high",
                    "support": 1000 + index,
                    "total": 5000 + index,
                    "method": "repeated long method text " * 12,
                    "detail": "repeated long scientific detail " * 20,
                }
            )
        return {
            "decision_id": "Experiment group 1:precursor-mass-tolerance",
            "decision_type": "mass_tolerance",
            "experiment_group": "Experiment group 1",
            "target_field": "comment[precursor mass tolerance]",
            "write_semantics": "evidence_only_no_original_sdrf",
            "target_rows": [],
            "target_runs": [f"run-{index:03d}.raw" for index in range(run_count)],
            "original": {
                "status": "unavailable",
                "reason": "no-original-sdrf",
                "column_present": False,
                "column_count": 0,
                "rows": [],
            },
            "candidate_values": ["12 ppm"],
            "evidence": {
                "unit": "ppm",
                "group_median": 10.5,
                "group_common_max": 11.5,
                "supporting_runs": run_count,
                "group_runs": run_count,
                "source_measurement_support": 111700,
                "selection_policy": "maximum-supported-per-run-estimate-rounded-up",
                "minimum_compatible_run_fraction": 0.8,
                "per_run_estimates": per_run,
                "method": "cohort tolerance synthesis",
                "detail": "The aggregate candidate is conservative and evidence derived.",
            },
            "evidence_semantics": {},
            "local_context": {"ptm_families": []},
        }

    def test_large_tolerance_request_has_bounded_model_projection(self) -> None:
        request = self._request(self._tolerance_decision(169))
        payload = build_llama_chat_payload(request)
        content = payload["messages"][1]["content"]
        model_input = json.loads(content)
        projected = model_input["decision"]

        self.assertEqual(
            set(model_input),
            {"request_id", "project_accession", "decision"},
        )
        self.assertEqual(projected["target_run_count"], 169)
        self.assertNotIn("target_runs", projected)
        self.assertNotIn("per_run_estimates", projected["evidence"])
        summary = projected["evidence"]["per_run_estimate_summary"]
        self.assertEqual(summary["run_count"], 169)
        self.assertEqual(summary["estimate"]["count"], 169)
        self.assertLess(len(content.encode("utf-8")), 12_000)

    def test_large_ptm_request_has_bounded_model_projection(self) -> None:
        per_run_support = []
        for index in range(169):
            per_run_support.append(
                {
                    "run_id": f"run-{index:03d}.raw",
                    "observations": [
                        {
                            "delta_mass_da": 14.0156,
                            "pair_support": 200 + index,
                            "median_spectral_similarity": 0.8,
                            "classification": "biological-ptm",
                            "confidence": "moderate",
                            "residual_da": 0.001,
                            "candidate_options": [
                                {
                                    "accession": "UniMod:34",
                                    "name": "Methylation",
                                    "category": "biological-ptm",
                                }
                            ],
                        }
                    ],
                }
            )

        decision = {
            "decision_id": "Experiment group 1:modification-family:14.015600",
            "decision_type": "modification",
            "experiment_group": "Experiment group 1",
            "target_field": "comment[modification parameters]",
            "write_semantics": "evidence_only_no_original_sdrf",
            "target_rows": [],
            "target_runs": [f"run-{index:03d}.raw" for index in range(169)],
            "original": {
                "status": "unavailable",
                "reason": "no-original-sdrf",
                "column_present": False,
                "column_count": 0,
                "rows": [],
            },
            "candidate_values": ["NT=Methylation;AC=UniMod:34"],
            "evidence": {
                "observed_delta_mass_da": 14.0156,
                "supporting_runs": 169,
                "group_runs": 169,
                "run_prevalence": 1.0,
                "high_support_run_fraction": 0.8,
                "candidate_run_fraction": 1.0,
                "recurrent_family_probability": 0.999,
                "mass_identity_ambiguous": False,
                "semantic_evidence_status": "supported",
                "candidate_options": [
                    {
                        "accession": "UniMod:34",
                        "name": "Methylation",
                        "semantic_evidence": {
                            "status": "supported",
                            "sources": ["study-metadata"],
                            "notes": ["methylation enrichment"],
                        },
                    }
                ],
                "per_run_support": per_run_support,
            },
            "evidence_semantics": {
                "recurrent_family_probability_is_identity_probability": False,
            },
            "local_context": {"ptm_families": []},
        }

        payload = build_llama_chat_payload(self._request(decision))
        content = payload["messages"][1]["content"]
        projected = json.loads(content)["decision"]

        self.assertNotIn("per_run_support", projected["evidence"])
        summary = projected["evidence"]["per_run_support_summary"]
        self.assertEqual(summary["run_count"], 169)
        self.assertEqual(summary["observation_count"], 169)
        self.assertLess(len(content.encode("utf-8")), 12_000)

    def test_empty_request_never_starts_llama(self) -> None:
        request = self._request(self._tolerance_decision(1))
        request["decisions"] = []
        adapter = LocalLlamaCppAdapter(
            server_path=Path("/definitely/missing/llama-server"),
            model_path=Path("/definitely/missing/model.gguf"),
        )

        with patch("prideqc.refinement_model._json_request") as request_mock:
            result = adapter.adjudicate(request)

        request_mock.assert_not_called()
        self.assertEqual(result.decisions["decisions"], [])
        self.assertEqual(result.decisions["input_mode"], "no-original-sdrf")
        self.assertEqual(result.audit["decision_counts"]["total"], 0)
        self.assertEqual(result.audit["decision_counts"]["policy_resolved"], 0)
        self.assertEqual(result.audit["decision_counts"]["model_called"], 0)
        self.assertFalse(result.audit["runtime"]["invoked"])
        self.assertEqual(result.audit["prompt_version"], SYSTEM_PROMPT_VERSION)
        self.assertEqual(
            result.audit["model_input_projection_version"],
            MODEL_INPUT_PROJECTION_VERSION,
        )


if __name__ == "__main__":
    unittest.main()
