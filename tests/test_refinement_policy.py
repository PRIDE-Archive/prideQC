from __future__ import annotations

import copy
import unittest
from typing import Any, cast

from prideqc.refinement_policy import (
    pre_adjudication_policy_metadata,
    resolve_pre_adjudication_policy,
)


class RefinementPolicyTests(unittest.TestCase):
    def _tolerance(
        self,
        original: str,
        *,
        scope_status: str = "accession-complete",
    ) -> dict[str, object]:
        estimates = [7.8, 8.0, 7.9, 8.0]
        return {
            "decision_id": "Experiment group 1:precursor-mass-tolerance",
            "decision_type": "mass_tolerance",
            "candidate_values": ["8 ppm"],
            "parameter_scope_status": scope_status,
            "original": {"rows": [{"row": 2, "values": [original]}]},
            "evidence": {
                "group_median": 7.95,
                "group_common_max": 8.0,
                "supporting_runs": 4,
                "group_runs": 4,
                "selection_policy": "maximum-supported-per-run-estimate-rounded-up",
                "per_run_estimates": [
                    {
                        "run_id": f"run-{index}.raw",
                        "status": "available",
                        "value": value,
                        "unit": "ppm",
                        "confidence": "high" if index < 3 else "moderate",
                        "resolution_regime": "high-resolution",
                    }
                    for index, value in enumerate(estimates)
                ],
            },
        }

    def _ptm(
        self,
        *,
        scope_status: str = "accession-complete",
    ) -> dict[str, object]:
        observations = []
        for index in range(10):
            observations.append(
                {
                    "run_id": f"run-{index}.raw",
                    "observations": [
                        {
                            "delta_mass_da": 14.0155,
                            "pair_support": 80,
                            "unique_spectrum_support": 50,
                            "median_spectral_similarity": 0.8,
                            "classification": "putative-ptm",
                            "confidence": "high-support",
                            "match_tolerance_da": 0.02,
                            "candidate_options": [
                                {
                                    "accession": "UniMod:34",
                                    "name": "Methylation",
                                    "category": "biological-ptm",
                                    "residual_da": 0.001,
                                }
                            ],
                        }
                    ],
                }
            )
        return {
            "decision_id": "Experiment group 1:modification-family:14.015500",
            "decision_type": "modification",
            "candidate_values": ["NT=Methylation;AC=UniMod:34"],
            "parameter_scope_status": scope_status,
            "original": {"rows": [{"row": 2, "values": []}]},
            "evidence": {
                "supporting_runs": 10,
                "group_runs": 10,
                "run_prevalence": 1.0,
                "high_support_run_fraction": 1.0,
                "candidate_run_fraction": 1.0,
                "recurrent_family_probability": 0.999999,
                "mass_identity_ambiguous": False,
                "candidate_options": [
                    {
                        "accession": "UniMod:34",
                        "name": "Methylation",
                        "semantic_evidence": {
                            "status": "not-evaluated",
                            "sources": [],
                            "notes": [],
                        },
                    }
                ],
                "semantic_evidence_status": "not-evaluated",
                "per_run_support": observations,
            },
        }

    def test_policy_v3_records_reconstruction_thresholds(self) -> None:
        metadata = pre_adjudication_policy_metadata()
        self.assertEqual(metadata["version"], "prideqc-pre-adjudication-policy-v3")
        self.assertTrue(
            metadata["raw_precision_can_fill_missing_reported_tolerance_when_confident"]
        )
        self.assertTrue(
            metadata["raw_ptm_identity_is_model_eligibility_not_automatic_acceptance"]
        )
        self.assertFalse(metadata["repository_summary_metadata_is_hard_negative_evidence"])

    def test_existing_reported_tolerance_matching_candidate_is_noop_accept(self) -> None:
        resolution = resolve_pre_adjudication_policy(self._tolerance("8 ppm"))
        self.assertIsNotNone(resolution)
        assert resolution is not None
        self.assertEqual(resolution.decision, "accept")
        self.assertEqual(resolution.selected_value, "8 ppm")
        self.assertEqual(resolution.rule, "tolerance-original-already-matches")

    def test_existing_reported_tolerance_is_preserved_when_candidate_differs(self) -> None:
        resolution = resolve_pre_adjudication_policy(self._tolerance("20 ppm"))
        self.assertIsNotNone(resolution)
        assert resolution is not None
        self.assertEqual(resolution.decision, "reject")
        self.assertIsNone(resolution.selected_value)
        self.assertEqual(resolution.rule, "tolerance-preserve-reported-value")

    def test_missing_stable_tolerance_is_reconstructed_for_verified_scope(self) -> None:
        resolution = resolve_pre_adjudication_policy(self._tolerance("not available"))
        self.assertIsNotNone(resolution)
        assert resolution is not None
        self.assertEqual(resolution.decision, "accept")
        self.assertEqual(resolution.selected_value, "8 ppm")
        self.assertEqual(
            resolution.rule,
            "tolerance-missing-high-confidence-reconstruction",
        )
        self.assertGreaterEqual(resolution.metrics["coverage"], 0.8)
        self.assertLessEqual(resolution.metrics["relative_mad"], 0.10)

    def test_missing_tolerance_abstains_for_unverified_partial_scope(self) -> None:
        decision = self._tolerance(
            "not available",
            scope_status="cohort-group-unverified",
        )
        resolution = resolve_pre_adjudication_policy(decision)
        self.assertIsNotNone(resolution)
        assert resolution is not None
        self.assertEqual(resolution.decision, "abstain")
        self.assertEqual(resolution.rule, "tolerance-missing-insufficient-confidence")
        self.assertEqual(
            resolution.metrics["parameter_scope_status"],
            "cohort-group-unverified",
        )

    def test_high_confidence_raw_ptm_is_model_eligible_for_verified_scope(self) -> None:
        self.assertIsNone(resolve_pre_adjudication_policy(self._ptm()))

    def test_high_confidence_raw_ptm_abstains_for_unverified_partial_scope(self) -> None:
        resolution = resolve_pre_adjudication_policy(
            self._ptm(scope_status="cohort-group-unverified")
        )
        self.assertIsNotNone(resolution)
        assert resolution is not None
        self.assertEqual(resolution.decision, "abstain")
        self.assertEqual(resolution.rule, "ptm-parameter-scope-unverified")

    def test_project_metadata_conflict_is_not_a_deterministic_veto(self) -> None:
        decision = self._ptm()
        decision["evidence"]["semantic_evidence_status"] = "conflicting"  # type: ignore[index]
        self.assertIsNone(resolve_pre_adjudication_policy(decision))

    def test_mass_only_ptm_below_identity_gate_abstains_without_model(self) -> None:
        decision = self._ptm()
        decision["evidence"]["high_support_run_fraction"] = 0.8  # type: ignore[index]
        resolution = resolve_pre_adjudication_policy(decision)
        self.assertIsNotNone(resolution)
        assert resolution is not None
        self.assertEqual(resolution.decision, "abstain")
        self.assertEqual(resolution.rule, "ptm-insufficient-identity-confidence")
        self.assertFalse(resolution.metrics["raw_identity_gate_met"])

    def test_semantically_supported_borderline_ptm_is_left_for_model(self) -> None:
        decision = self._ptm()
        decision["evidence"]["high_support_run_fraction"] = 0.8  # type: ignore[index]
        decision["evidence"]["semantic_evidence_status"] = "supported"  # type: ignore[index]
        self.assertIsNone(resolve_pre_adjudication_policy(decision))

    def test_ptm_already_reported_is_noop_accept(self) -> None:
        decision = self._ptm()
        decision["original"] = {
            "rows": [{"row": 2, "values": ["NT=Methylation;AC=UniMod:34"]}]
        }
        resolution = resolve_pre_adjudication_policy(decision)
        self.assertIsNotNone(resolution)
        assert resolution is not None
        self.assertEqual(resolution.decision, "accept")
        self.assertEqual(resolution.rule, "ptm-original-already-reported")

    def test_ptm_mass_ambiguity_always_abstains(self) -> None:
        decision = self._ptm()
        decision["evidence"]["mass_identity_ambiguous"] = True  # type: ignore[index]
        resolution = resolve_pre_adjudication_policy(decision)
        self.assertIsNotNone(resolution)
        assert resolution is not None
        self.assertEqual(resolution.decision, "abstain")
        self.assertEqual(resolution.rule, "ptm-mass-identity-ambiguous")

    def test_high_confidence_requires_biological_ptm_classification(self) -> None:
        decision = self._ptm()
        weakened = copy.deepcopy(decision)
        evidence = cast(dict[str, Any], weakened["evidence"])
        per_run = cast(list[dict[str, Any]], evidence["per_run_support"])
        observation = cast(dict[str, Any], per_run[0]["observations"][0])
        option = cast(dict[str, Any], observation["candidate_options"][0])
        option["category"] = "sample-prep-or-artifact"
        resolution = resolve_pre_adjudication_policy(weakened)
        self.assertIsNotNone(resolution)
        assert resolution is not None
        self.assertEqual(resolution.decision, "abstain")


if __name__ == "__main__":
    unittest.main()
