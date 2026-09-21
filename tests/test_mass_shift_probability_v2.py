from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path

SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "benchmarks"
    / "calibrate_mass_shift_probability_v2.py"
)
SPEC = importlib.util.spec_from_file_location("calibrate_mass_shift_probability_v2", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class MassShiftProbabilityV2Tests(unittest.TestCase):
    def cluster(self, *, key: str = "c1", delta: float = 79.9665):
        return MODULE.v1.ClusterEvidence(
            accession="PXDTEST",
            data_file="run.raw",
            cluster_key=key,
            cluster_rank=1,
            delta_mass_da=delta,
            cluster_sigma_da=0.001,
            cluster_min_da=delta - 0.002,
            cluster_max_da=delta + 0.002,
            pair_support=120,
            unique_spectrum_support=180,
            median_spectral_similarity=0.90,
            classification="putative-ptm",
            confidence="high-support",
            match_tolerance_da=0.02,
            support_fraction=0.02,
            diagnostic_candidate_count=2,
            source_tsv="x.tsv",
        )

    def candidate(self, cluster, accession, name, mass, residual):
        return MODULE.v1.CandidateObservation(
            cluster=cluster,
            unimod_accession=accession,
            unimod_name=name,
            theoretical_mass_da=mass,
            residual_da=residual,
            candidate_category="biological-ptm",
            score=0.0,
        )

    def test_one_best_annotation_per_cluster(self) -> None:
        cluster = self.cluster()
        rows = [
            self.candidate(cluster, "UniMod:21", "Phosphorylation", 79.966331, 0.000169),
            self.candidate(cluster, "UniMod:40", "O-Sulfonation", 79.9568, 0.0097),
        ]
        annotations = MODULE.biological_cluster_annotations(rows)
        self.assertEqual(len(annotations), 1)
        self.assertEqual(annotations[0].unimod_accession, "UniMod:21")
        self.assertGreater(annotations[0].candidate_dominance, 0.8)

    def test_exact_isobaric_candidates_are_ambiguous(self) -> None:
        cluster = self.cluster(delta=100.016044)
        annotation = MODULE._best_candidate(
            cluster=cluster,
            candidates=[
                ("UniMod:64", "Succinyl", 100.016044, 0.0),
                ("UniMod:914", "Methylmalonylation", 100.016044, 0.0),
            ],
        )
        self.assertAlmostEqual(annotation.candidate_dominance, 0.5)
        self.assertAlmostEqual(annotation.residual_margin_da or 0.0, 0.0)

    def test_dense_offsets_are_deterministic_and_away_from_zero(self) -> None:
        first = MODULE.deterministic_decoy_offsets(
            128,
            minimum_offset_da=0.05,
            maximum_offset_da=5.0,
        )
        second = MODULE.deterministic_decoy_offsets(
            128,
            minimum_offset_da=0.05,
            maximum_offset_da=5.0,
        )
        self.assertEqual(first, second)
        self.assertEqual(len(first), 128)
        self.assertTrue(all(abs(value) > 0.05 for value in first))
        self.assertTrue(any(value < 0 for value in first))
        self.assertTrue(any(value > 0 for value in first))

    def test_zero_decoy_bin_does_not_reach_v1_point_995_floor(self) -> None:
        model = MODULE.fit_probability_model(
            [10.0] * 25,
            [],
            decoy_replicates=128,
            decoy_min_offset_da=0.05,
            decoy_max_offset_da=5.0,
            catalog_hash="x",
            biological_catalog_records=10,
            training_accessions=("PXD1",),
            minimum_targets_per_bin=25,
        )
        probability = MODULE.probability_for_score(model, 10.0)
        self.assertAlmostEqual(probability, 0.98)
        self.assertLess(probability, 0.995)

    def test_q_values_are_not_worse_for_higher_scores(self) -> None:
        targets = [10.0] * 30 + [5.0] * 70
        decoys = [9.0] * 64 + [4.0] * 512
        model = MODULE.fit_probability_model(
            targets,
            decoys,
            decoy_replicates=128,
            decoy_min_offset_da=0.05,
            decoy_max_offset_da=5.0,
            catalog_hash="x",
            biological_catalog_records=10,
            training_accessions=("PXD1",),
            minimum_targets_per_bin=25,
        )
        self.assertLessEqual(
            MODULE.q_value_for_score(model, 10.0),
            MODULE.q_value_for_score(model, 5.0),
        )

    def test_shortlist_requires_candidate_dominance(self) -> None:
        base = {
            "pxd_accession": "PXDTEST",
            "data_file": "run1.raw",
            "unimod_accession": "UniMod:21",
            "unimod_name": "Phosphorylation",
            "mass_annotation_probability": 0.99,
            "annotation_q_value": 0.01,
            "candidate_mass_dominance": 0.90,
            "unimod_residual_da": 0.001,
            "pair_support": 100,
            "median_spectral_similarity": 0.90,
        }
        ambiguous = {
            **base,
            "data_file": "run3.raw",
            "unimod_accession": "UniMod:64",
            "unimod_name": "Succinyl",
            "candidate_mass_dominance": 0.50,
        }
        rows = [base, {**base, "data_file": "run2.raw"}, ambiguous]
        suggestions = MODULE.build_sdrf_suggestions(
            rows,
            successful_runs_by_accession={"PXDTEST": {"run1.raw", "run2.raw", "run3.raw"}},
            probability_threshold=0.95,
            q_value_threshold=0.05,
            minimum_candidate_dominance=0.80,
            minimum_confident_runs=2,
            minimum_confident_run_fraction=0.5,
            maximum_candidates_per_accession=3,
        )
        self.assertEqual(len(suggestions), 1)
        self.assertEqual(suggestions[0]["unimod_accession"], "UniMod:21")
        self.assertEqual(suggestions[0]["sdrf_status"], "review-required")

    def test_probability_semantics_reject_identity_claim(self) -> None:
        model = MODULE.fit_probability_model(
            [8.0, 7.0],
            [6.0],
            decoy_replicates=16,
            decoy_min_offset_da=0.05,
            decoy_max_offset_da=5.0,
            catalog_hash="x",
            biological_catalog_records=2,
            training_accessions=("PXD1",),
            minimum_targets_per_bin=1,
        )
        self.assertIn("not peptide/site/PTM-identity", model.probability_semantics)


if __name__ == "__main__":
    unittest.main()
