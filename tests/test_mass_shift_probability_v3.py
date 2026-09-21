from __future__ import annotations

import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path

SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "benchmarks"
    / "calibrate_mass_shift_probability_v3.py"
)
SPEC = importlib.util.spec_from_file_location("calibrate_mass_shift_probability_v3", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class MassShiftProbabilityV3Tests(unittest.TestCase):
    def cluster(
        self,
        accession: str,
        run: str,
        rank: int,
        delta: float,
        *,
        pair_support: int = 100,
        similarity: float = 0.9,
        classification: str = "putative-ptm",
    ):
        return MODULE.v1.ClusterEvidence(
            accession=accession,
            data_file=run,
            cluster_key=f"{accession}:{run}:{rank}:{delta}",
            cluster_rank=rank,
            delta_mass_da=delta,
            cluster_sigma_da=0.001,
            cluster_min_da=delta - 0.002,
            cluster_max_da=delta + 0.002,
            pair_support=pair_support,
            unique_spectrum_support=pair_support + 20,
            median_spectral_similarity=similarity,
            classification=classification,
            confidence="high-support",
            match_tolerance_da=0.02,
            support_fraction=0.02,
            diagnostic_candidate_count=2,
            source_tsv=f"{run}.mass-shifts.tsv",
        )

    def catalog(self, accession: str, name: str, mass: float):
        return MODULE.v1.CatalogRecord(accession, name, mass, "biological-ptm")

    def test_family_aggregation_counts_each_run_once(self) -> None:
        clusters = [
            self.cluster("PXD1", "a.raw", 1, 79.9662, pair_support=200),
            self.cluster("PXD1", "a.raw", 2, 79.9670, pair_support=10),
            self.cluster("PXD1", "b.raw", 1, 79.9664),
            self.cluster("PXD1", "c.raw", 1, 42.0106),
        ]
        families = MODULE.aggregate_mass_families(
            clusters,
            successful_runs_by_accession={"PXD1": {"a.raw", "b.raw", "c.raw"}},
            family_tolerance_da=0.02,
        )
        phospho = min(families, key=lambda item: abs(item.median_mass_da - 79.9663))
        self.assertEqual(phospho.run_count, 2)
        self.assertAlmostEqual(phospho.run_prevalence, 2 / 3)
        self.assertEqual(len(phospho.run_names), 2)

    def test_beta_superiority_is_half_for_equal_uniform_posteriors(self) -> None:
        probability = MODULE.beta_superiority_probability(0, 0, 0, 0)
        self.assertAlmostEqual(probability, 0.5, places=12)

    def test_beta_superiority_rewards_study_specific_prevalence(self) -> None:
        probability = MODULE.beta_superiority_probability(9, 10, 1, 30)
        self.assertGreater(probability, 0.999)
        reverse = MODULE.beta_superiority_probability(1, 30, 9, 10)
        self.assertLess(reverse, 0.001)

    def test_background_prevalence_excludes_target_accession(self) -> None:
        family = MODULE.MassFamily(
            accession="PXD1",
            family_id="PXD1:MF0001",
            median_mass_da=79.9663,
            mass_mad_da=0.0001,
            run_count=2,
            successful_runs=2,
            run_prevalence=1.0,
            median_cluster_rank=2.0,
            median_pair_support=100.0,
            median_unique_spectrum_support=120.0,
            median_spectral_similarity=0.9,
            median_support_fraction=0.02,
            median_cluster_sigma_da=0.001,
            run_names=("a.raw", "b.raw"),
            source_cluster_keys=("a", "b"),
        )
        counts = MODULE.prevalence_counts_for_family(
            family,
            run_masses_by_accession={
                "PXD1": {"a.raw": [79.9663], "b.raw": [79.9662]},
                "PXD2": {"c.raw": [14.0157], "d.raw": [79.9664]},
            },
            successful_runs_by_accession={
                "PXD1": {"a.raw", "b.raw"},
                "PXD2": {"c.raw", "d.raw"},
            },
            family_tolerance_da=0.02,
        )
        self.assertEqual(counts, (2, 2, 1, 2))

    def test_cross_run_centroid_separates_phospho_from_osulfonation(self) -> None:
        family = MODULE.MassFamily(
            accession="PXD1",
            family_id="PXD1:MF0001",
            median_mass_da=79.96635,
            mass_mad_da=0.00015,
            run_count=20,
            successful_runs=20,
            run_prevalence=1.0,
            median_cluster_rank=5.0,
            median_pair_support=500.0,
            median_unique_spectrum_support=700.0,
            median_spectral_similarity=0.9,
            median_support_fraction=0.02,
            median_cluster_sigma_da=0.001,
            run_names=tuple(f"r{i}.raw" for i in range(20)),
            source_cluster_keys=tuple(f"c{i}" for i in range(20)),
        )
        annotation = MODULE.annotate_family(
            family,
            [
                self.catalog("UniMod:21", "Phosphorylation", 79.966331),
                self.catalog("UniMod:40", "O-Sulfonation", 79.9568),
            ],
            annotation_window_da=0.02,
            systematic_mass_floor_da=0.0015,
        )
        assert annotation is not None
        self.assertEqual(annotation.unimod_accession, "UniMod:21")
        self.assertGreater(annotation.candidate_mass_dominance, 0.99)

    def test_exact_isobaric_family_remains_ambiguous(self) -> None:
        family = MODULE.MassFamily(
            accession="PXD1",
            family_id="PXD1:MF0001",
            median_mass_da=100.016044,
            mass_mad_da=0.0001,
            run_count=10,
            successful_runs=10,
            run_prevalence=1.0,
            median_cluster_rank=10.0,
            median_pair_support=50.0,
            median_unique_spectrum_support=80.0,
            median_spectral_similarity=0.8,
            median_support_fraction=0.01,
            median_cluster_sigma_da=0.001,
            run_names=tuple(f"r{i}.raw" for i in range(10)),
            source_cluster_keys=tuple(f"c{i}" for i in range(10)),
        )
        annotation = MODULE.annotate_family(
            family,
            [
                self.catalog("UniMod:64", "Succinyl", 100.016044),
                self.catalog("UniMod:914", "Methylmalonylation", 100.016044),
            ],
            annotation_window_da=0.02,
            systematic_mass_floor_da=0.0015,
        )
        assert annotation is not None
        self.assertAlmostEqual(annotation.candidate_mass_dominance, 0.5)

    def test_review_requires_effect_size_and_mass_dominance(self) -> None:
        base = {
            "pxd_accession": "PXD1",
            "family_id": "PXD1:MF0001",
            "median_mass_da": 79.9663,
            "family_runs": 9,
            "successful_runs": 10,
            "family_run_prevalence": 0.9,
            "background_prevalence": 0.05,
            "prevalence_difference": 0.85,
            "prevalence_ratio": 18.0,
            "study_enrichment_probability": 0.999,
            "study_enrichment_q_value": 0.001,
            "unimod_accession": "UniMod:21",
            "unimod_name": "Phosphorylation",
            "unimod_residual_da": 0.0001,
            "candidate_mass_dominance": 0.99,
            "biological_candidate_count": 2,
            "median_pair_support": 100.0,
            "median_unique_spectrum_support": 120.0,
            "median_spectral_similarity": 0.9,
        }
        ubiquitous = {
            **base,
            "family_id": "PXD1:MF0002",
            "unimod_accession": "UniMod:34",
            "unimod_name": "Methylation",
            "background_prevalence": 0.8,
            "prevalence_difference": 0.1,
            "prevalence_ratio": 1.125,
        }
        ambiguous = {
            **base,
            "family_id": "PXD1:MF0003",
            "unimod_accession": "UniMod:64",
            "unimod_name": "Succinyl",
            "candidate_mass_dominance": 0.5,
        }
        review = MODULE.build_sdrf_review(
            [base, ubiquitous, ambiguous],
            semantic_evidence={},
            enrichment_probability_threshold=0.95,
            enrichment_q_value_threshold=0.05,
            minimum_prevalence_difference=0.20,
            minimum_prevalence_ratio=2.0,
            minimum_accession_prevalence=0.10,
            minimum_family_runs=2,
            minimum_candidate_dominance=0.80,
            maximum_candidates_per_accession=3,
        )
        self.assertEqual(len(review), 1)
        self.assertEqual(review[0]["unimod_accession"], "UniMod:21")
        self.assertEqual(review[0]["sdrf_status"], "semantic-evidence-required")

    def test_semantic_evidence_is_independent_gate(self) -> None:
        row = {
            "pxd_accession": "PXD1",
            "family_id": "PXD1:MF0001",
            "median_mass_da": 79.9663,
            "family_runs": 9,
            "successful_runs": 10,
            "family_run_prevalence": 0.9,
            "background_prevalence": 0.05,
            "prevalence_difference": 0.85,
            "prevalence_ratio": 18.0,
            "study_enrichment_probability": 0.999,
            "study_enrichment_q_value": 0.001,
            "unimod_accession": "UniMod:21",
            "unimod_name": "Phosphorylation",
            "unimod_residual_da": 0.0001,
            "candidate_mass_dominance": 0.99,
            "biological_candidate_count": 2,
            "median_pair_support": 100.0,
            "median_unique_spectrum_support": 120.0,
            "median_spectral_similarity": 0.9,
        }
        evidence = {
            ("PXD1", "UniMod:21"): {
                "semantic_evidence_status": "supported",
                "semantic_evidence_sources": "publication",
                "semantic_evidence_notes": "phosphopeptide enrichment",
            }
        }
        review = MODULE.build_sdrf_review(
            [row],
            semantic_evidence=evidence,
            enrichment_probability_threshold=0.95,
            enrichment_q_value_threshold=0.05,
            minimum_prevalence_difference=0.20,
            minimum_prevalence_ratio=2.0,
            minimum_accession_prevalence=0.10,
            minimum_family_runs=2,
            minimum_candidate_dominance=0.80,
            maximum_candidates_per_accession=3,
        )
        self.assertEqual(review[0]["sdrf_status"], "review-required")
        self.assertEqual(review[0]["semantic_evidence_status"], "supported")

    def test_study_evidence_parser_detects_mixed_support(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "evidence.tsv"
            path.write_text(
                "pxd_accession\tunimod_accession\tevidence_status\tevidence_source\tevidence_note\n"
                "PXD1\tUniMod:21\tsupported\tpublication\tphospho enrichment\n"
                "PXD1\tUniMod:21\tconflicting\tsearch\tsearch file unclear\n"
            )
            evidence = MODULE.read_study_evidence(path)
        item = evidence[("PXD1", "UniMod:21")]
        self.assertEqual(item["semantic_evidence_status"], "mixed")
        self.assertIn("publication", item["semantic_evidence_sources"])
        self.assertIn("search", item["semantic_evidence_sources"])


if __name__ == "__main__":
    unittest.main()
