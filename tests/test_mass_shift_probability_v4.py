from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "benchmarks" / "calibrate_mass_shift_probability_v4.py"
SPEC = importlib.util.spec_from_file_location("mass_shift_probability_v4", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
v4 = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = v4
SPEC.loader.exec_module(v4)


class _FakeModification:
    def __init__(self, classification: str, origin: str = "K") -> None:
        self.classification = classification
        self.origin = origin

    def getUniModAccession(self):
        return "UniMod:1"

    def getDiffMonoMass(self):
        return 42.010565

    def getFullName(self):
        return "Acetyl"

    def getId(self):
        return "Acetyl"

    def getSourceClassification(self):
        return self.classification

    def getSourceClassificationName(self, value):
        return value

    def getOrigin(self):
        return self.origin

    def getTermSpecificity(self):
        return "Anywhere"

    def getTermSpecificityName(self, value):
        return value


class _FakeDB:
    def __init__(self) -> None:
        self.mods = [
            _FakeModification("Post-translational", "K"),
            _FakeModification("Chemical derivative", "N-term"),
        ]

    def getNumberOfModifications(self):
        return len(self.mods)

    def getModification(self, index):
        return self.mods[index]


class _FakeOMS:
    @staticmethod
    def ModificationsDB():
        return _FakeDB()


def _cluster(
    accession: str,
    data_file: str,
    mass: float,
    classification: str = "putative-ptm",
    rank: int = 1,
):
    return SimpleNamespace(
        accession=accession,
        data_file=data_file,
        delta_mass_da=mass,
        classification=classification,
        pair_support=50,
        unique_spectrum_support=80,
        cluster_rank=rank,
        cluster_key=f"{accession}:{data_file}:{rank}",
        median_spectral_similarity=0.8,
        support_fraction=0.01,
        cluster_sigma_da=0.001,
    )


class ProbabilityV4Tests(unittest.TestCase):
    def test_catalog_preserves_mixed_specificity_classifications(self):
        records = v4.catalog_from_openms_specificity_aware(_FakeOMS)
        self.assertEqual(len(records), 1)
        record = records[0]
        self.assertEqual(record.accession, "UniMod:1")
        self.assertEqual(
            set(record.source_classifications),
            {"Post-translational", "Chemical derivative"},
        )
        self.assertTrue(record.has_biological_specificity)
        self.assertEqual(set(record.origins), {"K", "N-term"})

    def test_beta_raw_prevalence_probability_is_monotone(self):
        low = v4.beta_tail_probability_above(2, 10, 0.10)
        high = v4.beta_tail_probability_above(8, 10, 0.10)
        self.assertGreater(high, low)
        self.assertGreater(high, 0.95)

    def test_family_aggregation_includes_sample_prep_classification(self):
        clusters = [
            _cluster("PXD1", "a.raw", 100.0160, "sample-prep-modification"),
            _cluster("PXD1", "b.raw", 100.0162, "sample-prep-modification"),
        ]
        families = v4.aggregate_recurrent_mass_families(
            clusters,
            successful_runs_by_accession={"PXD1": {"a.raw", "b.raw"}},
            family_tolerance_da=0.02,
        )
        self.assertEqual(len(families), 1)
        self.assertEqual(families[0].run_count, 2)
        self.assertAlmostEqual(families[0].median_mass_da, 100.0161, places=4)

    def test_duplicate_clusters_from_same_run_do_not_inflate_prevalence(self):
        clusters = [
            _cluster("PXD1", "a.raw", 79.9662, rank=1),
            _cluster("PXD1", "a.raw", 79.9664, rank=2),
            _cluster("PXD1", "b.raw", 79.9663, rank=1),
        ]
        families = v4.aggregate_recurrent_mass_families(
            clusters,
            successful_runs_by_accession={"PXD1": {"a.raw", "b.raw", "c.raw"}},
            family_tolerance_da=0.02,
        )
        self.assertEqual(families[0].run_count, 2)
        self.assertAlmostEqual(families[0].run_prevalence, 2 / 3)

    def test_semantic_supported_plus_raw_family_is_review_required(self):
        family = v4.v3.MassFamily(
            accession="PXD1",
            family_id="PXD1:RF0001",
            median_mass_da=42.0104,
            mass_mad_da=0.0002,
            run_count=9,
            successful_runs=10,
            run_prevalence=0.9,
            median_cluster_rank=3,
            median_pair_support=40,
            median_unique_spectrum_support=70,
            median_spectral_similarity=0.8,
            median_support_fraction=0.01,
            median_cluster_sigma_da=0.001,
            run_names=tuple(f"r{i}.raw" for i in range(9)),
            source_cluster_keys=(),
        )
        evidence = {
            ("PXD1", "UniMod:1"): v4.EvidenceRecord(
                "PXD1", "UniMod:1", "supported", ("publication",), ("Kac study",)
            )
        }
        catalog = [
            v4.CatalogRecord(
                "UniMod:1",
                "Acetyl",
                42.010565,
                ("Post-translational", "Chemical derivative"),
                ("K",),
                ("Anywhere",),
            )
        ]
        clusters = [_cluster("PXD1", f"r{i}.raw", 42.0104) for i in range(9)]
        rows = v4.build_semantic_review(
            evidence=evidence,
            families=[family],
            clusters=clusters,
            successful_runs_by_accession={"PXD1": {f"r{i}.raw" for i in range(10)}},
            catalog=catalog,
            raw_match_window_da=0.02,
            minimum_raw_prevalence=0.10,
            raw_prevalence_probability_threshold=0.95,
            minimum_family_runs=2,
            family_tolerance_da=0.02,
        )
        self.assertEqual(rows[0]["raw_confirmation_status"], "confirmed")
        self.assertEqual(rows[0]["sdrf_status"], "review-required")
        self.assertEqual(rows[0]["unimod_name"], "Acetyl")

    def test_semantic_supported_exact_isobar_keeps_ambiguity_visible(self):
        family = v4.v3.MassFamily(
            accession="PXD1",
            family_id="PXD1:RF0001",
            median_mass_da=100.016044,
            mass_mad_da=0.0001,
            run_count=5,
            successful_runs=10,
            run_prevalence=0.5,
            median_cluster_rank=10,
            median_pair_support=20,
            median_unique_spectrum_support=30,
            median_spectral_similarity=0.7,
            median_support_fraction=0.01,
            median_cluster_sigma_da=0.002,
            run_names=tuple(f"r{i}.raw" for i in range(5)),
            source_cluster_keys=(),
        )
        evidence = {
            ("PXD1", "UniMod:64"): v4.EvidenceRecord(
                "PXD1", "UniMod:64", "supported", ("publication",), ()
            )
        }
        catalog = [
            v4.CatalogRecord("UniMod:64", "Succinyl", 100.016044, ("Chemical derivative",), ("K",), ("Anywhere",)),
            v4.CatalogRecord("UniMod:914", "Methylmalonylation", 100.016044, ("Post-translational",), ("S",), ("Anywhere",)),
        ]
        clusters = [_cluster("PXD1", f"r{i}.raw", 100.016044, "sample-prep-modification") for i in range(5)]
        rows = v4.build_semantic_review(
            evidence=evidence,
            families=[family],
            clusters=clusters,
            successful_runs_by_accession={"PXD1": {f"r{i}.raw" for i in range(10)}},
            catalog=catalog,
            raw_match_window_da=0.02,
            minimum_raw_prevalence=0.10,
            raw_prevalence_probability_threshold=0.95,
            minimum_family_runs=2,
            family_tolerance_da=0.02,
        )
        self.assertEqual(rows[0]["sdrf_status"], "review-required")
        self.assertEqual(rows[0]["mass_identity_ambiguous"], 1)
        self.assertIn("UniMod:914", rows[0]["mass_compatible_alternative_accessions"])

    def test_supported_but_missing_raw_family_is_not_review_required(self):
        evidence = {
            ("PXD1", "UniMod:64"): v4.EvidenceRecord(
                "PXD1", "UniMod:64", "supported", ("publication",), ()
            )
        }
        catalog = [
            v4.CatalogRecord("UniMod:64", "Succinyl", 100.016044, (), (), ())
        ]
        rows = v4.build_semantic_review(
            evidence=evidence,
            families=[],
            clusters=[],
            successful_runs_by_accession={"PXD1": {"a.raw", "b.raw"}},
            catalog=catalog,
            raw_match_window_da=0.02,
            minimum_raw_prevalence=0.10,
            raw_prevalence_probability_threshold=0.95,
            minimum_family_runs=2,
            family_tolerance_da=0.02,
        )
        self.assertEqual(rows[0]["raw_confirmation_status"], "not-detected")
        self.assertEqual(rows[0]["sdrf_status"], "semantic-supported-raw-unconfirmed")

    def test_conflicting_semantics_holds_even_with_raw_support(self):
        family = v4.v3.MassFamily(
            accession="PXD1", family_id="PXD1:RF0001", median_mass_da=79.966331,
            mass_mad_da=0.0, run_count=10, successful_runs=10, run_prevalence=1.0,
            median_cluster_rank=1, median_pair_support=100, median_unique_spectrum_support=150,
            median_spectral_similarity=0.9, median_support_fraction=0.02,
            median_cluster_sigma_da=0.001, run_names=(), source_cluster_keys=(),
        )
        evidence = {("PXD1", "UniMod:21"): v4.EvidenceRecord("PXD1", "UniMod:21", "conflicting", (), ())}
        catalog = [v4.CatalogRecord("UniMod:21", "Phosphorylation", 79.966331, ("Post-translational",), ("S",), ("Anywhere",))]
        rows = v4.build_semantic_review(
            evidence=evidence, families=[family], clusters=[],
            successful_runs_by_accession={"PXD1": {f"r{i}" for i in range(10)}},
            catalog=catalog, raw_match_window_da=0.02, minimum_raw_prevalence=0.10,
            raw_prevalence_probability_threshold=0.95, minimum_family_runs=2,
            family_tolerance_da=0.02,
        )
        self.assertEqual(rows[0]["sdrf_status"], "hold-conflicting-evidence")

    def test_no_semantic_evidence_means_no_sdrf_candidates(self):
        rows = v4.build_semantic_review(
            evidence={}, families=[], clusters=[], successful_runs_by_accession={},
            catalog=[], raw_match_window_da=0.02, minimum_raw_prevalence=0.10,
            raw_prevalence_probability_threshold=0.95, minimum_family_runs=2,
            family_tolerance_da=0.02,
        )
        self.assertEqual(rows, [])


if __name__ == "__main__":
    unittest.main()
