from __future__ import annotations

import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path

SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "benchmarks"
    / "calibrate_mass_shift_probability_v1.py"
)
SPEC = importlib.util.spec_from_file_location("calibrate_mass_shift_probability_v1", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class MassShiftProbabilityTests(unittest.TestCase):
    def test_evidence_score_rewards_support_similarity_and_mass_fit(self) -> None:
        weak = MODULE.evidence_score(
            pair_support=8,
            unique_spectrum_support=7,
            median_spectral_similarity=0.35,
            support_fraction=0.001,
            cluster_sigma_da=0.01,
            match_tolerance_da=0.02,
            residual_da=0.015,
        )
        strong = MODULE.evidence_score(
            pair_support=120,
            unique_spectrum_support=180,
            median_spectral_similarity=0.90,
            support_fraction=0.02,
            cluster_sigma_da=0.002,
            match_tolerance_da=0.02,
            residual_da=0.001,
        )
        self.assertGreater(strong, weak)

    def test_target_decoy_probability_is_monotone_with_score(self) -> None:
        targets = [10.0, 9.8, 9.6, 9.4, 8.0, 7.8, 7.6, 6.0, 5.8, 5.6]
        decoys = [7.9, 7.7, 6.1, 5.9, 5.7, 5.5, 5.3, 5.1]
        model = MODULE.fit_probability_model(
            targets,
            decoys,
            decoy_offsets=(0.3, 0.7),
            catalog_hash="abc",
            biological_catalog_records=10,
            training_accessions=("PXD1", "PXD2"),
            minimum_targets_per_bin=2,
        )
        probabilities = [item.probability for item in model.calibration_bins]
        self.assertEqual(probabilities, sorted(probabilities, reverse=True))
        self.assertGreaterEqual(
            MODULE.probability_for_score(model, 10.0),
            MODULE.probability_for_score(model, 5.5),
        )
        self.assertLessEqual(
            MODULE.q_value_for_score(model, 10.0),
            MODULE.q_value_for_score(model, 5.5),
        )

    def test_model_round_trip_preserves_probabilities(self) -> None:
        model = MODULE.fit_probability_model(
            [10.0, 9.0, 8.0, 7.0],
            [7.5, 6.5, 5.5],
            decoy_offsets=(0.3, 0.7),
            catalog_hash="catalog",
            biological_catalog_records=20,
            training_accessions=("PXD1",),
            minimum_targets_per_bin=2,
        )
        restored = MODULE.model_from_json(MODULE.model_to_json(model))
        self.assertEqual(restored.model_version, MODULE.MODEL_VERSION)
        self.assertAlmostEqual(
            MODULE.probability_for_score(model, 8.5),
            MODULE.probability_for_score(restored, 8.5),
        )

    def test_sdrf_shortlist_filters_probability_q_and_ambiguity(self) -> None:
        base = {
            "pxd_accession": "PXDTEST",
            "data_file": "run1.raw",
            "diagnostic_unimod_candidate_count": 1,
            "ptm_candidate_probability": 0.99,
            "candidate_q_value": 0.01,
            "unimod_accession": "UNIMOD:21",
            "unimod_name": "Phospho",
            "unimod_residual_da": 0.001,
            "pair_support": 100,
            "median_spectral_similarity": 0.9,
        }
        rows = [base, {**base, "data_file": "run2.raw"}]
        # High-probability but ambiguous candidate must not enter the shortlist.
        rows.append(
            {
                **base,
                "data_file": "run3.raw",
                "unimod_accession": "UNIMOD:999",
                "unimod_name": "Ambiguous",
                "diagnostic_unimod_candidate_count": 8,
            }
        )
        suggestions = MODULE.build_sdrf_suggestions(
            rows,
            probability_threshold=0.95,
            q_value_threshold=0.05,
            maximum_candidate_ambiguity=3,
            minimum_confident_runs=2,
            minimum_confident_run_fraction=0.5,
            maximum_candidates_per_accession=5,
        )
        self.assertEqual(len(suggestions), 1)
        self.assertEqual(suggestions[0]["unimod_accession"], "UNIMOD:21")
        self.assertEqual(
            suggestions[0]["sdrf_value"],
            "NT=Phospho;AC=UNIMOD:21",
        )
        self.assertEqual(suggestions[0]["sdrf_status"], "review-required")
        self.assertIn("does not establish searched modification", suggestions[0]["sdrf_warning"])

    def test_read_tables_only_scores_biological_ptm_candidates(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "PXDTEST" / "repository-metadata" / "1_run.raw"
            root.mkdir(parents=True)
            path = root / "run.raw.mass-shifts.tsv"
            fields = [
                "data_file",
                "cluster_rank",
                "delta_mass_da",
                "cluster_sigma_da",
                "cluster_min_da",
                "cluster_max_da",
                "pair_support",
                "unique_spectrum_support",
                "median_spectral_similarity",
                "classification",
                "confidence",
                "match_tolerance_da",
                "support_fraction_of_accepted_pairs",
                "diagnostic_unimod_candidate_count",
                "unimod_accession",
                "unimod_name",
                "unimod_theoretical_delta_mass_da",
                "unimod_residual_da",
                "unimod_source_classification",
                "unimod_candidate_category",
            ]
            rows = [
                [
                    "run.raw", "1", "79.9665", "0.001", "79.965", "79.968",
                    "80", "120", "0.85", "putative-ptm", "high-support", "0.02",
                    "0.01", "2", "UNIMOD:21", "Phospho", "79.966331", "0.000169",
                    "Post-translational", "biological-ptm",
                ],
                [
                    "run.raw", "1", "79.9665", "0.001", "79.965", "79.968",
                    "80", "120", "0.85", "putative-ptm", "high-support", "0.02",
                    "0.01", "2", "UNIMOD:999", "Chemical", "79.9664", "0.0001",
                    "Chemical derivative", "sample-prep-or-artifact",
                ],
            ]
            with path.open("w", encoding="utf-8", newline="") as handle:
                handle.write("\t".join(fields) + "\n")
                for row in rows:
                    handle.write("\t".join(row) + "\n")

            clusters, candidates = MODULE.read_mass_shift_tables(Path(tmp))
            targets = MODULE._biological_targets(candidates)
            self.assertEqual(len(clusters), 1)
            self.assertEqual(len(candidates), 2)
            self.assertEqual(len(targets), 1)
            self.assertEqual(targets[0].unimod_accession, "UNIMOD:21")

    def test_metadata_semantics_explicitly_reject_identity_probability(self) -> None:
        model = MODULE.fit_probability_model(
            [8.0, 7.0],
            [6.0],
            decoy_offsets=(0.3,),
            catalog_hash="x",
            biological_catalog_records=2,
            training_accessions=("PXD1",),
            minimum_targets_per_bin=1,
        )
        payload = json.dumps(MODULE.model_to_json(model))
        self.assertIn("not peptide/site/PTM-identity probability", payload)


if __name__ == "__main__":
    unittest.main()

APPLY_SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "benchmarks"
    / "apply_confident_ptm_sdrf_review.py"
)
APPLY_SPEC = importlib.util.spec_from_file_location(
    "apply_confident_ptm_sdrf_review",
    APPLY_SCRIPT,
)
assert APPLY_SPEC and APPLY_SPEC.loader
APPLY_MODULE = importlib.util.module_from_spec(APPLY_SPEC)
sys.modules[APPLY_SPEC.name] = APPLY_MODULE
APPLY_SPEC.loader.exec_module(APPLY_MODULE)


class ConfidentPtmSdrfApplyTests(unittest.TestCase):
    def test_review_requires_explicit_acceptance(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            review = Path(tmp) / "review.tsv"
            review.write_text(
                "pxd_accession\tunimod_accession\tunimod_name\tsdrf_value\tsdrf_status\n"
                "PXD1\tUNIMOD:21\tPhospho\tNT=Phospho;AC=UNIMOD:21\treview-required\n",
                encoding="utf-8",
            )
            self.assertEqual(APPLY_MODULE.read_review(review), [])

    def test_apply_adds_only_accepted_new_accessions(self) -> None:
        columns = [
            "source name",
            "comment[data file]",
            "comment[modification parameters]",
        ]
        rows = [
            ["s1", "a.raw", "NT=Carbamidomethyl;AC=UNIMOD:4"],
            ["s2", "b.raw", "NT=Carbamidomethyl;AC=UNIMOD:4"],
        ]
        accepted = [
            APPLY_MODULE.AcceptedModification(
                "PXD1",
                "UNIMOD:21",
                "Phospho",
                "NT=Phospho;AC=UNIMOD:21",
            ),
            APPLY_MODULE.AcceptedModification(
                "PXD1",
                "UNIMOD:4",
                "Carbamidomethyl",
                "NT=Carbamidomethyl;AC=UNIMOD:4",
            ),
        ]
        new_columns, new_rows, audit = APPLY_MODULE.apply_review(columns, rows, accepted)
        self.assertEqual(new_columns.count("comment[modification parameters]"), 2)
        self.assertTrue(all(row[-1] == "NT=Phospho;AC=UNIMOD:21" for row in new_rows))
        self.assertEqual([item["action"] for item in audit], [
            "added",
            "skipped-existing-accession",
        ])

    def test_apply_keeps_factor_columns_last(self) -> None:
        columns = [
            "source name",
            "comment[modification parameters]",
            "comment[file uri]",
            "comment[data file]",
            "factor value[compound]",
            "factor value[enrichment process]",
        ]
        rows = [
            [
                "s1",
                "NT=Carbamidomethyl;AC=UNIMOD:4",
                "file:///a.raw",
                "a.raw",
                "none",
                "phosphoproteomics",
            ]
        ]
        accepted = [
            APPLY_MODULE.AcceptedModification(
                "PXD1",
                "UNIMOD:21",
                "Phosphorylation",
                "NT=Phosphorylation;AC=UNIMOD:21",
            )
        ]

        new_columns, new_rows, audit = APPLY_MODULE.apply_review(columns, rows, accepted)

        self.assertEqual(
            new_columns,
            [
                "source name",
                "comment[modification parameters]",
                "comment[modification parameters]",
                "comment[file uri]",
                "comment[data file]",
                "factor value[compound]",
                "factor value[enrichment process]",
            ],
        )
        self.assertEqual(new_rows[0][2], "NT=Phosphorylation;AC=UNIMOD:21")
        self.assertEqual(new_rows[0][-2:], ["none", "phosphoproteomics"])
        self.assertEqual(audit[0]["action"], "added")

    def test_accepted_review_value_must_contain_matching_accession(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            review = Path(tmp) / "review.tsv"
            review.write_text(
                "pxd_accession\tunimod_accession\tunimod_name\tsdrf_value\tsdrf_status\n"
                "PXD1\tUNIMOD:21\tPhospho\tNT=Phospho;AC=UNIMOD:35\taccepted\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "does not contain"):
                APPLY_MODULE.read_review(review)
