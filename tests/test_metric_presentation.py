"""Editorial QC tier registry must never redefine measurements or their CV terms."""

from __future__ import annotations

import csv
import tempfile
import unittest
from collections import Counter
from pathlib import Path

from scripts.generate_qc_metrics_reference import (
    CORE_COUNT,
    DEST,
    PRESENTATION,
    REGISTRY,
    load_catalog,
    render,
)


class MetricPresentationTests(unittest.TestCase):
    def test_exact_registry_bijection_and_counts(self):
        catalog = load_catalog()
        self.assertEqual(len(catalog), 104)
        self.assertEqual(
            Counter(policy["tier"] for _, policy in catalog),
            {"core": 18, "diagnostic": 39, "extended": 47},
        )
        self.assertEqual(CORE_COUNT, 18)
        self.assertEqual({metric["shape"] for metric, _ in catalog},
                         {"scalar", "tuple", "table"})

    def test_core_covers_essential_run_measurements(self):
        policies = {m["key"]: p for m, p in load_catalog()}
        for key in (
            "NumberOfSpectra", "NumberOfSpectra_MS1", "NumberOfSpectra_MS2",
            "ChromatographyDuration", "MS2_PrecursorCharge_Fractions",
            "TIC_MS1_Median", "TIC_MS1_CV", "TIC_MS2_Median",
            "IsolationWidth_MS2_Median", "AcquisitionCycle_MS2Count_Quantiles",
        ):
            with self.subTest(key=key):
                self.assertEqual(policies[key]["tier"], "core")
        self.assertEqual(policies["AcquisitionCycle_MS2Count_Quantiles"]["applicability"],
                         "cycles")
        self.assertEqual(policies["MS2_PrecursorCharge_Fractions"]["applicability"],
                         "ms2")
        self.assertEqual(policies["FAIMS_CV_Values"]["tier"], "diagnostic")
        self.assertEqual(policies["FAIMS_CV_Values"]["applicability"], "faims")

    def test_display_policy_cannot_change_source_accessions_or_definitions(self):
        with REGISTRY.open(encoding="utf-8", newline="") as handle:
            raw = list(csv.DictReader(handle, delimiter="\t"))
        catalog = load_catalog()
        self.assertEqual(raw, [metric for metric, _ in catalog])
        # All 104 definitions, including extended entries, remain in the reference.
        doc = render()
        for metric in raw:
            self.assertIn(f"`{metric['accession']}`", doc)
        self.assertIn("MS2/MS3", doc)
        self.assertIn("not comparable across uncalibrated runs", doc)

    def test_generated_document_is_current(self):
        self.assertEqual(DEST.read_text(encoding="utf-8"), render())
        self.assertIn("**18 core**", render())
        self.assertIn("**39 diagnostic**", render())
        self.assertIn("**47 extended**", render())

    def test_duplicate_or_missing_policy_fails_closed(self):
        with tempfile.TemporaryDirectory() as temp:
            policy_path = Path(temp) / "policy.tsv"
            policy_path.write_text(PRESENTATION.read_text(encoding="utf-8"))
            text = policy_path.read_text(encoding="utf-8")
            lines = text.splitlines(keepends=True)
            policy_path.write_text("".join(lines[:-1]), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "match the registry exactly"):
                load_catalog(presentation=policy_path)
            policy_path.write_text(text + lines[1], encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "Duplicate metric keys"):
                load_catalog(presentation=policy_path)

    def test_invalid_tier_and_missing_core_caveat_fail_closed(self):
        with tempfile.TemporaryDirectory() as temp:
            policy_path = Path(temp) / "policy.tsv"
            contents = PRESENTATION.read_text(encoding="utf-8")
            policy_path.write_text(contents.replace("\tcore\t", "\tunreviewed\t", 1),
                                   encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "Unrecognized presentation tier"):
                load_catalog(presentation=policy_path)
            rows = list(csv.DictReader(contents.splitlines(), delimiter="\t"))
            first_core = next(row for row in rows if row["tier"] == "core")
            first_core["caveat"] = ""
            with policy_path.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=rows[0].keys(), delimiter="\t")
                writer.writeheader()
                writer.writerows(rows)
            with self.assertRaisesRegex(ValueError, "Core metric lacks"):
                load_catalog(presentation=policy_path)


if __name__ == "__main__":
    unittest.main()
