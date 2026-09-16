from __future__ import annotations

import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

MODULE_PATH = (
    Path(__file__).resolve().parents[1]
    / "scripts"
    / "benchmarks"
    / "compare_historical_tolerances.py"
)
SPEC = importlib.util.spec_from_file_location("compare_historical_tolerances", MODULE_PATH)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class HistoricalToleranceComparisonTests(unittest.TestCase):
    def test_candidate_field_names_and_historical_tolerance(self):
        annotations = {
            "suggested_fragment_search_tolerance_da": {
                "value": '{"suggested_tolerance": 0.4}',
                "evidence": "inferred",
            }
        }
        unit, value, payload = MODULE.candidate_from_annotations(annotations)
        self.assertEqual(unit, "Da")
        self.assertEqual(value, 0.4)
        self.assertEqual(payload["suggested_tolerance"], 0.4)

    def test_parse_historical_tolerance(self):
        self.assertEqual(MODULE.parse_historical_tolerance("20 ppm"), (20.0, "ppm"))
        self.assertEqual(MODULE.parse_historical_tolerance("0.8 Da"), (0.8, "Da"))
        self.assertEqual(MODULE.parse_historical_tolerance("not available"), (None, ""))
        self.assertEqual(MODULE.parse_historical_tolerance("27 % NCE"), (None, ""))

    def test_mass_context_conversion_da_to_ppm_is_not_primary_comparison(self):
        result = MODULE.mass_context_conversion(
            0.05,
            "Da",
            10.0,
            "ppm",
            (400.0, 500.0, 600.0),
        )
        self.assertEqual(result["mass_context_unit"], "ppm")
        self.assertEqual(float(result["mass_context_mz_min"]), 400.0)
        self.assertEqual(float(result["mass_context_mz_midpoint"]), 500.0)
        self.assertEqual(float(result["mass_context_mz_max"]), 600.0)
        self.assertEqual(float(result["historical_equivalent_min"]), 125.0)
        self.assertEqual(float(result["historical_equivalent_midpoint"]), 100.0)
        self.assertAlmostEqual(float(result["historical_equivalent_max"]), 83.3333333333, places=6)
        self.assertEqual(float(result["candidate_equivalent_midpoint"]), 10.0)
        self.assertEqual(
            float(result["mass_context_ratio_midpoint"]),
            0.1,
        )

    def test_mass_context_conversion_requires_positive_observed_range(self):
        result = MODULE.observed_ms2_mass_context(
            {"ObservedMzRange_MS2": json.dumps([400.0, 600.0])}
        )
        self.assertEqual(result, (400.0, 500.0, 600.0))
        self.assertEqual(
            MODULE.observed_ms2_mass_context({"ObservedMzRange_MS2": "not-json"}),
            (None, None, None),
        )

    def test_same_unit_ratio_uses_matching_units_only(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            result_dir = root / "result"
            result_dir.mkdir()
            (result_dir / "hpc-run-info.txt").write_text(
                "task_id=1\npxd=PXDTEST\narchive_file=test.raw\ntask_outcome=success\n",
                encoding="utf-8",
            )
            annotation_rows = [
                {
                    "data_file": "test.raw",
                    "field": "suggested_fragment_search_tolerance_ppm",
                    "value": json.dumps({"suggested_tolerance": 12.0}),
                    "evidence": "inferred",
                },
                {
                    "data_file": "test.raw",
                    "field": "estimated_fragment_mass_error_ppm",
                    "value": json.dumps({"robust_inlier_fraction": 0.9}),
                    "evidence": "inferred",
                },
            ]
            with (result_dir / "annotations.tsv").open("w", encoding="utf-8") as handle:
                handle.write("data_file\tfield\tvalue\tevidence\n")
                for row in annotation_rows:
                    handle.write(
                        "\t".join(
                            row[key]
                            for key in ("data_file", "field", "value", "evidence")
                        )
                        + "\n"
                    )
            with (result_dir / "metrics.tsv").open("w", encoding="utf-8") as handle:
                handle.write("data_file\tObservedMzRange_MS2\t[400,600]\n")
                handle.write("test.raw\tObservedMzRange_MS2\t[400,600]\n")
            sdrf = result_dir / "input.sdrf.tsv"
            sdrf.write_text(
                "comment[data file]\tcomment[fragment mass tolerance]\n"
                "test.raw\t20 ppm\n",
                encoding="utf-8",
            )
            row = MODULE.build_row(result_dir, {})
            self.assertEqual(row["comparison"], "same-unit")
            self.assertAlmostEqual(float(row["ratio_candidate_over_historical"]), 0.6)
            self.assertEqual(row["mass_context_ratio_midpoint"], "0.6")
            self.assertEqual(row["mass_context_ratio_min"], "0.6")
            self.assertEqual(row["mass_context_ratio_max"], "0.6")


if __name__ == "__main__":
    unittest.main()
