from __future__ import annotations

import importlib.util
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
        self.assertEqual(MODULE.parse_historical_tolerance("NT=5 ppm;AC=foo"), (5.0, "ppm"))
        self.assertEqual(MODULE.parse_historical_tolerance("27 % NCE"), (None, ""))


if __name__ == "__main__":
    unittest.main()
