"""MS-level coverage: absent levels are explicit and never read as measured zeros."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from prideqc.mzqc import MzQCWriter
from prideqc.pipeline import Analyzer
from tests.helpers import MemoryReader, spectrum


def _coverage(result):
    return next(a for a in result.annotations if a.field == "ms_level_coverage")


class MsLevelCoverageTests(unittest.TestCase):
    def test_ms2_only_run_marks_ms1_metrics_not_applicable_in_mzqc(self):
        spectra = [spectrum(i, [10, 20], 2, charge=2) for i in range(5)]
        result = Analyzer(MemoryReader(spectra)).analyze("prm.raw")
        coverage = _coverage(result)
        self.assertEqual(coverage.value["status"], "MS1_ABSENT")
        self.assertEqual(coverage.value["absent_levels"], ["MS1"])
        self.assertEqual(coverage.value["spectra_per_level"], {"MS2": 5})
        self.assertIn("RT_TIC_Quantiles", coverage.value["not_applicable_metrics"])
        self.assertIn("MS1_to_MS2_Ratio", coverage.value["not_applicable_metrics"])
        self.assertNotIn("ChargeMean", coverage.value["not_applicable_metrics"])
        self.assertIn("MS1: no spectra acquired; MS1 metrics are not applicable.", result.warnings)
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "prm.mzQC"
            MzQCWriter().write(result, out)
            run = MzQCWriter().build(result, out.with_suffix(".obo"))["mzQC"]["runQualities"][0]
            named = {m["name"]: m["value"] for m in run["qualityMetrics"]}
            self.assertEqual(named["MS level coverage and not-applicable metrics"]["status"],
                             "MS1_ABSENT")
            self.assertEqual(named["number of MS1 spectra"], 0)
            self.assertNotIn("MS1 median cycle time", named)
            self.assertIn("MS level coverage", out.with_suffix(".obo").read_text())

    def test_full_run_has_no_not_applicable_metrics(self):
        spectra = [spectrum(0, [10]), spectrum(1, [10], 2, charge=2)]
        coverage = _coverage(Analyzer(MemoryReader(spectra)).analyze("dda.mzML"))
        self.assertEqual(coverage.value["status"], "MS1_AND_MS2_PRESENT")
        self.assertEqual(coverage.value["not_applicable_metrics"], [])

    def test_empty_run_is_no_spectra(self):
        result = Analyzer(MemoryReader([])).analyze("empty.mzML")
        self.assertEqual(_coverage(result).value["status"], "NO_SPECTRA")
        self.assertNotIn("MS1: no spectra acquired; MS1 metrics are not applicable.",
                         result.warnings)


if __name__ == "__main__":
    unittest.main()
