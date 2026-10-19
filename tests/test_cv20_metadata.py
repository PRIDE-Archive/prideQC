"""Official PSI-MS value shapes, FAIMS metadata, and retired public metrics."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from prideqc.metrics import QCMetricCalculator, RunSummary
from prideqc.models import AnalysisResult, Metric
from prideqc.mzqc import MzQCWriter, definition
from prideqc.pipeline import Analyzer
from tests.helpers import MemoryReader, spectrum


class Cv20Tests(unittest.TestCase):
    def test_existing_qc_terms_match_value_contract(self):
        for key, cv, shape in (
            ("ChromatographyDuration", "MS:4000067", "scalar"),
            ("ObservedMzRange_MS1", "MS:4000222", "tuple"),
            ("MedianCycleTime_MS1", "MS:4000192", "scalar"),
        ):
            with self.subTest(key=key):
                term = definition(Metric(key, 1))
                self.assertEqual((term.accession, term.shape), (cv, shape))

    def test_faims_are_scalar_cv_metadata_and_roundtrip(self):
        from dataclasses import replace

        examples = [replace(spectrum(t, [1]), faims_cv=v)
                    for t, v in [(0, -45.0), (1, -65.0), (2, -45.0)]]
        result = Analyzer(MemoryReader(examples)).analyze("faims.mzML")
        self.assertEqual(result.metadata.faims_compensation_voltages, [-65.0, -45.0])
        self.assertEqual(AnalysisResult.from_dict(result.to_dict()).metadata.faims_compensation_voltages,
                         [-65.0, -45.0])
        with tempfile.TemporaryDirectory() as temp:
            doc = MzQCWriter().build(result, Path(temp) / "local.obo")
        run = doc["mzQC"]["runQualities"][0]
        props = run["metadata"]["inputFiles"][0]["fileProperties"]
        faims = [p for p in props if p["accession"] == "MS:1001581"]
        self.assertEqual([p["value"] for p in faims], [-65.0, -45.0])
        self.assertTrue(all(p["unit"]["accession"] == "UO:0000218" for p in faims))
        self.assertFalse(any(p["accession"] in {"MS:1003450", "MS:1003451"} for p in props))
        qc_keys = {m.key for m in result.metrics}
        self.assertFalse(any(key.startswith(("FAIMS_CV_", "AcquisitionCycle_")) for key in qc_keys))
        self.assertNotIn("NumberOfMSLevels", qc_keys)

    def test_per_spectrum_ms_level_must_not_be_misused_as_count(self):
        for key in ("NumberOfMSLevels", "FAIMS_CV_Values", "AcquisitionCycle_Count"):
            with self.subTest(key=key), self.assertRaises(ValueError):
                definition(Metric(key, 2))
        active = {metric.key for metric in QCMetricCalculator().calculate(RunSummary())}
        self.assertNotIn("NumberOfMSLevels", active)


if __name__ == "__main__":
    unittest.main()
