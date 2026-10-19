"""PSI-MS DIA window metrics: actual geometry, acquisition gating, and fail-closed cases."""

from __future__ import annotations

import unittest
from pathlib import Path

import numpy as np

from prideqc.metrics import QCMetricCalculator, RunSummary
from prideqc.models import CVTerm, Metric, Precursor, RunMetadata, Spectrum
from prideqc.mzqc import definition
from prideqc.pipeline import Analyzer
from tests.helpers import MemoryReader


def scan(rt, target=425.0, lower=10.0, upper=10.0, *, cv=None,
         mobility=None, level=2, extra=False, unresolved=False):
    first = Precursor(target, 2, 1.0, lower + upper if lower is not None and upper is not None else None,
                      isolation_lower_offset=lower, isolation_upper_offset=upper,
                      ion_mobility_lower=mobility[0] if mobility else None,
                      ion_mobility_upper=mobility[1] if mobility else None)
    return Spectrum(level, float(rt), np.array([100.]), np.array([1.]),
                    precursors=(first, Precursor(450.0)) if extra else (first,),
                    faims_cv=cv, unresolved_ion_mobility=unresolved)


def values(items, *, dia=True):
    run = RunSummary()
    for item in items:
        run.consume_spectrum(item)
    return {m.key: m.value for m in QCMetricCalculator().calculate(run, dia=dia)}


class DIAIsolationWindowTests(unittest.TestCase):
    def test_same_target_different_offsets_are_two_windows(self):
        got = values([scan(1, lower=5, upper=5), scan(2, lower=10, upper=10),
                      scan(3, lower=5, upper=5)])
        self.assertEqual(got["DIAIsolationWindowCount"], 2)
        self.assertEqual(got["DIAIsolationWindowWidths"], [10.0, 20.0])
        self.assertEqual(definition(Metric("DIAIsolationWindowCount", 2)).accession, "MS:4000194")
        self.assertEqual(definition(Metric("DIAIsolationWindowWidths", [10, 20])).accession,
                         "MS:4000195")
        self.assertEqual(definition(Metric("DIAIsolationWindowWidths", [10, 20])).shape, "tuple")

    def test_fixed_width_repeats_both_endpoints(self):
        got = values([scan(1, target=400), scan(2, target=500), scan(3, target=400)])
        self.assertEqual(got["DIAIsolationWindowCount"], 2)
        self.assertEqual(got["DIAIsolationWindowWidths"], [20.0, 20.0])

    def test_no_dia_proof_does_not_publish_dia_metrics(self):
        self.assertNotIn("DIAIsolationWindowCount", values([scan(1)], dia=False))
        self.assertNotIn("DIAIsolationWindowWidths", values([scan(1)], dia=False))

    def test_faims_and_mobility_are_distinct_window_dimensions(self):
        got = values([scan(1, cv=-45, mobility=(0.8, 1.0)),
                      scan(2, cv=-50, mobility=(0.8, 1.0)),
                      scan(3, cv=-45, mobility=(0.9, 1.1))])
        self.assertEqual(got["DIAIsolationWindowCount"], 3)

    def test_mixed_missing_faims_or_im_dimensions_abstains(self):
        self.assertNotIn("DIAIsolationWindowCount", values([
            scan(1, cv=-45), scan(2, cv=None)]))
        self.assertNotIn("DIAIsolationWindowCount", values([
            scan(1, mobility=(0.8, 1.0)), scan(2)]))

    def test_unresolved_mobility_and_incomplete_window_abstain(self):
        for items in (
            [scan(1, unresolved=True)],
            [scan(1, lower=None, upper=None)],
            [scan(1, extra=True)],
            [scan(1, lower=-1, upper=11)],
            [scan(1, mobility=(1.2, 0.8))],
            [scan(1), scan(2, lower=None, upper=None)],
        ):
            with self.subTest(items=items):
                self.assertNotIn("DIAIsolationWindowCount", values(items))
                self.assertNotIn("DIAIsolationWindowWidths", values(items))

    def test_rounding_noise_in_window_boundaries_is_suppressed(self):
        got = values([scan(1), scan(2, target=425.00000002)])
        self.assertEqual(got["DIAIsolationWindowCount"], 1)

    def test_analyzer_gates_using_actual_acquisition_annotation(self):
        observed = RunMetadata(acquisition_methods=[
            CVTerm("PRIDE:0000450", "Data-independent acquisition")])
        result = Analyzer(reader=MemoryReader([scan(1), scan(2)], observed)).analyze(Path("test.mzML"))
        self.assertIn("DIAIsolationWindowCount", {m.key for m in result.metrics})
        other = RunMetadata(acquisition_methods=[
            CVTerm("PRIDE:0000627", "Data-dependent acquisition")])
        result = Analyzer(reader=MemoryReader([scan(1), scan(2)], other)).analyze(Path("test.mzML"))
        self.assertNotIn("DIAIsolationWindowCount", {m.key for m in result.metrics})

    def test_all_retired_metrics_no_longer_serializable(self):
        removed = (
            "IsolationWidth_MS2_Median", "IsolationWidth_MS2_Quantiles",
            "IsolationWidth_MS2_Min", "IsolationWidth_MS2_Max", "IsolationWidth_MS2_Count",
            "IsolationPrecursorMz_MS2_UniqueCount", "IsolationPrecursorMz_MS2_UniqueFraction",
            "IsolationPrecursorMz_MS2_MaxRepeatFraction", "IsolationPrecursorMz_MS2_RepeatedFraction",
            "MultiplePrecursors_MS2_Count", "MS2_PrecursorCharge_Fractions",
        )
        for key in removed:
            with self.subTest(key=key), self.assertRaisesRegex(ValueError, "Unregistered QC metric"):
                definition(Metric(key, 1))
