"""Numerical oracles for the shared metrics, including upstream edge cases."""

import math
import unittest

import numpy as np

from prideqc.metrics import QCMetricCalculator, RunSummary
from prideqc.models import Precursor, Spectrum
from tests.helpers import spectrum


def calculate(spectra):
    summary = RunSummary()
    for item in spectra:
        summary.consume_spectrum(item)
    return {m.key: m.value for m in QCMetricCalculator().calculate(summary)}, summary


class MetricTests(unittest.TestCase):
    def test_id_free_middle_and_shortest_half_proxies(self):
        items = [spectrum(rt, [tic]) for rt, tic in zip(
            [0, 1, 2, 10, 20], [1, 100, 3, 4, 5], strict=True
        )]
        metrics, _ = calculate(items)
        self.assertEqual(metrics["MedianTIC_in_RT_MS1_IQR"], 3.5)
        self.assertEqual(metrics["TIC_MS1_MedianInHalfRange"], 3)
        self.assertAlmostEqual(metrics["RT_MS1_IQRRate"], 1 / 3)

    def test_precursor_extent_uses_only_recorded_intensity(self):
        metrics, _ = calculate([
            spectrum(0, [999], 2, charge=2, precursor_intensity=0),
            spectrum(1, [999], 2, charge=2, precursor_intensity=10),
            spectrum(2, [999], 2, charge=2, precursor_intensity=30),
        ])
        self.assertAlmostEqual(metrics["ExtentPrecursorIntensity_95over5_MS2"], 29 / 11)

    def test_irregular_tic_integral_and_exact_quartile_boundaries(self):
        metrics, _ = calculate([spectrum(0, [0]), spectrum(10, [20])])
        self.assertEqual(metrics["TIC_MS1_Area"], 100)
        np.testing.assert_allclose(metrics["TIC_MS1_Area_RTQuantiles"], [6.25, 18.75, 31.25, 43.75])
        self.assertEqual(sum(metrics["TIC_MS1_Area_RTQuantiles"]), 100)

    def test_irregular_spacing_and_unsorted_input(self):
        metrics, _ = calculate([spectrum(7, [5]), spectrum(0, [1]), spectrum(2, [3])])
        self.assertEqual(metrics["TIC_MS1_Area"], 24)
        self.assertAlmostEqual(sum(metrics["TIC_MS1_Area_RTQuantiles"]), 24)

    def test_level_specific_rate_and_all_level_duration(self):
        metrics, _ = calculate([spectrum(0, [1]), spectrum(60, [1]),
                                spectrum(100, [1], 2), spectrum(110, [1], 2),
                                spectrum(200, [999], 3)])
        self.assertEqual(metrics["ScanRate_MS1"], 2)
        self.assertEqual(metrics["ScanRate_MS2"], 12)
        self.assertEqual(metrics["ChromatographyDuration"], 200)
        self.assertEqual(metrics["BasePeak_All_Max"], 999)
        self.assertEqual(metrics["NumberOfSpectra_MS3"], 1)

    def test_finite_rt_filtering_does_not_drop_scan_counts(self):
        metrics, summary = calculate([spectrum(0, [1]), spectrum(10, [1]), spectrum(math.nan, [3])])
        self.assertEqual(metrics["NumberOfSpectra_MS1"], 3)
        self.assertEqual(metrics["ScanRate_MS1"], 12)
        self.assertEqual(metrics["TIC_MS1_Area"], 10)
        self.assertTrue(summary.warnings())

    def test_empty_scans_count_toward_density(self):
        metrics, _ = calculate([spectrum(0, []), spectrum(1, [1, 1]), spectrum(2, [0, 0])])
        self.assertEqual(metrics["EmptyScans_MS1"], 2)
        self.assertEqual(metrics["PeakDensity_MS1_Quantiles"], [1, 2, 2])

    def test_invalid_peaks_are_excluded_and_reported(self):
        metrics, summary = calculate(
            [spectrum(0, [1, math.nan, -1, 5], mz=[100, 200, 300, math.inf])],
        )
        self.assertEqual(metrics["InvalidPeakCount_MS1"], 3)
        self.assertEqual(metrics["NumberOfSpectralPeaks"], 1)
        self.assertEqual(metrics["TIC_MS1_Median"], 1)
        self.assertTrue(summary.warnings())

    def test_fastest_frequency_uses_closed_minute_not_minimum_gap(self):
        metrics, _ = calculate([spectrum(t, [1]) for t in [0, 0.001, 60, 120]])
        self.assertEqual(metrics["FastestFrequency_MS1"], 3 / 60)

    def test_duplicate_rt_windows_include_all_scans(self):
        metrics, _ = calculate([spectrum(t, [1]) for t in [0, 0, 60]])
        self.assertEqual(metrics["FastestFrequency_MS1"], 3 / 60)
        self.assertAlmostEqual(sum(metrics["TIC_MS1_Area_RTQuantiles"]), 60)

    def test_sample_cv(self):
        metrics, _ = calculate([spectrum(0, [1]), spectrum(1, [3])])
        self.assertAlmostEqual(metrics["TIC_MS1_CV"], math.sqrt(2) / 2)

    def test_tenfold_changes_zero_denominator_and_exact_boundary(self):
        metrics, _ = calculate([spectrum(i, [x]) for i, x in enumerate([0, 1, 10, 1, 0])])
        self.assertEqual(metrics["TIC_MS1_SignalJump10x_Count"], 1)
        self.assertEqual(metrics["TIC_MS1_SignalFall10x_Count"], 2)

    def test_absolute_tic_changes_produce_finite_ratios(self):
        metrics, _ = calculate([spectrum(i, [x]) for i, x in enumerate([1, 9, 3, 7, 5])])
        np.testing.assert_allclose(
            metrics["TIC_MS1_ChangeQuartileLogRatios"],
            np.log([5 / 3.5, 6.5 / 5, 8 / 6.5]),
        )

    def test_rt_interval_fractions_sum_to_one(self):
        metrics, _ = calculate([spectrum(t, [1]) for t in [0, 1, 4, 10, 20]])
        self.assertAlmostEqual(sum(metrics["RT_MS1_Quantiles"]), 1)
        self.assertAlmostEqual(sum(metrics["RT_TIC_Quantiles"]), 1)

    def test_charge_denominator_includes_unknown_and_negative(self):
        metrics, _ = calculate(
            [spectrum(i, [1], 2, charge=c) for i, c in enumerate([2, 3, 0, -1, None, 7])],
        )
        table = metrics["MS2_PrecursorCharge_Fractions"]
        self.assertEqual(table["count"], [0, 1, 1, 0, 0, 1, 3])
        self.assertAlmostEqual(sum(table["fraction"]), 1)
        self.assertEqual(metrics["ChargeMin"], 2)
        self.assertEqual(metrics["ChargeMean"], 4)

    def test_missing_precursor_intensity_is_never_fragment_tic(self):
        metrics, _ = calculate([spectrum(1, [999], 2, charge=2, precursor_intensity=0)])
        self.assertIsNone(metrics["PrecursorIntensity_MS2_Mean"])
        self.assertEqual(metrics["PrecursorIntensity_MS2_MissingCount"], 1)

    def test_multiple_precursor_rule_is_visible(self):
        item = Spectrum(2, 0, np.array([100.]), np.array([1.]),
                        (Precursor(400, 2, 10), Precursor(600, 4, 20)))
        metrics, summary = calculate([item])
        self.assertEqual(metrics["ChargeMean"], 2)
        self.assertEqual(metrics["MzRange_MS2"], [400, 400])
        self.assertTrue(summary.warnings())

    def test_no_spectra_and_single_scan_undefined_values(self):
        metrics, _ = calculate([])
        self.assertEqual(metrics["NumberOfSpectra_MS1"], 0)
        self.assertIsNone(metrics["ChromatographyDuration"])
        metrics, _ = calculate([spectrum(0, [1])])
        self.assertIsNone(metrics["ScanRate_MS1"])
        self.assertIsNone(metrics["TIC_MS1_Area"])

    def test_chromatograms_are_counts_of_points_not_features(self):
        summary = RunSummary()
        summary.consume_chromatogram(np.array([math.nan, 0, 10]), 1)
        metrics = {m.key: m.value for m in QCMetricCalculator().calculate(summary)}
        self.assertEqual(metrics["NumberOfChromatogramDataPoints"], 3)
        self.assertEqual(metrics["ChromatogramRTRange"], [0, 10])

    def test_faims_is_not_drift_time(self):
        item = Spectrum(1, 0, np.array([]), np.array([]), faims_cv=-45)
        metrics, _ = calculate([item])
        self.assertEqual(metrics["FAIMS_CV_Values"], [-45])

    def test_summary_does_not_retain_peak_arrays(self):
        import gc
        import weakref

        summary = RunSummary()
        item = spectrum(0, [1] * 10000)
        reference = weakref.ref(item.intensity)
        summary.consume_spectrum(item)
        del item
        gc.collect()
        self.assertIsNone(reference())
        self.assertEqual(len(summary.levels[1].tic), 1)
