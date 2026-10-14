"""Offline tests for isolated experimental profile TMT reporter evidence."""

from __future__ import annotations

import unittest
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from prideqc.models import EvidenceKind
from prideqc.mzqc import annotation_definition
from prideqc.pipeline import Analyzer, WorkflowOptions
from prideqc.tmt import TMTReporterCollector
from prideqc.tmt_profile import TMTProfileReporterCollector
from tests.helpers import MemoryReader, spectrum


class _MockTMT6:
    def getChannelInformation(self):
        return [
            SimpleNamespace(name=name, center=mz)
            for name, mz in (
                ("126", 126.127725),
                ("127", 127.124760),
                ("128", 128.134433),
                ("129", 129.131468),
            )
        ]


class _MockSpectrum:
    def __init__(self):
        self.peaks = (np.array([]), np.array([]))

    def set_peaks(self, peaks):
        self.peaks = peaks

    def get_peaks(self):
        return self.peaks


class _MockPicker:
    """Deterministic fixture picker, not a claim about real OpenMS performance."""

    def pick(self, source, target):
        mz, intensity = source.get_peaks()
        maxima = (
            np.flatnonzero(
                (intensity[1:-1] > intensity[:-2])
                & (intensity[1:-1] >= intensity[2:])
                & (intensity[1:-1] > 0)
            )
            + 1
        )
        target.set_peaks((mz[maxima], intensity[maxima]))


class _FailPicker:
    def pick(self, source, target):
        raise RuntimeError("mock peak-picker failure")


class _MockOMS:
    __version__ = "mock-profile-oms"
    TMTSixPlexQuantitationMethod = _MockTMT6
    MSSpectrum = _MockSpectrum
    PeakPickerHiRes = _MockPicker


def _profile(peaks: list[float], level: int = 3):
    """Five points per Gaussian to exercise scan accounting and peak picking."""
    grid = np.linspace(125.5, 130.0, 10001)
    y = np.zeros_like(grid)
    for center in peaks:
        y += 1000 * np.exp(-0.5 * ((grid - center) / 0.0003) ** 2)
    return spectrum(1, y.tolist(), level, mz=grid.tolist(), representation="profile")


class TMTProfileTests(unittest.TestCase):
    def test_three_reporter_peaks_are_profile_only(self):
        collector = TMTProfileReporterCollector(oms=_MockOMS())
        collector.consume_spectrum(_profile([126.127725, 127.124760, 128.134433]))
        report = collector.annotations()[0]
        self.assertEqual(report.field, "tmt_profile_reporter_evidence")
        self.assertEqual(report.kind, EvidenceKind.INFERRED)
        self.assertIsNone(report.sdrf_column)
        self.assertIsNone(report.sdrf_value)
        level = report.value["levels"]["MS3"]
        self.assertEqual(level["profile_spectra_evaluable"], 1)
        self.assertEqual(level["methods"]["TMT6"]["reporter_positive_spectra"], 1)
        self.assertEqual(level["methods"]["TMT6"]["positive_fraction_of_evaluable_profile"], 1)
        channel = level["methods"]["TMT6"]["channels"][0]
        self.assertIn("profile_spectra_with_local_mz_sampling", channel)
        self.assertNotIn("fraction_of_centroid_spectra", channel)
        self.assertEqual(report.total, 1)
        self.assertIsNotNone(annotation_definition(report))

    def test_centroid_and_profile_outputs_remain_independent_and_one_pass(self):
        profile = _profile([126.127725, 127.124760, 128.134433])
        centroid = spectrum(
            2,
            [1000, 900, 800],
            3,
            mz=[126.127725, 127.124760, 128.134433],
        )
        native = TMTReporterCollector(oms=_MockOMS())
        picked = TMTProfileReporterCollector(oms=_MockOMS())
        reader = MemoryReader([profile, centroid])
        result = Analyzer(reader).analyze(Path("test.mzML"), collectors=[native, picked])
        self.assertEqual(reader.read_count, 1)
        native_evidence = next(a for a in result.annotations if a.field == "tmt_reporter_evidence")
        profile_evidence = next(
            a for a in result.annotations if a.field == "tmt_profile_reporter_evidence"
        )
        self.assertEqual(native_evidence.value["levels"]["MS3"]["profile_skipped"], 1)
        self.assertEqual(native_evidence.value["levels"]["MS3"]["centroid_spectra"], 1)
        self.assertEqual(profile_evidence.value["levels"]["MS3"]["profile_spectra_evaluable"], 1)
        self.assertEqual(profile_evidence.value["levels"]["MS3"]["non_profile_spectra_excluded"], 1)
        from prideqc.mzqc import MzQCWriter

        doc = MzQCWriter().build(result, Path("/tmp/test-tmt-profile-cv.obo"))
        values = [
            item["value"]
            for item in doc["mzQC"]["runQualities"][0]["qualityMetrics"]
            if "tmt_profile_reporter_evidence" in item.get("description", "")
        ]
        self.assertEqual(values, [profile_evidence.value])

    def test_no_reporter_region_and_unknown_are_not_negatives(self):
        picker = TMTProfileReporterCollector(oms=_MockOMS())
        picker.consume_spectrum(
            spectrum(0, [5, 2, 1], 2, mz=[200, 201, 202], representation="profile")
        )
        picker.consume_spectrum(
            spectrum(0, [1, 2, 3], 2, mz=[126, 127, 128], representation="unknown")
        )
        report = picker.annotations()[0]
        level = report.value["levels"]["MS2"]
        self.assertEqual(report.kind, EvidenceKind.UNAVAILABLE)
        self.assertEqual(level["profile_spectra_without_reporter_region"], 1)
        self.assertEqual(level["profile_spectra_evaluable"], 0)
        self.assertIsNone(level["methods"]["TMT6"]["positive_fraction_of_evaluable_profile"])

    def test_flat_profile_has_zero_reporter_support(self):
        collector = TMTProfileReporterCollector(oms=_MockOMS())
        scan = spectrum(
            0,
            [10] * 12001,
            3,
            mz=np.linspace(125.5, 136.5, 12001).tolist(),
            representation="profile",
        )
        collector.consume_spectrum(scan)
        report = collector.annotations()[0]
        self.assertEqual(report.kind, EvidenceKind.OBSERVED)
        self.assertEqual(report.value["levels"]["MS3"]["profile_spectra_with_no_picked_peaks"], 1)
        self.assertEqual(
            report.value["levels"]["MS3"]["methods"]["TMT6"]["reporter_positive_spectra"],
            0,
        )

    def test_failures_and_invalid_arrays_are_distinct(self):
        class FailingOMS(_MockOMS):
            PeakPickerHiRes = _FailPicker

        collector = TMTProfileReporterCollector(oms=FailingOMS())
        collector.consume_spectrum(_profile([126.127725]))
        collector.consume_spectrum(
            spectrum(0, [1, 2, 3], 2, mz=[128, 127, 126], representation="profile")
        )
        levels = collector.annotations()[0].value["levels"]
        self.assertEqual(levels["MS3"]["profile_peak_picker_failures"], 1)
        self.assertEqual(levels["MS2"]["profile_spectra_invalid_samples"], 1)
        self.assertEqual(collector.annotations()[0].kind, EvidenceKind.UNAVAILABLE)

    def test_profile_threshold_is_relative_to_original_full_scan(self):
        collector = TMTProfileReporterCollector(oms=_MockOMS())
        peaks = _profile([126.127725, 127.124760, 128.134433])
        # A very strong non-reporter peak should suppress weak reporter-region
        # matches by the same 0.1%-of-full-scan rule as native centroid data.
        mz = np.append(peaks.mz, [200.0])
        intensity = np.append(peaks.intensity, [2e6])
        collector.consume_spectrum(
            spectrum(0, intensity.tolist(), 2, mz=mz.tolist(), representation="profile")
        )
        report = collector.annotations()[0]
        self.assertEqual(report.kind, EvidenceKind.OBSERVED)
        self.assertEqual(
            report.value["levels"]["MS2"]["methods"]["TMT6"]["reporter_positive_spectra"],
            0,
        )

    def test_profile_cli_and_workflow_are_opt_in(self):
        from prideqc.cli import parser

        self.assertFalse(WorkflowOptions().tmt_profile_reporter_evidence)
        args = parser().parse_args(
            ["analyze", "run.mzML", "-o", "test-output", "--tmt-profile-reporter-evidence"]
        )
        self.assertTrue(args.tmt_profile_reporter_evidence)
        self.assertFalse(args.tmt_reporter_evidence)


if __name__ == "__main__":
    unittest.main()
