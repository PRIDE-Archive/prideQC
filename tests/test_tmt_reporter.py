"""Offline regression contract for one-pass TMT channel evidence."""

from __future__ import annotations

import math
import unittest
from types import SimpleNamespace

from prideqc.models import EvidenceKind
from prideqc.mzqc import annotation_definition
from prideqc.pipeline import Analyzer
from prideqc.tmt import TMTReporterCollector
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


class _MockTMT16:
    def getChannelInformation(self):
        return [
            SimpleNamespace(name=name, center=mz)
            for name, mz in (
                ("126", 126.127725),
                ("127N", 127.124760),
                ("127C", 127.131081),
                ("134N", 134.148245),
            )
        ]


class _MockOMS:
    __version__ = "test-oms"
    TMTSixPlexQuantitationMethod = _MockTMT6
    TMTSixteenPlexQuantitationMethod = _MockTMT16


def _report(collector: TMTReporterCollector):
    result = collector.annotations()[0]
    return result, result.value["levels"]


class TMTReporterTests(unittest.TestCase):
    def test_ms2_ms3_are_independent_and_catalogs_are_not_plex_calls(self):
        collector = TMTReporterCollector(oms=_MockOMS(), ppm=10)
        peaks = [126.127725, 127.124760, 127.131081, 128.134433]
        collector.consume_spectrum(spectrum(0, [1000, 900, 800, 700], 3, mz=peaks))
        collector.consume_spectrum(spectrum(1, [1000, 900, 800, 700], 2, mz=peaks))
        collector.consume_spectrum(
            spectrum(2, [1000, 900, 800, 700], 2, mz=peaks, representation="profile")
        )
        result, levels = _report(collector)
        self.assertEqual(result.kind, EvidenceKind.INFERRED)
        self.assertIsNone(result.sdrf_column)
        self.assertIsNone(result.sdrf_value)
        self.assertIsNone(result.support)  # Catalogs overlap, no sum of independent spectra.
        self.assertEqual(result.total, 2)
        self.assertEqual(levels["MS2"]["centroid_spectra"], 1)
        self.assertEqual(levels["MS2"]["profile_skipped"], 1)
        self.assertEqual(levels["MS3"]["methods"]["TMT6"]["reporter_positive_spectra"], 1)
        self.assertEqual(levels["MS2"]["methods"]["TMTpro16"]["reporter_positive_spectra"], 1)
        self.assertIn("TMT11", result.value["unavailable_methods"])
        self.assertEqual(result.value["pyopenms_version"], "test-oms")
        self.assertEqual(
            result.value["scientific_status"],
            "EXPLORATORY_NOT_FOR_AUTOMATIC_SDRF_WRITEBACK",
        )

    def test_peak_selection_uses_highest_in_window_and_records_ppm_error(self):
        collector = TMTReporterCollector(oms=_MockOMS(), min_relative_intensity=0.001)
        collector.consume_spectrum(
            spectrum(
                0,
                [100, 500, 300, 200, 1],
                3,
                mz=[126.127725, 126.127725 * (1 + 2 / 1e6), 127.124760, 128.134433, 140.0],
            )
        )
        _, levels = _report(collector)
        channels = levels["MS3"]["methods"]["TMT6"]["channels"]
        first = channels[0]
        self.assertEqual(first["spectra_with_hit"], 1)
        self.assertEqual(first["uncorrected_total_matched_intensity"], 500)
        self.assertAlmostEqual(first["mean_absolute_mass_error_ppm"], 2.0)
        self.assertEqual(levels["MS3"]["methods"]["TMT6"]["reporter_positive_spectra"], 1)

    def test_weak_extra_channel_does_not_imply_an_unused_channel(self):
        collector = TMTReporterCollector(oms=_MockOMS(), min_relative_intensity=0.001)
        for _ in range(10):
            collector.consume_spectrum(
                spectrum(
                    0, [1000, 900, 800, 1], 2, mz=[126.127725, 127.124760, 127.131081, 134.148245]
                )
            )
        _, levels = _report(collector)
        channels = levels["MS2"]["methods"]["TMTpro16"]["channels"]
        weak = next(ch for ch in channels if ch["channel"] == "134N")
        self.assertEqual(weak["spectra_with_hit"], 10)
        self.assertEqual(weak["uncorrected_total_matched_intensity"], 10)
        self.assertIsNone(collector.annotations()[0].sdrf_value)

    def test_profile_unknown_and_empty_centroid_are_not_tmt_negative(self):
        collector = TMTReporterCollector(oms=_MockOMS())
        collector.consume_spectrum(
            spectrum(
                0,
                [100, 100, 100],
                2,
                mz=[126.127725, 127.124760, 128.134433],
                representation="profile",
            )
        )
        collector.consume_spectrum(spectrum(1, [100], 3, mz=[126.127725], representation="unknown"))
        result, levels = _report(collector)
        self.assertEqual(result.kind, EvidenceKind.UNAVAILABLE)
        self.assertIsNone(levels["MS2"]["methods"]["TMT6"]["positive_fraction_of_centroid"])
        self.assertEqual(levels["MS2"]["profile_skipped"], 1)
        self.assertEqual(levels["MS3"]["unknown_skipped"], 1)

    def test_bad_parameters_fail_closed(self):
        for params in (
            {"ppm": 0},
            {"ppm": math.nan},
            {"min_relative_intensity": math.nan},
            {"min_relative_intensity": 2},
            {"minimum_channel_hits": 0},
        ):
            with self.subTest(params=params), self.assertRaises(ValueError):
                TMTReporterCollector(oms=_MockOMS(), **params)
        with self.assertRaises(RuntimeError):
            TMTReporterCollector(oms=SimpleNamespace(__version__="empty"))

    def test_analyzer_is_one_pass_and_mzqc_definition_exists(self):
        reader = MemoryReader(
            [spectrum(0, [100, 90, 80], 3, mz=[126.127725, 127.124760, 128.134433])]
        )
        result = Analyzer(reader).analyze(
            "test.mzML", collectors=[TMTReporterCollector(oms=_MockOMS())]
        )
        self.assertEqual(reader.read_count, 1)
        evidence = next(a for a in result.annotations if a.field == "tmt_reporter_evidence")
        self.assertIn("MS3", evidence.value["levels"])
        term = annotation_definition(evidence)
        self.assertIsNotNone(term)
        self.assertEqual(term.shape, "table")
        self.assertTrue(term.accession.startswith("QCPRIDE:"))

    def test_pipeline_default_remains_opt_in(self):
        from prideqc.cli import parser
        from prideqc.pipeline import WorkflowOptions

        self.assertFalse(WorkflowOptions().tmt_reporter_evidence)
        self.assertTrue(WorkflowOptions(tmt_reporter_evidence=True).tmt_reporter_evidence)
        args = parser().parse_args(
            ["analyze", "run.mzML", "-o", "output", "--tmt-reporter-evidence"]
        )
        self.assertTrue(args.tmt_reporter_evidence)


if __name__ == "__main__":
    unittest.main()
