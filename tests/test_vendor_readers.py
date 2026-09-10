"""Capability and archive tests for native pyOpenMS vendor readers."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
from zipfile import ZipFile

from prideqc.models import AnalysisResult, RunMetadata, SpectrumSink
from prideqc.mzqc import MzQCWriter
from prideqc.pipeline import WorkflowOptions, _analyze_file
from prideqc.readers import (
    PyOpenMSReader,
    VendorReaderUnavailable,
    _metadata_from_experiment,
    direct_vendor_support,
)


class _EmptyExperiment:
    def getSpectra(self):
        return []

    def getChromatograms(self):
        return []


class _Enums:
    POSITIVE = 1
    NEGATIVE = 2


class _FakeOMS:
    __version__ = "3.6.0.dev"
    MSExperiment = _EmptyExperiment

    class Precursor:
        ActivationMethod = _Enums

    class IonSource:
        Polarity = _Enums

    class SpectrumSettings:
        SpectrumType = _Enums

    class DriftTimeUnit:
        FAIMS_COMPENSATION_VOLTAGE = 1


class _Sink(SpectrumSink):
    def consume_spectrum(self, spectrum):
        pass

    def consume_chromatogram(self, rt, kind):
        pass


class VendorReaderTests(unittest.TestCase):
    def test_thermo_reader_is_feature_detected_and_used(self):
        calls = []

        class ThermoRawFile:
            def load(self, path, experiment):
                calls.append(Path(path))

        oms = type("OMS", (_FakeOMS,), {"ThermoRawFile": ThermoRawFile})
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "sample.raw"
            path.touch()
            reader = PyOpenMSReader(oms=oms)
            reader.read(path, _Sink())
        self.assertEqual(calls, [path])
        self.assertTrue(direct_vendor_support(oms, path))

    def test_bruker_directory_reader_is_feature_detected(self):
        calls = []

        class BrukerTimsFile:
            def load(self, path):
                calls.append(Path(path))
                return _EmptyExperiment()

        oms = type("OMS", (_FakeOMS,), {"BrukerTimsFile": BrukerTimsFile})
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "sample.d"
            path.mkdir()
            PyOpenMSReader(oms=oms).read(path, _Sink())
        self.assertEqual(calls, [path])

    def test_d_zip_falls_back_to_a_temporary_extracted_directory(self):
        calls = []

        class BrukerTimsFile:
            def load(self, path):
                calls.append(Path(path))
                return _EmptyExperiment()

        oms = type("OMS", (_FakeOMS,), {"BrukerTimsFile": BrukerTimsFile})
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            archive = root / "sample.d.zip"
            with ZipFile(archive, "w") as handle:
                handle.writestr("sample.d/analysis.tdf", "data")
            PyOpenMSReader(oms=oms).read(archive, _Sink())
        self.assertEqual(len(calls), 1)
        self.assertTrue(calls[0].name == "sample.d")
        self.assertFalse(calls[0].exists(), "temporary extraction must be cleaned up")

    def test_missing_reader_has_actionable_upgrade_message(self):
        oms = type("OMS", (_FakeOMS,), {"__version__": "3.5.0"})
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "sample.raw"
            path.touch()
            with self.assertRaises(VendorReaderUnavailable) as raised:
                PyOpenMSReader(oms=oms).read(path, _Sink())
            message = str(raised.exception)
            self.assertIn("ThermoRawFile", message)
            self.assertIn("development", message)
            self.assertIn("converter", message)
            self.assertIn("https://pypi.openms.de/simple/pyopenms/", message)

    def test_version_alone_never_claims_thermo_support(self):
        oms = type("OMS", (_FakeOMS,), {"__version__": "3.6.0"})
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "sample.raw"
            path.touch()
            self.assertFalse(PyOpenMSReader(oms=oms).supports_direct(path))

    def test_bruker_reader_unavailable_is_reported(self):
        oms = type("OMS", (_FakeOMS,), {"__version__": "3.5.0"})
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "sample.d"
            path.mkdir()
            with self.assertRaisesRegex(VendorReaderUnavailable, "BrukerTimsFile"):
                PyOpenMSReader(oms=oms).read(path, _Sink())

    def test_generic_filehandler_does_not_claim_vendor_support(self):
        class FileHandler:
            def load(self, path, experiment):
                raise AssertionError("must not be called")

        oms = type("OMS", (_FakeOMS,), {"FileHandler": FileHandler})
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "sample.raw"
            path.touch()
            self.assertFalse(PyOpenMSReader(oms=oms).supports_direct(path))

    def test_vendor_experiment_metadata_is_preserved_without_cv_guessing(self):
        class SourceFile:
            def getNameOfFile(self):
                return "original.raw"

        class Instrument:
            def getName(self):
                return "Orbitrap"

            def getVendor(self):
                return "Thermo Fisher"

            def getModel(self):
                return "Exploris"

            def getSerialNumber(self):
                return "SN-1"

            def getMassAnalyzers(self):
                return []

            def getIonSources(self):
                return []

        class Settings:
            def getSourceFiles(self):
                return [SourceFile()]

            def getInstrument(self):
                return Instrument()

        class Experiment:
            def getExperimentalSettings(self):
                return Settings()

        metadata = _metadata_from_experiment(Experiment(), Path("sample.raw"))
        self.assertEqual(metadata.source_files, ["sample.raw", "original.raw"])
        self.assertEqual(metadata.serial_numbers, ["SN-1"])
        self.assertEqual(metadata.instrument_details["model"], "Exploris")
        self.assertEqual(metadata.instruments, [])

    def test_mzqc_file_format_follows_original_input(self):
        for filename, accession, name in (
            ("sample.mzML", "MS:1000584", "mzML format"),
            ("sample.raw", "MS:1000563", "Thermo RAW format"),
            ("sample.d", "MS:1002817", "Bruker TDF format"),
            ("sample.d.zip", "MS:1002817", "Bruker TDF format"),
        ):
            with self.subTest(filename=filename):
                result = AnalysisResult(
                    Path(filename), RunMetadata(), [], [], [], "test", 0.0,
                )
                input_file = MzQCWriter().build(result, Path("terms.obo"))["mzQC"][
                    "runQualities"
                ][0]["metadata"]["inputFiles"][0]
                self.assertEqual(input_file["fileFormat"], {"accession": accession, "name": name})

    def test_workflow_uses_converter_only_when_native_reader_is_unavailable(self):
        class Reader:
            engine_version = "3.5.0"

            def supports_direct(self, path):
                return False

            def unavailable_error(self, path):
                return VendorReaderUnavailable("upgrade pyOpenMS or use a converter")

        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            source = root / "sample.raw"
            source.touch()
            converted = root / "sample.mzML"
            converted.touch()
            with patch("prideqc.readers.PyOpenMSReader", return_value=Reader()), patch(
                "prideqc.pipeline.ExternalConverter",
            ) as converter_type, patch("prideqc.pipeline.Analyzer") as analyzer_type:
                converter_type.return_value.convert.return_value = converted
                analyzer_type.return_value.analyze.return_value = SimpleNamespace(source_path=None)
                outcome = _analyze_file((source, root / "out", WorkflowOptions(converter="msconvert")))
        self.assertIsNone(outcome.error)
        converter_type.return_value.convert.assert_called_once()
        analyzer_type.return_value.analyze.assert_called_once()


if __name__ == "__main__":
    unittest.main()
