"""Capability and archive tests for native pyOpenMS vendor readers."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
from zipfile import ZipFile

from prideqc.models import SpectrumSink
from prideqc.pipeline import WorkflowOptions, _analyze_file
from prideqc.readers import PyOpenMSReader, VendorReaderUnavailable, direct_vendor_support


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
            def load(self, path, experiment):
                calls.append(Path(path))

        oms = type("OMS", (_FakeOMS,), {"BrukerTimsFile": BrukerTimsFile})
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "sample.d"
            path.mkdir()
            PyOpenMSReader(oms=oms).read(path, _Sink())
        self.assertEqual(calls, [path])

    def test_d_zip_falls_back_to_a_temporary_extracted_directory(self):
        calls = []

        class BrukerTimsFile:
            def load(self, path, experiment):
                calls.append(Path(path))
                if str(path).endswith(".zip"):
                    raise ValueError("archive not supported")

        oms = type("OMS", (_FakeOMS,), {"BrukerTimsFile": BrukerTimsFile})
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            archive = root / "sample.d.zip"
            with ZipFile(archive, "w") as handle:
                handle.writestr("sample.d/analysis.tdf", "data")
            PyOpenMSReader(oms=oms).read(archive, _Sink())
        self.assertEqual(calls[0], archive)
        self.assertTrue(calls[1].name == "sample.d")
        self.assertFalse(calls[1].exists(), "temporary extraction must be cleaned up")

    def test_missing_reader_has_actionable_upgrade_message(self):
        oms = type("OMS", (_FakeOMS,), {"__version__": "3.5.0"})
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "sample.raw"
            path.touch()
            with self.assertRaisesRegex(VendorReaderUnavailable, r">=3\.6.*converter"):
                PyOpenMSReader(oms=oms).read(path, _Sink())

    def test_filehandler_capability_is_supported_on_new_builds(self):
        calls = []

        class FileHandler:
            def loadExperiment(self, path, experiment):
                calls.append(Path(path))

        oms = type("OMS", (_FakeOMS,), {"FileHandler": FileHandler})
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "sample.raw"
            path.touch()
            PyOpenMSReader(oms=oms).read(path, _Sink())
        self.assertEqual(calls, [path])

    def test_generic_old_filehandler_does_not_claim_vendor_support(self):
        class FileHandler:
            def loadExperiment(self, path, experiment):
                raise AssertionError("must not be called")

        oms = type("OMS", (_FakeOMS,), {"__version__": "3.5.0", "FileHandler": FileHandler})
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "sample.raw"
            path.touch()
            self.assertFalse(PyOpenMSReader(oms=oms).supports_direct(path))

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
