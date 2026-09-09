"""Real binding and mzML round-trip gates, run on minimum and current OpenMS."""

import gzip
import importlib.util
import json
import os
import shutil
import tempfile
import unittest
from enum import Enum
from pathlib import Path

import numpy as np

from prideqc.pipeline import Analyzer, Workflow, WorkflowOptions
from prideqc.readers import _enum_int, _enum_names

HAS_OPENMS = importlib.util.find_spec("pyopenms") is not None


class EnumCompatibilityTests(unittest.TestCase):
    def test_integer_and_python_enum_bindings(self):
        class Activation(Enum):
            HCD = 1
            EThcD = 2

        self.assertEqual(_enum_int(1), 1)
        self.assertEqual(_enum_int(Activation.EThcD), 2)
        self.assertEqual(_enum_names(Activation)[2], "EThcD")


@unittest.skipUnless(HAS_OPENMS, "pyOpenMS not installed")
class PyOpenMSIntegrationTests(unittest.TestCase):
    def setUp(self):
        import pyopenms as oms

        self.oms = oms
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)

    def fixture(self, name="run.mzML", indexed=True):
        oms = self.oms
        experiment = oms.MSExperiment()
        instrument = experiment.getInstrument()
        instrument.setName("Orbitrap Fusion")
        experiment.setInstrument(instrument)
        for i in range(4):
            spectrum = oms.MSSpectrum()
            spectrum.setNativeID(f"scan={i + 1}")
            spectrum.setMSLevel(1 if i % 2 == 0 else 2)
            spectrum.setRT(float(i * 10))
            spectrum.setType(oms.SpectrumSettings.SpectrumType.CENTROID)
            spectrum.set_peaks((np.array([100., 200.]), np.array([10., 20.])))
            settings = spectrum.getInstrumentSettings()
            settings.setPolarity(oms.IonSource.Polarity.POSITIVE)
            spectrum.setInstrumentSettings(settings)
            if i % 2:
                precursor = oms.Precursor()
                precursor.setMZ(500.)
                precursor.setCharge(2)
                precursor.setIntensity(123.)
                precursor.setIsolationWindowLowerOffset(1.)
                precursor.setIsolationWindowUpperOffset(1.)
                precursor.setActivationMethods({oms.Precursor.ActivationMethod.HCD})
                precursor.setActivationEnergy(30.)
                spectrum.setPrecursors([precursor])
            else:
                spectrum.setDriftTime(-45.)
                spectrum.setDriftTimeUnit(oms.DriftTimeUnit.FAIMS_COMPENSATION_VOLTAGE)
            experiment.addSpectrum(spectrum)
        path = self.root / name
        writer = oms.MzMLFile()
        options = writer.getOptions()
        options.setWriteIndex(indexed)
        writer.setOptions(options)
        writer.store(str(path), experiment)
        return path

    def test_indexed_and_nonindexed_mzml(self):
        for indexed in (True, False):
            with self.subTest(indexed=indexed):
                path = self.fixture(f"indexed-{indexed}.mzML", indexed=indexed)
                result = Analyzer().analyze(path)
                metrics = {m.key: m.value for m in result.metrics}
                self.assertEqual(metrics["NumberOfSpectra_MS1"], 2)
                self.assertEqual(metrics["NumberOfSpectra_MS2"], 2)
                self.assertEqual(metrics["TIC_MS1_Area"], 600.)
                self.assertEqual(metrics["FAIMS_CV_Values"], [-45.])
                self.assertEqual(metrics["Polarity_MS1"]["polarity"], ["positive"])
                evidence = next(a for a in result.annotations if a.field == "dissociation_ms2")
                self.assertIn("MS:1000422", evidence.sdrf_value)

    def test_gzip_mzml(self):
        path = self.fixture()
        compressed = path.with_suffix(".mzML.gz")
        with path.open("rb") as source, gzip.open(compressed, "wb") as target:
            shutil.copyfileobj(source, target)
        result = Analyzer().analyze(compressed)
        self.assertEqual(result.input_path, compressed)
        self.assertEqual({m.key: m.value for m in result.metrics}["NumberOfSpectra"], 4)

    def test_serial_and_parallel_results_match(self):
        paths = [self.fixture("a.mzML"), self.fixture("b.mzML")]
        serial = Workflow().run(paths, self.root / "serial")
        parallel = Workflow(WorkflowOptions(workers=2)).run(paths, self.root / "parallel")
        self.assertTrue(serial["success"])
        self.assertTrue(parallel["success"])
        for path in paths:
            def metrics(folder, path=path):
                content = json.loads((self.root / folder / f"{path.name}.summary.json").read_text())
                return content["metrics"]
            self.assertEqual(metrics("serial"), metrics("parallel"))

    def test_corrupt_mzml_is_failed_not_zero_quality(self):
        path = self.root / "bad.mzML"
        path.write_text("not XML")
        manifest = Workflow().run([path], self.root / "out")
        self.assertFalse(manifest["success"])
        self.assertFalse((self.root / "out" / "bad.mzML.mzQC").exists())


class RequiredDependencyGate(unittest.TestCase):
    def test_ci_dependencies_cannot_silently_skip(self):
        if os.environ.get("PRIDE_QC_REQUIRE_INTEGRATION") == "1":
            self.assertTrue(HAS_OPENMS, "pyOpenMS missing from CI environment")
            self.assertIsNotNone(
                importlib.util.find_spec("jsonschema"),
                "jsonschema missing from CI environment",
            )
