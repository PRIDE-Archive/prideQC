"""Real binding and mzML reader gates, run on minimum and current OpenMS."""

import base64
import gzip
import hashlib
import importlib.util
import json
import os
import shutil
import struct
import tempfile
import unittest
import xml.etree.ElementTree as ET
from enum import Enum
from pathlib import Path

from prideqc.pipeline import Analyzer, Workflow, WorkflowOptions
from prideqc.readers import _enum_int, _enum_names

HAS_OPENMS = importlib.util.find_spec("pyopenms") is not None
FIXTURE_DIRECTORY = Path(__file__).parent / "data"
MZML_NAMESPACE = {"mz": "http://psi.hupo.org/ms/mzml"}


class MzMLFixtureTests(unittest.TestCase):
    """Check the input oracle even when pyOpenMS is unavailable locally."""

    def test_explicit_negative_faims_values_and_peak_arrays(self):
        for kind in ("indexed", "nonindexed"):
            with self.subTest(kind=kind):
                root = ET.parse(FIXTURE_DIRECTORY / f"faims-{kind}.mzML").getroot()
                spectra = root.findall(".//mz:spectrum", MZML_NAMESPACE)
                self.assertEqual(len(spectra), 4)
                for index, spectrum in enumerate(spectra):
                    values = spectrum.findall(
                        ".//mz:cvParam[@accession='MS:1001581']", MZML_NAMESPACE
                    )
                    self.assertEqual(
                        [float(node.attrib["value"]) for node in values],
                        [-45.0] if index % 2 == 0 else [],
                    )
                    self.assertTrue(all(node.get("unitAccession") == "UO:0000218" for node in values))
                    arrays = spectrum.findall(".//mz:binaryDataArray", MZML_NAMESPACE)
                    self.assertEqual(len(arrays), 2)
                    for array, expected in zip(arrays, ((100., 200.), (10., 20.)), strict=True):
                        encoded = array.findtext("mz:binary", namespaces=MZML_NAMESPACE)
                        self.assertEqual(len(encoded), int(array.attrib["encodedLength"]))
                        self.assertEqual(struct.unpack("<2d", base64.b64decode(encoded)), expected)

    def test_index_offsets_and_checksum_match_fixture_bytes(self):
        data = (FIXTURE_DIRECTORY / "faims-indexed.mzML").read_bytes()
        root = ET.fromstring(data)
        self.assertEqual(root.tag, "{http://psi.hupo.org/ms/mzml}indexedmzML")
        offsets = root.findall("mz:indexList/mz:index/mz:offset", MZML_NAMESPACE)
        self.assertEqual(len(offsets), 4)
        for offset in offsets:
            start = int(offset.text)
            opening_tag = data[start:data.index(b">", start) + 1]
            self.assertTrue(opening_tag.startswith(b"<spectrum "))
            self.assertIn(f'id="{offset.attrib["idRef"]}"'.encode(), opening_tag)
        index_start = int(root.findtext("mz:indexListOffset", namespaces=MZML_NAMESPACE))
        self.assertTrue(data[index_start:].startswith(b"<indexList "))
        checksum_end = data.index(b"<fileChecksum>") + len(b"<fileChecksum>")
        self.assertEqual(
            hashlib.sha1(data[:checksum_end], usedforsecurity=False).hexdigest(),
            root.findtext("mz:fileChecksum", namespaces=MZML_NAMESPACE),
        )

    def test_indexed_and_nonindexed_spectrum_content_is_identical(self):
        plain = ET.parse(FIXTURE_DIRECTORY / "faims-nonindexed.mzML").getroot()
        indexed = ET.parse(FIXTURE_DIRECTORY / "faims-indexed.mzML").getroot()
        self.assertEqual(plain.tag, "{http://psi.hupo.org/ms/mzml}mzML")
        nested = indexed.find("mz:mzML", MZML_NAMESPACE)
        self.assertIsNotNone(nested)
        nested.tail = None  # Whitespace outside mzML belongs to the index wrapper.
        self.assertEqual(ET.tostring(plain), ET.tostring(nested))


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
        # The 3.5.0 writer omits the negative FAIMS term in the former fixture.
        # Read explicit mzML bytes so the reader oracle is independent of store().
        path = self.root / name
        kind = "indexed" if indexed else "nonindexed"
        shutil.copyfile(FIXTURE_DIRECTORY / f"faims-{kind}.mzML", path)
        return path

    def test_indexed_and_nonindexed_mzml(self):
        for indexed in (True, False):
            with self.subTest(indexed=indexed):
                path = self.fixture(f"indexed-{indexed}.mzML", indexed=indexed)
                # Localize failures to native parsing before testing our adapter.
                experiment = self.oms.MSExperiment()
                self.oms.MzMLFile().load(str(path), experiment)
                native_values = [
                    float(spectrum.getDriftTime()) for spectrum in experiment
                    if _enum_int(spectrum.getDriftTimeUnit())
                    == _enum_int(self.oms.DriftTimeUnit.FAIMS_COMPENSATION_VOLTAGE)
                ]
                self.assertEqual(native_values, [-45., -45.], "Native mzML reader lost FAIMS CV")
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
        metrics = {m.key: m.value for m in result.metrics}
        self.assertEqual(metrics["NumberOfSpectra"], 4)
        self.assertEqual(metrics["FAIMS_CV_Values"], [-45.])

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
