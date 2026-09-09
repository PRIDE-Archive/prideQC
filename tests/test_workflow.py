"""Evidence, SDRF, serialization and orchestration contracts."""

import importlib.util
import json
import math
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np

from prideqc.annotations import DiagnosticIonCollector
from prideqc.cli import main
from prideqc.conversion import ExternalConverter
from prideqc.io import atomic_text, json_safe
from prideqc.models import Annotation, EvidenceKind, Metric, RunMetadata
from prideqc.mzqc import MzQCWriter, definition
from prideqc.pipeline import Analyzer, FileOutcome, Workflow, WorkflowOptions
from prideqc.readers import read_header
from prideqc.sdrf import SDRFDocument
from tests.helpers import FUSION, MemoryReader, spectrum


def analyze(spectra=None, name="run.mzML", metadata=None):
    reader = MemoryReader(spectra or [spectrum(0, [10]), spectrum(10, [20])], metadata)
    return Analyzer(reader).analyze(name)


class AnnotationTests(unittest.TestCase):
    def test_one_pass_shared_by_metrics_and_diagnostics(self):
        collector = DiagnosticIonCollector()
        reader = MemoryReader(
            [spectrum(0, [10, 10, 10], 2, mz=[126.127726, 127.131081, 128.134436])],
        )
        result = Analyzer(reader).analyze("a.mzML", collectors=[collector])
        self.assertEqual(reader.read_count, 1)
        evidence = next(a for a in result.annotations if a.field == "signature_TMT_family")
        self.assertEqual(evidence.support, 1)
        self.assertIsNone(evidence.sdrf_value)

    def test_profile_and_unknown_spectra_do_not_become_reporter_evidence(self):
        collector = DiagnosticIonCollector()
        collector.consume_spectrum(
            spectrum(
                0,
                [10, 10, 10],
                2,
                mz=[126.127726, 127.131081, 128.134436],
                representation="profile",
            ),
        )
        self.assertEqual(collector.total, 0)
        self.assertEqual(collector.skipped_profile, 1)

    def test_acquisition_heuristic_is_not_observed_fact(self):
        result = analyze([spectrum(i, [1], 2, charge=2, width=20) for i in range(100)])
        evidence = next(a for a in result.annotations if a.field == "acquisition_method")
        self.assertEqual(evidence.kind, EvidenceKind.INFERRED)
        self.assertIn("PRIDE:0000450", evidence.sdrf_value)

    def test_sparse_or_mixed_widths_abstain(self):
        for widths in ([20] * 10, [2] * 50 + [20] * 50):
            result = analyze([spectrum(i, [1], 2, charge=2, width=w) for i, w in enumerate(widths)])
            evidence = next(a for a in result.annotations if a.field == "acquisition_method")
            self.assertIsNone(evidence.value)

    def test_tolerance_is_not_guessed_from_instrument_or_width(self):
        result = analyze()
        evidence = next(a for a in result.annotations if a.field == "precursor_mass_tolerance")
        self.assertEqual(evidence.kind, EvidenceKind.UNAVAILABLE)
        self.assertIsNone(evidence.sdrf_value)


class SDRFTests(unittest.TestCase):
    def test_not_applicable_is_preserved_as_an_existing_assertion(self):
        document = SDRFDocument(
            ["comment[data file]", "comment[instrument]"],
            [["run.mzML", "not applicable"]],
        )
        self.assertEqual(document.annotate([analyze()])[0].status, "conflict")
        self.assertEqual(document.rows[0][1], "not applicable")

    def test_duplicate_modification_columns_and_repeated_sample_rows_survive(self):
        columns = ["source name", "comment[data file]", "comment[modification parameters]", "comment[modification parameters]"]
        document = SDRFDocument(columns.copy(), [["s1", "run.raw", "Oxidation", "Carbamidomethyl"],
                                                 ["s2", "run.raw", "", ""]], ["# meta\n"])
        result = analyze(metadata=RunMetadata(instruments=[FUSION], source_files=["run.raw"]))
        changes = document.annotate([result])
        self.assertEqual(document.columns[:4], columns)
        self.assertEqual(document.rows[0][:4], ["s1", "run.raw", "Oxidation", "Carbamidomethyl"])
        self.assertEqual(sum(c.status == "filled" for c in changes), 2)
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "out.tsv"
            document.write(path)
            again = SDRFDocument.read(path)
            self.assertEqual(again.columns, document.columns)
            self.assertEqual(again.rows, document.rows)
            self.assertEqual(again.preamble, ["# meta\n"])

    def test_conflict_preserved_until_explicit_overwrite(self):
        document = SDRFDocument(
            ["comment[data file]", "comment[instrument]"],
            [["run.mzML", "existing model"]],
        )
        changes = document.annotate([analyze()])
        self.assertEqual(changes[0].status, "conflict")
        self.assertEqual(document.rows[0][1], "existing model")
        document.annotate([analyze()], overwrite=True)
        self.assertIn("MS:1002416", document.rows[0][1])

    def test_cv_accession_match_preserves_original_label(self):
        document = SDRFDocument(
            ["comment[data file]", "comment[instrument]"],
            [["run.mzML", "AC=MS:1002416;NT=old spelling"]],
        )
        self.assertEqual(document.annotate([analyze()])[0].status, "match")
        self.assertEqual(document.rows[0][1], "AC=MS:1002416;NT=old spelling")

    def test_no_implicit_stem_matching(self):
        document = SDRFDocument(["comment[data file]"], [["run.raw"]])
        self.assertEqual(document.annotate([analyze()])[0].status, "unmatched")
        changes = document.annotate([analyze()], aliases={"run.raw": "run.mzML"})
        self.assertEqual(changes[0].status, "filled")

    def test_ambiguous_source_alias_is_rejected_before_mutation(self):
        metadata = RunMetadata(instruments=[FUSION], source_files=["run.raw"])
        document = SDRFDocument(["comment[data file]"], [["run.raw"]])
        with self.assertRaisesRegex(ValueError, "Ambiguous"):
            document.annotate(
                [analyze(
                    name="a.mzML",
                    metadata=metadata,
                ), analyze(
                    name="b.mzML",
                    metadata=metadata,
                )],
            )
        self.assertEqual(len(document.columns), 1)

    def test_inferred_annotation_requires_explicit_option(self):
        result = analyze()
        result.annotations = [Annotation(
            "mode",
            "DIA",
            EvidenceKind.INFERRED,
            "test",
            sdrf_column="comment[proteomics data acquisition method]",
            sdrf_value="NT=Data-independent acquisition;AC=PRIDE:0000450",
        )]
        document = SDRFDocument(["comment[data file]"], [["run.mzML"]])
        self.assertEqual(document.annotate([result])[0].status, "suggestion")
        self.assertEqual(len(document.columns), 1)
        self.assertEqual(document.annotate([result], include_inferred=True)[0].status, "filled")

    def test_duplicate_target_columns_are_reported_without_change(self):
        document = SDRFDocument(
            ["comment[data file]", "comment[instrument]", "comment[instrument]"],
            [["run.mzML", "", ""]],
        )
        self.assertEqual(document.annotate([analyze()])[0].status, "ambiguous_column")
        self.assertEqual(document.rows[0][1:], ["", ""])

    def test_invalid_tsv_width_or_duplicate_data_file_column_rejected(self):
        with self.assertRaises(ValueError):
            SDRFDocument(["comment[data file]"], [["a", "extra"]])
        with self.assertRaises(ValueError):
            SDRFDocument(["comment[data file]", "comment[data file]"], [])


class SerializationTests(unittest.TestCase):
    def test_json_nonfinite_values_are_null_recursively(self):
        value = json_safe({"a": np.array([1, math.nan]), "b": {"x": math.inf}})
        self.assertEqual(value, {"a": [1.0, None], "b": {"x": None}})
        json.dumps(value, allow_nan=False)

    def test_mzqc_one_run_has_complete_provenance_and_no_duplicate_accessions(self):
        result = analyze()
        data = MzQCWriter().build(result, Path("run.obo"))["mzQC"]
        self.assertEqual(len(data["runQualities"]), 1)
        run = data["runQualities"][0]
        self.assertTrue(run["metadata"]["inputFiles"][0]["location"].startswith("file://"))
        self.assertEqual(run["metadata"]["analysisSoftware"][0]["name"], "prideqc")
        accessions = [m["accession"] for m in run["qualityMetrics"]]
        self.assertEqual(len(accessions), len(set(accessions)))
        self.assertNotIn("setQualities", data)
        self.assertTrue(data["creationDate"].endswith("+00:00"))

    @unittest.skipUnless(importlib.util.find_spec("jsonschema"), "jsonschema not installed")
    def test_official_mzqc_schema_with_format_checks(self):
        import jsonschema

        schema = json.loads((Path(__file__).parent / "data" / "mzqc_schema.json").read_text())
        for items in ([], [spectrum(0, [1])], [spectrum(0, [1], 2, charge=2)]):
            result = Analyzer(MemoryReader(items)).analyze("run.mzML")
            jsonschema.Draft7Validator(schema, format_checker=jsonschema.FormatChecker()).validate(
                MzQCWriter().build(result, Path("local.obo")))

    def test_cv_shape_stable_when_unavailable(self):
        self.assertEqual(definition(Metric("MzRange_MS2", None)).shape, "tuple")
        self.assertEqual(definition(Metric("MzRange_MS2", [400, 600])).shape, "tuple")

    def test_unknown_metric_cannot_be_mislabeled(self):
        with self.assertRaisesRegex(ValueError, "Unregistered"):
            definition(Metric("made_up", 1))

    def test_atomic_write_failure_leaves_previous_file(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "result.json"
            path.write_text("old")
            with self.assertRaises(RuntimeError), atomic_text(path) as handle:
                handle.write("partial")
                raise RuntimeError("failure")
            self.assertEqual(path.read_text(), "old")
            self.assertEqual(len(list(Path(folder).iterdir())), 1)


class ReaderHeaderTests(unittest.TestCase):
    def test_preserves_cv_and_reference_groups_without_reading_payload(self):
        xml = '''<mzML xmlns="http://psi.hupo.org/ms/mzml">
        <referenceableParamGroupList><referenceableParamGroup id="inst">
        <cvParam accession="MS:1002416" name="Orbitrap Fusion" value=""/>
        </referenceableParamGroup></referenceableParamGroupList>
        <fileDescription><sourceFileList><sourceFile id="x" name="run.raw" location="file:///data"/></sourceFileList></fileDescription>
        <instrumentConfigurationList><instrumentConfiguration id="IC1">
        <referenceableParamGroupRef ref="inst"/>
        <cvParam accession="MS:1000529" name="instrument serial number" value="ABC"/>
        <componentList><analyzer><cvParam accession="MS:1000484" name="orbitrap"/></analyzer></componentList>
        </instrumentConfiguration></instrumentConfigurationList>
        <run startTimeStamp="2026-01-01T00:00:00Z"><spectrumList count="1">'''
        # Malformed content is far beyond the header parser's read-ahead buffer.
        xml += " " * 65536 + "<invalid payload"
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "header.mzML"
            path.write_text(xml)
            metadata = read_header(path)
        self.assertEqual(metadata.instruments, [FUSION])
        self.assertEqual(metadata.serial_numbers, ["ABC"])
        self.assertEqual(metadata.source_files, ["run.raw"])
        self.assertEqual(metadata.analyzers[0].accession, "MS:1000484")


class WorkflowTests(unittest.TestCase):
    def test_partial_failure_writes_manifest_and_successful_mzqc(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            good, bad = root / "good.mzML", root / "bad.mzML"
            good.touch()
            bad.touch()
            result = analyze(name=str(good))
            with patch(
                "prideqc.pipeline._analyze_file",
                side_effect=[FileOutcome(good, result), FileOutcome(bad, error="broken file")],
            ):
                manifest = Workflow(WorkflowOptions(continue_on_error=True)).run(
                    [good, bad],
                    root / "out",
                )
            self.assertFalse(manifest["success"])
            self.assertTrue((root / "out/good.mzML.mzQC").exists())
            self.assertFalse((root / "out/bad.mzML.mzQC").exists())
            self.assertEqual(
                json.loads((root / "out/manifest.json").read_text())["files"][1]["status"],
                "failed",
            )

    def test_output_collision_and_existing_results_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            source = root / "run.mzML"
            source.touch()
            with self.assertRaises(ValueError):
                Workflow().run([source, source], root / "out")
            (root / "out").mkdir()
            (root / "out" / "existing").touch()
            with self.assertRaises(FileExistsError):
                Workflow().run([source], root / "out")

    def test_fail_fast_records_unprocessed_files(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            first, second = root / "a.mzML", root / "b.mzML"
            first.touch()
            second.touch()
            with patch(
                "prideqc.pipeline._analyze_file",
                return_value=FileOutcome(first, error="bad"),
            ) as worker:
                manifest = Workflow().run([first, second], root / "out")
            self.assertEqual(worker.call_count, 1)
            self.assertEqual(manifest["unprocessed"], [str(second)])

    def test_cli_failed_run_returns_nonzero(self):
        with patch(
            "prideqc.pipeline.Workflow.run",
            return_value={"files": [], "errors": [], "success": False},
        ):
            self.assertEqual(main(["analyze", "a.mzML", "-o", "out"]), 1)

    def test_help_does_not_import_runtime_dependencies(self):
        code = (
            "from prideqc.cli import parser; import sys; parser(); "
            "assert not {'pyopenms', 'numpy', 'pridepy', 'sdrf_pipelines'} & sys.modules.keys()"
        )
        subprocess.run([sys.executable, "-c", code], check=True)

    def test_converter_arguments_preserve_paths_and_never_use_shell(self):
        converter = ExternalConverter("msconvert")
        source = Path("/data/run; echo secret.wiff")
        command = converter.command(source, Path("/data/output folder"))
        self.assertEqual(command[1], str(source))
        self.assertEqual(command[-1], "/data/output folder")

    def test_input_validation(self):
        with self.assertRaises(ValueError):
            WorkflowOptions(workers=0)
        with self.assertRaises(ValueError):
            DiagnosticIonCollector(ppm=math.nan)
