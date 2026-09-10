"""Maintained dependency delegation, acquisition failures and workflow contracts."""

import importlib.util
import io
import json
import os
import tempfile
import types
import unittest
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path
from unittest.mock import MagicMock, create_autospec, patch

from prideqc.cli import main
from prideqc.pipeline import FileOutcome, Workflow
from prideqc.pride import DownloadOptions, PrideRepository, selected_files
from prideqc.validation import SDRFPipelinesValidator, ValidationReport
from tests.test_workflow import analyze


def report(path, *issues):
    return ValidationReport(str(path), "ms-proteomics", False, tuple(issues), "test")


def write_download(**arguments):
    (Path(arguments["output_folder"]) / arguments["file_name"]).write_bytes(
        b"test mzML bytes"
    )


class ValidatorTests(unittest.TestCase):
    def native_modules(self, reader):
        module = types.ModuleType("sdrf_pipelines.sdrf.sdrf")
        module.read_sdrf = reader
        return {
            "sdrf_pipelines": types.ModuleType("sdrf_pipelines"),
            "sdrf_pipelines.sdrf": types.ModuleType("sdrf_pipelines.sdrf"),
            "sdrf_pipelines.sdrf.sdrf": module,
        }

    def test_native_parser_template_and_ontology_are_used(self):
        reader = MagicMock()
        reader.return_value.validate_sdrf.return_value = [
            types.SimpleNamespace(message="missing instrument"),
            "bad value",
        ]
        with patch.dict("sys.modules", self.native_modules(reader)):
            result = SDRFPipelinesValidator("dia-acquisition", ontology=True).validate(
                Path("study.tsv")
            )
        reader.assert_called_once_with("study.tsv")
        reader.return_value.validate_sdrf.assert_called_once_with(
            template="dia-acquisition", skip_ontology=False
        )
        self.assertEqual(result.issues, ("missing instrument", "bad value"))
        self.assertFalse(result.valid)
        self.assertEqual(result.to_dict()["engine"], "sdrf-pipelines")

    def test_default_validation_skips_ontology_and_does_not_mutate_input(self):
        reader = MagicMock()
        reader.return_value.validate_sdrf.return_value = []
        with patch.dict("sys.modules", self.native_modules(reader)):
            result = SDRFPipelinesValidator().validate(Path("study.tsv"))
        self.assertTrue(result.valid)
        reader.return_value.validate_sdrf.assert_called_once_with(
            template="ms-proteomics", skip_ontology=True
        )

    def test_missing_dependency_and_validator_failure_never_pass(self):
        with patch.dict("sys.modules", {"sdrf_pipelines.sdrf.sdrf": None}):
            with self.assertRaisesRegex(RuntimeError, "uv sync"):
                SDRFPipelinesValidator().validate(Path("study.tsv"))
        reader = MagicMock(side_effect=ValueError("malformed"))
        with patch.dict("sys.modules", self.native_modules(reader)):
            with self.assertRaisesRegex(RuntimeError, "could not validate"):
                SDRFPipelinesValidator().validate(Path("study.tsv"))

    def test_cli_validation_findings_exit_nonzero_and_write_report(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "report.json"
            with patch(
                "prideqc.validation.SDRFPipelinesValidator.validate",
                return_value=report("in.tsv", "missing column"),
            ):
                with redirect_stderr(io.StringIO()):
                    code = main(["validate-sdrf", "in.tsv", "--json", str(path)])
            self.assertEqual(code, 1)
            self.assertFalse(json.loads(path.read_text())["valid"])

    def test_workflow_records_input_and_output_and_fails_invalid_refinement(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            source = root / "run.mzML"
            source.touch()
            sdrf = root / "study.tsv"
            sdrf.write_text("comment[data file]\nrun.mzML\n")
            validator = MagicMock()
            validator.validate.side_effect = [
                report(sdrf, "missing instrument"),
                report("refined.tsv", "missing source"),
            ]
            with patch(
                "prideqc.pipeline._analyze_file",
                return_value=FileOutcome(source, analyze(name=str(source))),
            ):
                result = Workflow(validator=validator).run(
                    [source], root / "qc", sdrf=sdrf
                )
            self.assertFalse(result["success"])
            self.assertEqual(validator.validate.call_args_list[0].args, (sdrf,))
            self.assertEqual(
                validator.validate.call_args_list[1].args,
                (root / "qc/refined.sdrf.tsv",),
            )
            self.assertTrue((root / "qc/run.mzML.mzQC").exists())
            data = json.loads((root / "qc/sdrf-validation.json").read_text())
            self.assertFalse(data["refined"]["valid"])

    def test_repairable_input_does_not_fail_valid_refinement(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            source = root / "run.mzML"
            source.touch()
            sdrf = root / "study.tsv"
            sdrf.write_text("comment[data file]\nrun.mzML\n")
            validator = MagicMock()
            validator.validate.side_effect = [report(sdrf, "missing instrument"), report("refined.tsv")]
            with patch(
                "prideqc.pipeline._analyze_file",
                return_value=FileOutcome(source, analyze(name=str(source))),
            ):
                result = Workflow(validator=validator).run(
                    [source], root / "qc", sdrf=sdrf
                )
            self.assertTrue(result["success"])
            self.assertTrue((root / "qc/sdrf-changes.tsv").exists())

    def test_output_validator_exception_retains_evidence_and_records_failure(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            source = root / "run.mzML"
            source.touch()
            sdrf = root / "study.tsv"
            sdrf.write_text("comment[data file]\nrun.mzML\n")
            validator = MagicMock()
            validator.validate.side_effect = [
                report(sdrf),
                RuntimeError("validator failed"),
            ]
            outcome = FileOutcome(source, analyze(name=str(source)))
            with patch("prideqc.pipeline._analyze_file", return_value=outcome):
                result = Workflow(validator=validator).run(
                    [source], root / "qc", sdrf=sdrf
                )
            self.assertFalse(result["success"])
            self.assertTrue((root / "qc/sdrf-changes.tsv").exists())
            data = json.loads((root / "qc/sdrf-validation.json").read_text())
            self.assertIn("validator failed", data["refined_error"])


class PrideTests(unittest.TestCase):
    def test_current_client_is_loaded_lazily(self):
        module = types.ModuleType("pridepy.download.client")
        module.Client = MagicMock()
        repository = PrideRepository()
        module.Client.assert_not_called()
        with patch.dict(
            "sys.modules",
            {
                "pridepy": types.ModuleType("pridepy"),
                "pridepy.download": types.ModuleType("pridepy.download"),
                "pridepy.download.client": module,
            },
        ):
            self.assertIs(repository._get_client(), module.Client.return_value)
            self.assertIs(repository._get_client(), module.Client.return_value)
        module.Client.assert_called_once_with()

    def test_explicit_subset_is_deduplicated_and_delegated_with_checksums(self):
        client = MagicMock()
        client.download_file_by_name.side_effect = write_download
        with tempfile.TemporaryDirectory() as folder:
            destination = Path(folder) / "downloads"
            paths = PrideRepository(client).download_files(
                "pxd008644", ["run.mzML", "run.mzML"], destination
            )
            self.assertEqual(paths, [destination / "run.mzML"])
            args = client.download_file_by_name.call_args.kwargs
            self.assertEqual(args["accession"], "PXD008644")
            self.assertTrue(args["checksum_check"])
            self.assertFalse(args["skip_if_downloaded_already"])
            self.assertEqual(client.download_file_by_name.call_count, 1)
            receipt = json.loads((destination / "download-manifest.json").read_text())
            self.assertTrue(receipt["success"])
            self.assertGreater(receipt["files"][0]["size_bytes"], 0)
            self.assertEqual(len(list(destination.iterdir())), 2)

    def test_protocol_credentials_and_checksum_override_are_forwarded_without_serializing_secrets(self):
        client = MagicMock()
        client.download_file_by_name.side_effect = write_download
        with tempfile.TemporaryDirectory() as folder:
            PrideRepository(client).download_files(
                "PXD008644",
                ["run.mzML"],
                folder,
                options=DownloadOptions("globus", False, "50M"),
                username="private-user",
                password="secret-value",
            )
            args = client.download_file_by_name.call_args.kwargs
            self.assertEqual(args["protocol"], "globus")
            self.assertFalse(args["checksum_check"])
            self.assertEqual(args["password"], "secret-value")
            text = (Path(folder) / "download-manifest.json").read_text()
            self.assertNotIn("secret-value", text)
            self.assertNotIn("private-user", text)

    def test_incomplete_download_is_removed_and_not_reported_as_success(self):
        for failure in ("empty", "missing", "exception", "reported"):
            with self.subTest(failure=failure), tempfile.TemporaryDirectory() as folder:

                def fail(*, _failure=failure, **arguments):
                    path = Path(arguments["output_folder"]) / arguments["file_name"]
                    if _failure == "empty":
                        path.touch()
                    if _failure == "exception":
                        path.write_bytes(b"partial")
                        raise ValueError("https://user:secret@example.org")
                    if _failure == "reported":
                        path.write_bytes(b"partial")
                        return False

                client = MagicMock()
                client.download_file_by_name.side_effect = fail
                with self.assertRaises(RuntimeError) as raised:
                    PrideRepository(client).download_files(
                        "PXD008644", ["a.mzML", "b.mzML"], folder
                    )
                self.assertNotIn("secret", str(raised.exception))
                self.assertFalse((Path(folder) / "a.mzML").exists())
                receipt = json.loads(
                    (Path(folder) / "download-manifest.json").read_text()
                )
                self.assertFalse(receipt["success"])
                self.assertEqual(receipt["unprocessed"], ["b.mzML"])
                self.assertEqual(len(list(Path(folder).iterdir())), 1)

    def test_existing_files_unsafe_selection_and_collisions_fail_before_download(self):
        client = MagicMock()
        with tempfile.TemporaryDirectory() as folder:
            (Path(folder) / "old.mzML").touch()
            with self.assertRaises(FileExistsError):
                PrideRepository(client).download_files(
                    "PXD008644", ["a.mzML"], folder
                )
        for names in (
            ["../escape.raw"],
            ["a.raw", "A.raw"],
            ["not available"],
            ["download-manifest.json"],
            [],
        ):
            with self.subTest(names=names), self.assertRaises(ValueError):
                selected_files(names)
        client.download_file_by_name.assert_not_called()

    def test_project_workflow_downloads_sdrf_selection_then_runs_local_analysis(self):
        client = MagicMock()
        client.download_file_by_name.side_effect = write_download
        validator = MagicMock()
        validator.validate.side_effect = lambda path: report(path)
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            sdrf = root / "study.tsv"
            sdrf.write_text("comment[data file]\nrun.mzML\nrun.mzML\n")
            source = root / "downloads/run.mzML"
            with patch(
                "prideqc.pipeline._analyze_file",
                return_value=FileOutcome(source, analyze(name=str(source))),
            ):
                manifest = Workflow(validator=validator).run_project(
                    "PXD008644",
                    root / "qc",
                    download_directory=root / "downloads",
                    sdrf=sdrf,
                    repository=PrideRepository(client),
                )
            self.assertTrue(manifest["success"])
            self.assertEqual(client.download_file_by_name.call_count, 1)
            self.assertEqual(manifest["acquisition"]["accession"], "PXD008644")
            self.assertTrue((root / "qc/run.mzML.mzQC").exists())
            self.assertIn(
                "Orbitrap Fusion", (root / "qc/refined.sdrf.tsv").read_text()
            )

    def test_project_preflight_rejects_nested_outputs_and_unsupported_vendor_data_without_converter(self):
        repository = MagicMock()
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            with self.assertRaisesRegex(ValueError, "non-nested"):
                Workflow().run_project(
                    "PXD008644",
                    root,
                    download_directory=root / "downloads",
                    filenames=["a.mzML"],
                    repository=repository,
                )
            # Native vendor support is intentionally environment-dependent.
            # Force the unsupported-reader branch so this test remains valid
            # with the pyOpenMS 3.6 development build used by the container.
            with patch("prideqc.readers.PyOpenMSReader") as reader_type:
                reader_type.return_value.supports_direct.return_value = False
                with self.assertRaisesRegex(ValueError, "converter"):
                    Workflow().run_project(
                        "PXD008644",
                        root / "qc",
                        download_directory=root / "downloads",
                        filenames=["a.raw"],
                        repository=repository,
                    )
        repository.download_files.assert_not_called()

    def test_cli_fetch_and_remote_analysis_dispatch(self):
        with patch(
            "prideqc.pride.PrideRepository.download_files",
            return_value=[Path("downloads/a.mzML")],
        ) as download:
            with redirect_stdout(io.StringIO()):
                self.assertEqual(
                    main(
                        [
                            "fetch",
                            "PXD008644",
                            "--file",
                            "a.mzML",
                            "-o",
                            "downloads",
                        ]
                    ),
                    0,
                )
            self.assertEqual(download.call_args.args[:2], ("PXD008644", ["a.mzML"]))
        with patch(
            "prideqc.pipeline.Workflow.run_project",
            return_value={"files": [], "errors": [], "success": True},
        ) as run:
            with redirect_stdout(io.StringIO()):
                self.assertEqual(
                    main(
                        [
                            "analyze",
                            "--accession",
                            "PXD008644",
                            "--file",
                            "a.mzML",
                            "--download-dir",
                            "downloads",
                            "-o",
                            "qc",
                        ]
                    ),
                    0,
                )
            self.assertEqual(run.call_args.kwargs["filenames"], ["a.mzML"])


class InstalledDependencyTests(unittest.TestCase):
    def test_ci_requires_maintained_dependencies(self):
        if os.environ.get("PRIDE_QC_REQUIRE_INTEGRATION") == "1":
            for module in ("pridepy", "sdrf_pipelines"):
                self.assertIsNotNone(
                    importlib.util.find_spec(module),
                    f"Required CI dependency missing: {module}",
                )

    @unittest.skipUnless(
        importlib.util.find_spec("sdrf_pipelines"), "sdrf-pipelines not installed"
    )
    def test_installed_validator_returns_upstream_findings(self):
        from sdrf_pipelines.sdrf.sdrf import read_sdrf

        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "invalid.sdrf.tsv"
            path.write_text("source name\tcomment[data file]\nsample\trun.mzML\n")
            expected = read_sdrf(str(path)).validate_sdrf(
                template="ms-proteomics", skip_ontology=True
            )
            actual = SDRFPipelinesValidator().validate(path)
            self.assertTrue(
                expected, "Deliberately incomplete SDRF must have validation findings"
            )
            self.assertEqual(
                actual.issues,
                tuple(str(getattr(item, "message", item)) for item in expected),
            )
            self.assertFalse(actual.valid)

    @unittest.skipUnless(importlib.util.find_spec("pridepy"), "pridepy not installed")
    def test_installed_client_signature_accepts_download_adapter(self):
        from pridepy.download.client import Client

        # Enforce the real public API signature without a network transfer.
        client = create_autospec(Client, instance=True)
        client.download_file_by_name.side_effect = write_download
        with tempfile.TemporaryDirectory() as folder:
            paths = PrideRepository(client).download_files(
                "PXD008644", ["run.mzML"], folder
            )
            self.assertEqual(len(paths), 1)
