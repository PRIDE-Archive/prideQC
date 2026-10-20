"""Pilot isolation: RAW must stay node-local, subset never becomes a complete study."""
from __future__ import annotations

import importlib
import json
import sys
import tempfile
import unittest
import zipfile
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "nextflow" / "bin"))

analyze_unit = importlib.import_module("analyze_unit")
finalize_accession = importlib.import_module("finalize_accession")
source_qc = importlib.import_module("source_qc")


class ScratchPilotTests(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.root = Path(tmp.name)
        self.task = self.root / "task.json"

    def _task(self, archive: Path):
        self.task.write_text(json.dumps({
            "accession": "PXD041271", "task_id": 17, "archive_file": archive.name,
            "expected_runs": 1, "source_manifest": "pilot.tsv",
        }))

    def _run(self, archive: Path):
        self._task(archive)
        with patch.object(source_qc.Path, "cwd", return_value=self.root):
            return source_qc.run(self.task, "ftp", self.root / "node-local-raw",
                                 self.root / "ready", source=archive)

    def test_one_run_only_qc_outputs_persist(self):
        archive = self.root / "one.mzML"
        archive.write_bytes(b"not a real MS run - external QC is stubbed")

        def mocked_analysis(cmd, **_kwargs):
            output = Path(cmd[cmd.index("-o") + 1])
            output.mkdir()
            (output / "manifest.json").write_text(json.dumps({"files": [{"status": "success"}]}))
            (output / "one.mzML.mzQC").write_text(json.dumps({"mzQC": {
                "version": "1.0.0", "runQualities": [{"qualityMetrics": [
                    {"accession": "MS:4000059", "name": "number of MS1 spectra", "value": 1},
                    {"accession": "MS:4000060", "name": "number of MS2 spectra", "value": 2},
                ]}]
            }}))
            (output / "one.mzML.obo").write_text("format-version: 1.2\n")
            return SimpleNamespace(returncode=0, stdout="QC success", stderr="")

        with patch.object(analyze_unit.subprocess, "run", side_effect=mocked_analysis):
            result = self._run(archive)
        self.assertEqual(result["status"], "success")
        self.assertEqual(result["requested_estimators"], {"mass_error": False,
                                                         "mass_shifts": False})
        self.assertFalse((self.root / "node-local-raw").exists())
        staged = list((self.root / "ready").rglob("*"))
        self.assertFalse(any(p.name.lower().endswith((".raw", ".mzml", ".mzml.gz")) for p in staged))
        self.assertEqual(len(list((self.root / "ready").rglob("*.mzQC"))), 1)
        receipt_root = self.root / "ready" / "units"
        status = finalize_accession.finalize("PXD041271", 1, receipt_root, self.root / "final")
        self.assertEqual(status["status"], "partial")
        self.assertEqual(status["study_merge_performed"], False)
        self.assertFalse((self.root / "final" / "_study").exists())

    def test_zero_spectra_qc_is_not_success_even_with_zero_exit(self):
        archive = self.root / "empty.mzML"
        archive.write_bytes(b"stub empty mzML")

        def mocked_analysis(cmd, **_kwargs):
            output = Path(cmd[cmd.index("-o") + 1])
            output.mkdir()
            (output / "manifest.json").write_text(json.dumps({"files": [{"status": "success"}]}))
            (output / "empty.mzML.mzQC").write_text(json.dumps({"mzQC": {
                "runQualities": [{"qualityMetrics": [
                    {"accession": "MS:4000059", "name": "number of MS1 spectra", "value": 0},
                    {"accession": "MS:4000060", "name": "number of MS2 spectra", "value": 0},
                ]}]
            }}))
            (output / "empty.mzML.obo").write_text("format-version: 1.2\n")
            return SimpleNamespace(returncode=0, stdout="QC ran", stderr="")

        with patch.object(analyze_unit.subprocess, "run", side_effect=mocked_analysis):
            result = self._run(archive)
        self.assertEqual(result["status"], "failed")
        self.assertEqual(result["error_type"], "EmptySpectrumRun")
        self.assertEqual(result["failure_stage"], "qc_validation")
        self.assertFalse((self.root / "node-local-raw").exists())
        finalized = finalize_accession.finalize("PXD041271", 1, self.root / "ready" / "units", self.root / "final")
        self.assertEqual(finalized["status"], "incomplete")
        self.assertEqual(finalized["successful_runs"], 0)

    def test_ms3_only_is_not_rejected_as_empty(self):
        mzqc = self.root / "ms3.mzQC"
        mzqc.write_text(json.dumps({"mzQC": {"runQualities": [{"qualityMetrics": [
            {"accession": "MS:4000059", "name": "number of MS1 spectra", "value": 0},
            {"accession": "MS:4000060", "name": "number of MS2 spectra", "value": 0},
            {"accession": "QCPRIDE:EXAMPLE", "name": "NumberOfSpectra_MS3", "value": 4},
        ]}]}}))
        self.assertEqual(analyze_unit.require_observed_spectra(mzqc), 4)

    def test_bad_spectrum_count_is_rejected(self):
        mzqc = self.root / "invalid.mzQC"
        mzqc.write_text(json.dumps({"mzQC": {"runQualities": [{"qualityMetrics": [
            {"accession": "MS:4000059", "name": "number of MS1 spectra", "value": -1},
            {"accession": "MS:4000060", "name": "number of MS2 spectra", "value": 4},
        ]}]}}))
        with self.assertRaisesRegex(ValueError, "Invalid or duplicate"):
            analyze_unit.require_observed_spectra(mzqc)

    def test_multi_run_archive_fails_closed_without_qc_or_raw(self):
        archive = self.root / "many.zip"
        with zipfile.ZipFile(archive, "w") as zf:
            zf.writestr("one.raw", b"raw 1")
            zf.writestr("two.raw", b"raw 2")
        result = self._run(archive)
        self.assertEqual(result["status"], "failed")
        self.assertEqual(result["error_type"], "MultipleLogicalRunsUnsupported")
        self.assertFalse((self.root / "node-local-raw").exists())
        self.assertEqual(len(list((self.root / "ready").rglob("*.raw"))), 0)
        self.assertEqual(len(list((self.root / "ready").rglob("*.mzQC"))), 0)
        status = finalize_accession.finalize("PXD041271", 1, self.root / "ready" / "units", self.root / "final")
        self.assertEqual(status["status"], "incomplete")
        self.assertFalse((self.root / "final" / "_study").exists())

    def test_bad_archive_has_failure_receipt_and_cleans_scratch(self):
        archive = self.root / "empty.zip"
        with zipfile.ZipFile(archive, "w") as zf:
            zf.writestr("README.txt", b"not mass spectrometry")
        result = self._run(archive)
        self.assertEqual(result["status"], "failed")
        self.assertFalse((self.root / "node-local-raw").exists())
        self.assertEqual(len(list((self.root / "ready").rglob("result.json"))), 1)

    def test_rejects_scratch_outside_task_cwd(self):
        archive = self.root / "one.raw"
        archive.write_bytes(b"raw")
        self._task(archive)
        with patch.object(source_qc.Path, "cwd", return_value=self.root):
            with self.assertRaisesRegex(ValueError, "task-local cwd"):
                source_qc.run(self.task, "ftp", self.root / "nested" / "raw",
                              self.root / "ready", source=archive)


if __name__ == "__main__":
    unittest.main()
