"""Regression: finalizers must see each QC receipt and publish partial audits.

No RAW inputs, external services, or native MS dependencies are required.
"""
from __future__ import annotations

import importlib.util
import json
import shutil
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MODULE = ROOT / "nextflow/modules/local/finalize/main.nf"
FINALIZER = ROOT / "nextflow/bin/finalize_accession.py"


def _finalizer():
    spec = importlib.util.spec_from_file_location("prideqc_pilot_finalizer", FINALIZER)
    if spec is None or spec.loader is None:
        raise RuntimeError("Missing finalizer script")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


class FinalizerReceiptHandoffTests(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.root = Path(temp.name)

    def _run_finalization(self, status: str) -> dict:
        producer = self.root / "source_task" / "ready" / "units" / "archive-000000001-unit-00000"
        producer.mkdir(parents=True)
        receipt = {
            "accession": "PXD036869", "task_id": 1, "expected_file_tasks": 1,
            "unit_id": producer.name, "source_archive": "one.raw", "status": status,
            "failure_stage": None if status == "success" else "prideqc_analysis",
            "error_type": None if status == "success" else "PrideqcNonzeroExit",
        }
        (producer / "result.json").write_text(json.dumps(receipt), encoding="utf-8")
        staged = self.root / "finalizer_task" / "runs" / producer.name
        staged.parent.mkdir(parents=True)
        # Simulate Nextflow stageInMode copy for QC-only output directories.
        shutil.copytree(producer, staged)
        out = self.root / "finalizer_task" / "final"
        result = _finalizer().finalize("PXD036869", 1, staged.parent, out)
        self.assertTrue((out / "_audit" / "accession-status.json").is_file())
        self.assertTrue((out / "_audit" / "file-status.tsv").is_file())
        self.assertEqual(result["observed_receipts"], 1)
        self.assertEqual(result["source_archive_tasks_observed"], 1)
        self.assertFalse(result["study_merge_performed"])
        self.assertFalse((out / "_study").exists())
        return result

    def test_successful_single_run_is_partial_not_study_complete(self):
        result = self._run_finalization("success")
        self.assertEqual(result["status"], "partial")
        self.assertEqual(result["successful_runs"], 1)
        self.assertEqual(result["failed_runs"], 0)
        self.assertEqual(result["failures"], [])

    def test_failed_qc_receipt_keeps_accession_incomplete(self):
        result = self._run_finalization("failed")
        self.assertEqual(result["status"], "incomplete")
        self.assertEqual(result["successful_runs"], 0)
        self.assertEqual(result["failed_runs"], 1)
        self.assertIn("unit_failed:archive-000000001-unit-00000", result["failures"])

    def test_symlink_staging_is_reason_for_zero_receipts(self):
        producer = self.root / "producer"
        producer.mkdir()
        (producer / "result.json").write_text("{}", encoding="utf-8")
        stage = self.root / "runs"
        stage.mkdir()
        (stage / "unit").symlink_to(producer, target_is_directory=True)
        self.assertEqual(list(stage.rglob("result.json")), [])
        self.assertEqual(len(list((stage / "unit").rglob("result.json"))), 1)

    def test_process_copies_qc_directories_and_publishes_both_audits(self):
        module = MODULE.read_text(encoding="utf-8")
        self.assertIn("stageInMode 'copy'", module)
        self.assertIn("${params.outdir}/${accession}/_audit", module)
        self.assertIn("path('final/_audit/accession-status.json')", module)
        self.assertIn("path('final/_audit/file-status.tsv')", module)
        self.assertIn("filename.tokenize('/').last()", module)
        self.assertNotIn("pattern: 'final/*'", module)


if __name__ == "__main__":
    unittest.main()
