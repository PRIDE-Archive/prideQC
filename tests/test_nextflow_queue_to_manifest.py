"""Batch cutting from the inventory queue: deterministic, bounded, accepted by plan_manifest."""
from __future__ import annotations

import csv
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

BIN = Path(__file__).resolve().parents[1] / "nextflow/bin"


def load(name, file):
    spec = importlib.util.spec_from_file_location(name, BIN / file)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


Q = load("prideqc_nf_queue", "queue_to_manifest.py")
PLAN = load("prideqc_nf_plan_q", "plan_manifest.py")


class QueueToManifestTests(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.root = Path(tmp.name)
        rows = []
        for acc, sizes in (("PXD000001", [5e9, 1e8, 2e8]), ("PXD000002", [3e8, 4e8]), ("PXD000003", [7e9])):
            for i, size in enumerate(sizes):
                rows.append({"accession": acc, "file_name": f"{acc}_{i}.raw", "file_format": "thermo_raw",
                             "file_size_bytes": str(int(size)), "checksum": "",
                             "memory_gb": str(PLAN.memory_tier_gb(int(size))), "decision": "QUEUE_NATIVE_RAW"})
        rows.append({"accession": "PXD000002", "file_name": "x.mzML", "file_format": "mzml",
                     "file_size_bytes": "10", "checksum": "", "memory_gb": "16",
                     "decision": "QUEUE_MZML_NO_VENDOR"})
        self.queue = self.root / "queue.tsv"
        with self.queue.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)

    def test_large_first_then_round_robin_and_plan_manifest_accepts(self):
        index = self.root / "index.jsonl"
        index.write_text(json.dumps({"accession": "PXD000002", "source_file": "PXD000002_0.raw",
                                     "role": "canonical"}) + "\n")
        out = self.root / "batch"
        summary = Q.build(self.queue, out, max_files=4, max_gb=100, large_files=1, processed_index=index)
        with (out / "file-manifest.tsv").open() as handle:
            names = [r["archive_file"] for r in csv.DictReader(handle, delimiter="\t")]
        self.assertIn("PXD000003_0.raw", names)          # largest first
        self.assertNotIn("PXD000002_0.raw", names)       # already canonical
        self.assertNotIn("x.mzML", names)                # mzML lane not qualified
        self.assertEqual(summary["files"], 4)
        self.assertFalse(summary["complete_study_asserted"])
        plan = PLAN.prepare(out / "file-manifest.tsv", out / "report-manifest.tsv", self.root / "plan", set())
        self.assertEqual(plan["tasks"], 4)
        task = json.loads(next((self.root / "plan/tasks").glob("PXD000003__*.json")).read_text())
        self.assertEqual(task["memory_gb"], 128)
        again = Q.build(self.queue, self.root / "batch2", max_files=4, max_gb=100, large_files=1,
                        processed_index=index)
        self.assertEqual(again["file_manifest_sha256"], summary["file_manifest_sha256"])  # deterministic

    def test_limits_and_refusals(self):
        summary = Q.build(self.queue, self.root / "b", max_files=10, max_gb=1.0)
        self.assertLessEqual(summary["total_gb"], 1.0)
        one = Q.build(self.queue, self.root / "c", max_files=10, max_gb=100, per_accession=1)
        self.assertEqual(one["files"], 3)
        with self.assertRaisesRegex(ValueError, "Refusing existing"):
            Q.build(self.queue, self.root / "b", max_files=1, max_gb=1)
        with self.assertRaisesRegex(ValueError, "Empty batch"):
            Q.build(self.queue, self.root / "d", max_files=1, max_gb=0.01)

    def test_file_cap_and_per_tier(self):
        capped = Q.build(self.queue, self.root / "e", max_files=10, max_gb=100, max_file_gb=6,
                         large_files=1)
        self.assertEqual(capped["largest_gb"], 5.0)      # 7 GB file above the cap is never taken
        tiered = Q.build(self.queue, self.root / "f", max_files=2, max_gb=100, per_tier=1)
        self.assertEqual(tiered["memory_gb_tiers"], {"128": 1, "16": 1})  # one per tier before round-robin

    def test_download_failures_retry_only_on_request(self):
        index = self.root / "index.jsonl"
        index.write_text("\n".join(json.dumps(r) for r in (
            {"accession": "PXD000001", "source_file": "PXD000001_1.raw", "role": "failed",
             "failure_stage": "archive_preparation"},
            {"accession": "PXD000001", "source_file": "PXD000001_2.raw", "role": "failed",
             "failure_stage": "prideqc_analysis"})) + "\n")

        def names(out):
            with (out / "file-manifest.tsv").open() as handle:
                return {r["archive_file"] for r in csv.DictReader(handle, delimiter="\t")}
        Q.build(self.queue, self.root / "g", max_files=10, max_gb=100, processed_index=index)
        self.assertFalse({"PXD000001_1.raw", "PXD000001_2.raw"} & names(self.root / "g"))
        summary = Q.build(self.queue, self.root / "h", max_files=10, max_gb=100, processed_index=index,
                          retry_download_failures=True)
        self.assertIn("PXD000001_1.raw", names(self.root / "h"))      # download failure retried
        self.assertNotIn("PXD000001_2.raw", names(self.root / "h"))   # reader failure stays excluded
        self.assertTrue(summary["rules"]["retry_download_failures"])

    def test_accession_cap_fills_chosen_studies_and_retry_comes_first(self):
        index = self.root / "index.jsonl"
        index.write_text(json.dumps({"accession": "PXD000002", "source_file": "PXD000002_1.raw",
                                     "role": "failed", "failure_stage": "archive_preparation"}) + "\n")
        summary = Q.build(self.queue, self.root / "i", max_files=10, max_gb=100, max_accessions=1,
                          processed_index=index, retry_download_failures=True)
        self.assertEqual(summary["accessions"], 1)
        with (self.root / "i" / "file-manifest.tsv").open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual({r["pxd_accession"] for r in rows}, {"PXD000002"})   # retried study chosen first
        self.assertEqual(len(rows), 2)                                        # and filled before others


if __name__ == "__main__":
    unittest.main()
