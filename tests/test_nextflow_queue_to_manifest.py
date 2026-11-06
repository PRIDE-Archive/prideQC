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


if __name__ == "__main__":
    unittest.main()
