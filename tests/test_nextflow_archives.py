"""Regression tests for archive discovery and fail-closed run ledgers."""
from __future__ import annotations

import gzip
import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def _load(name: str, filename: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / "nextflow/bin" / filename)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load test module: {filename}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


arch = _load("nf_arch", "archive_units.py")
final = _load("nf_final", "finalize_accession.py")


class NextflowArchiveTests(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.root = Path(temp.name)

    def _zip(self, name, members):
        p = self.root / name
        with zipfile.ZipFile(p, "w") as zf:
            for member, content in members.items():
                zf.writestr(member, content)
        return p

    def test_big_zip_creates_multiple_independent_runs(self):
        archive = self._zip("all_runs.zip", {
            "experiments/one.raw": b"abc", "experiments/two.raw": b"xyz",
            "notes/readme.txt": b"information",
        })
        result = arch.expand(archive, self.root / "out")
        self.assertEqual(result["logical_run_count"], 2)
        self.assertEqual([r["input_basename"] for r in result["units"]], ["one.raw", "two.raw"])
        for item in result["units"]:
            self.assertEqual(len(list((self.root / "out" / item["unit_id"]).glob("*.raw"))), 1)

    def test_raw_zip_one_run(self):
        result = arch.expand(self._zip("one.raw.zip", {"one.raw": b"raw"}), self.root / "out")
        self.assertEqual(result["logical_run_count"], 1)
        self.assertEqual(result["units"][0]["format"], "raw")

    def test_d_zip_preserves_vendor_folder(self):
        result = arch.expand(self._zip("run.d.zip", {
            "run.d/analysis.tdf": b"tdf", "run.d/analysis.tdf_bin": b"binary",
        }), self.root / "out")
        self.assertEqual(result["logical_run_count"], 1)
        unit = self.root / "out" / result["units"][0]["unit_id"] / "run.d"
        self.assertTrue((unit / "analysis.tdf_bin").is_file())

    def test_nested_zip_two_run_expansion(self):
        inner = self._zip("run.raw.zip", {"run.raw": b"raw"})
        outer = self._zip("group.zip", {"nested/run.raw.zip": inner.read_bytes()})
        result = arch.expand(outer, self.root / "out")
        self.assertEqual(result["logical_run_count"], 1)

    def test_mzml_gz_is_one_pass_through_run(self):
        source = self.root / "data.mzML.gz"
        source.write_bytes(gzip.compress(b"<mzML>test</mzML>"))
        result = arch.expand(source, self.root / "out")
        self.assertEqual(result["logical_run_count"], 1)
        self.assertEqual(result["units"][0]["format"], "mzml_gz")

    def test_corrupt_gzip_fails_closed(self):
        source = self.root / "bad.mzML.gz"
        source.write_bytes(b"not gzip")
        with self.assertRaises(OSError):
            arch.expand(source, self.root / "out")

    def test_gzip_raw_is_materialized(self):
        source = self.root / "run.raw.gz"
        source.write_bytes(gzip.compress(b"binary raw data"))
        result = arch.expand(source, self.root / "out")
        self.assertEqual(result["logical_run_count"], 1)
        self.assertEqual(result["units"][0]["input_basename"], "run.raw")

    def test_zip_containing_individual_raw_gzip(self):
        source = self._zip("bundle.zip", {"nested/run.raw.gz": gzip.compress(b"raw")})
        result = arch.expand(source, self.root / "out")
        self.assertEqual(result["logical_run_count"], 1)
        self.assertEqual(result["units"][0]["input_basename"], "run.raw")

    def test_zip_slip_rejected(self):
        archive = self._zip("malicious.zip", {"../evil.raw": b"x"})
        with self.assertRaisesRegex(arch.ArchiveError, "Unsafe"):
            arch.expand(archive, self.root / "out")
        self.assertFalse((self.root / "evil.raw").exists())

    def test_empty_archive_fails_closed(self):
        archive = self._zip("report.zip", {"notes.txt": b"not data"})
        with self.assertRaisesRegex(arch.ArchiveError, "No supported"):
            arch.expand(archive, self.root / "out")

    def test_duplicate_member_rejected(self):
        archive = self._zip("duplicate.zip", {"A.raw": b"a", "a.raw": b"b"})
        with self.assertRaisesRegex(arch.ArchiveError, "Duplicate case-insensitive ZIP member"):
            arch.expand(archive, self.root / "out")

    def test_oversize_rejected(self):
        archive = self._zip("large.zip", {"a.raw": b"123456"})
        with self.assertRaisesRegex(arch.ArchiveError, "Invalid ZIP member size"):
            arch.expand(archive, self.root / "out", max_bytes=3)

    def test_archive_preparation_emits_failure_sentinel(self):
        source = self._zip("no-raw.zip", {"README.txt": b"no raw here"})
        task = self.root / "task.json"
        task.write_text(json.dumps({
            "accession": "PXD041271", "task_id": 99, "archive_file": source.name,
            "expected_runs": 1, "source_manifest": "test.tsv",
        }))
        command = [sys.executable, str(ROOT / "nextflow/bin/prepare_archive.py"),
                   "--task", str(task), "--protocol", "ftp",
                   "--local-source", str(source), "--outdir", str(self.root / "prepared")]
        run = subprocess.run(command, capture_output=True, text=True, check=False)
        self.assertEqual(run.returncode, 0, run.stderr)
        receipt = json.loads((self.root / "prepared/prepare-status.json").read_text())
        self.assertEqual(receipt["status"], "failed")
        units = list((self.root / "prepared/units").glob("*/unit.json"))
        self.assertEqual(len(units), 1)
        self.assertEqual(json.loads(units[0].read_text())["discovery_status"], "failed")

    def test_failed_discovery_generates_per_unit_qc_failure(self):
        unit = self.root / "unit"
        unit.mkdir()
        (unit / "unit.json").write_text(json.dumps({
            "accession": "PXD041271", "unit_id": "failed-archive-1",
            "task_id": 1, "source_archive": "broken.zip", "expected_file_tasks": 1,
            "discovery_status": "failed", "failure_stage": "archive_preparation",
            "error_type": "BadZipFile",
        }))
        cmd = [sys.executable, str(ROOT / "nextflow/bin/analyze_unit.py"),
               "--unit", str(unit), "--output", str(self.root / "result")]
        run = subprocess.run(cmd, capture_output=True, text=True, check=False)
        self.assertEqual(run.returncode, 0, run.stderr)
        result = json.loads((self.root / "result/result.json").read_text())
        self.assertEqual(result["failure_stage"], "archive_preparation")
        self.assertEqual(result["status"], "failed")

    def test_one_failed_unit_blocks_study_and_creates_ledger(self):
        root = self.root / "results"
        root.mkdir()
        statuses = [("unit1", "success"), ("unit2", "failed")]
        for identifier, status in statuses:
            item = root / identifier
            item.mkdir()
            (item / "result.json").write_text(json.dumps({
                "unit_id": identifier, "task_id": 1, "expected_file_tasks": 1,
                "accession": "PXD041271", "source_archive": "bundle.zip",
                "status": status, "failure_stage": "prideqc_analysis" if status == "failed" else None,
            }))
        summary = final.finalize("PXD041271", 2, root, self.root / "out")
        self.assertEqual(summary["status"], "incomplete")
        self.assertEqual(summary["failed_runs"], 1)
        self.assertFalse((self.root / "out/_study").exists())
        self.assertEqual(len((self.root / "out/_audit/file-status.tsv").read_text().splitlines()), 3)


if __name__ == "__main__":
    unittest.main()
