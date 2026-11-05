"""Offline unit checks for the new Nextflow input and Mongo exporter contracts."""

import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def load(name, file):
    spec = importlib.util.spec_from_file_location(name, file)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


planner = load("prideqc_nf_plan", ROOT / "nextflow/bin/plan_manifest.py")


class NextflowManifestTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.files = self.root / "files.tsv"
        self.reports = self.root / "reports.tsv"
        self.files.write_text(
            "task_id\tpxd_accession\tarchive_file\tsdrf_path\tneeds_file_map\n"
            "506\tPXD041271\ta.raw\t\t0\n507\tPXD041271\tb.raw\t\t0\n"
        )
        self.reports.write_text("report_id\tpxd_accession\texpected_file_tasks\n4\tPXD041271\t2\n")

    def test_plan_and_expected_counts(self):
        info = planner.prepare(self.files, self.reports, self.root / "out", set())
        self.assertEqual(info["tasks"], 2)
        first = json.loads((self.root / "out/tasks/PXD041271__000000506.json").read_text())
        self.assertEqual(first["expected_runs"], 2)
        self.assertEqual(first["archive_file"], "a.raw")

    def test_memory_tier_from_file_size(self):
        self.files.write_text(
            "task_id\tpxd_accession\tarchive_file\tsdrf_path\tneeds_file_map\tfile_size_bytes\n"
            "506\tPXD041271\ta.raw\t\t0\t999999999\n507\tPXD041271\tb.raw\t\t0\t2614135068\n"
        )
        planner.prepare(self.files, self.reports, self.root / "out", set())
        a = json.loads((self.root / "out/tasks/PXD041271__000000506.json").read_text())
        b = json.loads((self.root / "out/tasks/PXD041271__000000507.json").read_text())
        self.assertEqual((a["memory_gb"], b["memory_gb"]), (16, 32))
        self.assertEqual(b["file_size_bytes"], 2614135068)
        self.assertEqual(planner.memory_tier_gb(None), 32)
        self.assertEqual(planner.memory_tier_gb(5_000_000_000), 64)
        self.assertEqual(planner.memory_tier_gb(9_000_000_000), 128)

    def test_missing_size_column_defaults_to_32gb_and_bad_size_fails(self):
        planner.prepare(self.files, self.reports, self.root / "out", set())
        a = json.loads((self.root / "out/tasks/PXD041271__000000506.json").read_text())
        self.assertEqual((a["file_size_bytes"], a["memory_gb"]), (None, 32))
        self.files.write_text(
            "task_id\tpxd_accession\tarchive_file\tsdrf_path\tneeds_file_map\tfile_size_bytes\n"
            "506\tPXD041271\ta.raw\t\t0\t-5\n507\tPXD041271\tb.raw\t\t0\t10\n"
        )
        with self.assertRaisesRegex(ValueError, "invalid file_size_bytes"):
            planner.prepare(self.files, self.reports, self.root / "out2", set())

    def test_wrong_report_count_fails(self):
        self.reports.write_text(self.reports.read_text().replace("\t2\n", "\t3\n"))
        with self.assertRaisesRegex(ValueError, "Count mismatch"):
            planner.prepare(self.files, self.reports, self.root / "out", set())

    def test_duplicate_and_traversal_fail_closed(self):
        self.files.write_text(self.files.read_text().replace("b.raw", "../evil.raw"))
        with self.assertRaisesRegex(ValueError, "Unsafe filename"):
            planner.prepare(self.files, self.reports, self.root / "out", set())

    def test_sdrf_rows_explicitly_unsupported(self):
        self.files.write_text(
            self.files.read_text().replace("a.raw\t\t0", "a.raw\tdata.sdrf.tsv\t0")
        )
        with self.assertRaisesRegex(ValueError, "SDRF mapping not supported"):
            planner.prepare(self.files, self.reports, self.root / "out", set())


if __name__ == "__main__":
    unittest.main()
