from __future__ import annotations

import importlib.util
import tempfile
import unittest
from pathlib import Path

SCRIPT = (
    Path(__file__).parents[1]
    / "scripts"
    / "slurm"
    / "prepare_ground_truth_task_subset.py"
)
SPEC = importlib.util.spec_from_file_location("prepare_ground_truth_task_subset", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
subset = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(subset)


class GroundTruthTaskSubsetTests(unittest.TestCase):
    def test_resolve_and_write_uses_unix_line_endings(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            root = Path(tmp_dir)
            gt = root / "gt.tsv"
            manifest = root / "manifest.tsv"
            output = root / "tasks.tsv"

            gt.write_text(
                "gt_id\tpxd_accession\traw_file\n"
                "GT001\tPXD000001\ta.raw\n"
                "GT002\tPXD000002\tb.raw\n",
                encoding="utf-8",
            )
            manifest.write_text(
                "task_id\tpxd_accession\tarchive_file\n"
                "7\tPXD000001\ta.raw\n"
                "11\tPXD000002\tb.raw\n",
                encoding="utf-8",
            )

            rows = subset.resolve_tasks(gt, manifest)
            subset.write_task_subset(output, rows)

            self.assertEqual(subset.slurm_array(rows), "7,11")
            self.assertEqual([row["task_id"] for row in rows], ["7", "11"])

            data = output.read_bytes()
            self.assertNotIn(b"\r\n", data)
            self.assertNotIn(b"\r", data)
            self.assertEqual(data.count(b"\n"), 3)

    def test_duplicate_physical_task_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            root = Path(tmp_dir)
            gt = root / "gt.tsv"
            manifest = root / "manifest.tsv"

            gt.write_text(
                "gt_id\tpxd_accession\traw_file\n"
                "GT001\tPXD000001\ta.raw\n"
                "GT002\tPXD000001\ta.raw\n",
                encoding="utf-8",
            )
            manifest.write_text(
                "task_id\tpxd_accession\tarchive_file\n"
                "7\tPXD000001\ta.raw\n",
                encoding="utf-8",
            )

            with self.assertRaisesRegex(ValueError, "duplicate physical task IDs"):
                subset.resolve_tasks(gt, manifest)


if __name__ == "__main__":
    unittest.main()
