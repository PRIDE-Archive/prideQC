from __future__ import annotations

import csv
import importlib.util
import tempfile
import unittest
from pathlib import Path
from unittest import mock

SCRIPT = (
    Path(__file__).parents[1]
    / "scripts"
    / "slurm"
    / "prepare_prideqc_accession_manifest.py"
)
SPEC = importlib.util.spec_from_file_location("prepare_prideqc_accession_manifest", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
manifest = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(manifest)


def raw(name: str, size: int = 100) -> dict[str, object]:
    return {
        "fileName": name,
        "fileSizeBytes": size,
        "fileCategory": {"value": "RAW"},
        "publicFileLocations": [
            {"name": "FTP Protocol", "value": f"ftp://example.invalid/{name}"},
            {"name": "Aspera Protocol", "value": f"aspera://example.invalid/{name}"},
        ],
    }


class PrideAccessionManifestTests(unittest.TestCase):
    def test_page_records_accepts_current_v3_bare_list(self) -> None:
        records, number, pages = manifest._page_records([raw("1.raw")], "PXD000138")
        self.assertEqual(records[0]["fileName"], "1.raw")
        self.assertIsNone(number)
        self.assertIsNone(pages)

    def test_project_files_pages_using_total_records_header(self) -> None:
        payloads = [
            ([raw("1.raw")], {"total_records": "2"}),
            ([raw("2.raw")], {"total_records": "2"}),
        ]
        calls = 0

        def fake_request(url: str, timeout: int = 300):
            nonlocal calls
            response = payloads[calls]
            calls += 1
            return response

        with mock.patch.object(manifest, "request_json_with_headers", fake_request):
            records = manifest.pride_project_files("PXD000138")

        self.assertEqual([record["fileName"] for record in records], ["1.raw", "2.raw"])
        self.assertEqual(calls, 2)

    def test_explicit_smoke_subset_keeps_schema_and_leaves_sdrf_blank(self) -> None:
        records = [
            raw("1.raw", 10),
            raw("2.raw", 20),
            raw("10.raw", 30),
            {
                "fileName": "search.mzid",
                "fileSizeBytes": 5,
                "fileCategory": {"value": "RESULT"},
                "publicFileLocations": [],
            },
        ]
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            output = root / "smoke.tsv"
            meta = root / "smoke.meta.txt"
            inventory = root / "smoke.inventory.tsv"
            manifest.write_outputs(
                "PXD000138",
                records,
                ["10.raw", "1.raw"],
                output,
                meta,
                inventory,
                manifest.PRIDE_API,
            )

            with output.open(encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(list(rows[0]), manifest.MANIFEST_FIELDS)
            self.assertEqual([row["archive_file"] for row in rows], ["10.raw", "1.raw"])
            self.assertTrue(all(row["sdrf_path"] == "" for row in rows))
            self.assertTrue(all(row["sdrf_data_file"] == "" for row in rows))
            self.assertTrue(all(row["needs_file_map"] == "0" for row in rows))
            self.assertEqual(rows[0]["file_uri"], "ftp://example.invalid/10.raw")
            text = meta.read_text(encoding="utf-8")
            self.assertIn("manifest_mode=pride-project-metadata", text)
            self.assertIn("project_file_count=4", text)
            self.assertIn("raw_file_count=3", text)
            self.assertIn("selected_file_count=2", text)
            self.assertIn("trusted_sdrf_used=0", text)

    def test_full_raw_selection_is_naturally_sorted_and_excludes_results(self) -> None:
        records = [
            raw("10.raw"),
            raw("2.raw"),
            raw("1.raw"),
            {"fileName": "x.mgf", "fileCategory": {"value": "PEAK"}},
        ]
        raw_records, selected = manifest.select_raw_files(records, [])
        self.assertEqual(len(raw_records), 3)
        self.assertEqual(
            [record["fileName"] for record in selected], ["1.raw", "2.raw", "10.raw"]
        )

    def test_requested_non_raw_or_missing_file_is_rejected(self) -> None:
        records = [raw("1.raw"), {"fileName": "x.mgf", "fileCategory": {"value": "PEAK"}}]
        with self.assertRaisesRegex(ValueError, "not current PRIDE RAW files"):
            manifest.select_raw_files(records, ["x.mgf"])


if __name__ == "__main__":
    unittest.main()
