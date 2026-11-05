"""PRIDE-wide inventory: listing verification, representation choice, resume. Offline."""
from __future__ import annotations

import csv
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "nextflow/bin/pride_inventory.py"
SPEC = importlib.util.spec_from_file_location("prideqc_nf_inventory", SCRIPT)
INV = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(INV)


def record(name, size=1000, checksum=""):
    return {"fileName": name, "fileSizeBytes": size, "checksum": checksum,
            "fileCategory": {"value": "RAW"}}


class FakeApi:
    def __init__(self, projects, counts=None):
        self.projects = projects
        self.counts = counts or {}
        self.calls = []

    def __call__(self, url):
        self.calls.append(url)
        path = url.split("/v3/", 1)[1]
        if path.startswith("projects?"):
            page = int(path.split("page=")[1].split("&")[0])
            items = [{"accession": a} for a in sorted(self.projects)]
            return items[page * INV.PAGE:(page + 1) * INV.PAGE]
        accession = path.split("/")[1]
        files = self.projects[accession]
        if path.endswith("/files/count"):
            return self.counts.get(accession, len(files))
        page = int(path.split("page=")[1])
        return files[page * INV.PAGE:(page + 1) * INV.PAGE]


def rows(out):
    with (out / "inventory.tsv").open() as handle:
        return {r["file_name"]: r for r in csv.DictReader(handle, delimiter="\t")}


class InventoryTests(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.root = Path(tmp.name)

    def test_formats_and_representation_choice(self):
        api = FakeApi({"PXD000001": [
            record("a.raw", 2_000_000_000), record("a.mzML"),            # RAW preferred
            record("b.mzML.gz"),                                           # no vendor: mzML lane
            record("c.wiff"), record("c.wiff.scan"), record("c.mzML"),    # WIFF: policy pending
            record("d.zip"), record("e.mzXML"), record("f.d.zip"),
            record("g.raw.gz"), record("search.mzid.gz"), record("peaks.mgf"),
            record("bad name.raw"), record("zero.raw", 0),
        ]})
        summary = INV.build(["PXD000001"], api, self.root / "out")
        r = rows(self.root / "out")
        self.assertEqual(r["a.raw"]["decision"], "QUEUE_NATIVE_RAW")
        self.assertEqual(r["a.raw"]["memory_gb"], "32")
        self.assertEqual(r["a.mzML"]["reason"], "SKIP_MZML_VENDOR_RAW_PREFERRED")
        self.assertEqual(r["a.mzML"]["representation_evidence"], "stem_match_unverified")
        self.assertEqual(r["b.mzML.gz"]["decision"], "QUEUE_MZML_NO_VENDOR")
        self.assertEqual(r["c.wiff"]["reason"], "SKIP_SCIEX_WIFF_NO_QUALIFIED_READER")
        self.assertEqual(r["c.wiff.scan"]["reason"], "SKIP_SCIEX_WIFF_SCAN_COMPANION")
        self.assertEqual(r["c.mzML"]["reason"], "SKIP_MZML_VENDOR_PRESENT_BUT_UNSUPPORTED_POLICY_PENDING")
        self.assertEqual(r["d.zip"]["reason"], "SKIP_ARCHIVE_CONTENTS_NOT_INSPECTED")
        self.assertEqual(r["e.mzXML"]["reason"], "SKIP_MZXML_READER_NOT_QUALIFIED")
        self.assertEqual(r["f.d.zip"]["reason"], "SKIP_BRUKER_D_NO_QUALIFIED_READER")
        self.assertEqual(r["g.raw.gz"]["reason"], "SKIP_COMPRESSED_RAW_NOT_QUALIFIED")
        self.assertEqual(r["search.mzid.gz"]["reason"], "SKIP_NOT_A_SPECTRUM_FILE")
        self.assertEqual(r["peaks.mgf"]["reason"], "SKIP_PEAK_LIST_NOT_A_RUN_REPRESENTATION")
        self.assertEqual(r["bad name.raw"]["reason"], "SKIP_UNSAFE_FILE_NAME")
        self.assertEqual(r["zero.raw"]["reason"], "SKIP_SIZE_UNKNOWN_OR_ZERO")
        self.assertEqual(summary["queue_files"], 2)
        self.assertFalse(summary["complete_study_asserted"])
        with (self.root / "out/queue.tsv").open() as handle:
            queue = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual({q["file_name"] for q in queue}, {"a.raw", "b.mzML.gz"})

    def test_count_mismatch_blocks_queueing_and_pages_are_followed(self):
        many = [record(f"r{i:03d}.raw") for i in range(INV.PAGE + 5)]
        api = FakeApi({"PXD000002": many, "PXD000003": [record("x.raw")]},
                      counts={"PXD000003": 2})
        summary = INV.build(["PXD000002", "PXD000003"], api, self.root / "out")
        r = rows(self.root / "out")
        self.assertEqual(summary["projects_verified"], 1)
        self.assertEqual(r["x.raw"]["reason"], "SKIP_PROJECT_LISTING_NOT_VERIFIED")
        self.assertEqual(sum(1 for k in r if k.startswith("r")), INV.PAGE + 5)
        self.assertEqual(summary["queue_files"], INV.PAGE + 5)

    def test_resume_from_run_index_and_cache(self):
        index = self.root / "index.jsonl"
        index.write_text("\n".join(json.dumps(x) for x in [
            {"accession": "PXD000004", "source_file": "done.raw", "role": "canonical"},
            {"accession": "PXD000004", "source_file": "bad.raw", "role": "failed"},
            {"accession": "PXD000004", "source_file": "fixed.raw", "role": "failed"},
            {"accession": "PXD000004", "source_file": "fixed.raw", "role": "canonical"},
        ]) + "\n")
        api = FakeApi({"PXD000004": [record("done.raw"), record("bad.raw"), record("fixed.raw"),
                                     record("new.raw")]})
        cache = self.root / "cache"
        summary = INV.build(["PXD000004"], api, self.root / "out", cache=cache, processed_index=index)
        r = rows(self.root / "out")
        self.assertEqual(r["done.raw"]["decision"], "DONE")
        self.assertEqual(r["fixed.raw"]["decision"], "DONE")
        self.assertEqual(r["bad.raw"]["decision"], "FAILED_PREVIOUSLY")
        self.assertEqual(r["new.raw"]["decision"], "QUEUE_NATIVE_RAW")
        self.assertEqual(summary["queue_files"], 1)
        calls = len(api.calls)
        INV.build(["PXD000004"], api, self.root / "out2", cache=cache)
        self.assertEqual(len(api.calls), calls)  # verified listing served from cache
        with self.assertRaisesRegex(ValueError, "Refusing existing"):
            INV.build(["PXD000004"], api, self.root / "out", cache=cache)

    def test_fetch_error_is_recorded_not_fatal(self):
        def broken(url):
            if "PXD000006" in url:
                raise OSError("boom")
            return FakeApi({"PXD000005": [record("ok.raw")]})(url)
        summary = INV.build(["PXD000005", "PXD000006"], broken, self.root / "out")
        self.assertEqual(summary["projects_with_error"], 1)
        self.assertEqual(summary["queue_files"], 1)

    def test_list_projects_paginates_and_filters(self):
        projects = {f"PXD{i:06d}": [] for i in range(INV.PAGE + 3)}
        api = FakeApi(projects)
        self.assertEqual(len(INV.list_projects(api)), INV.PAGE + 3)
        self.assertEqual(len(INV.list_projects(api, limit=5)), 5)


if __name__ == "__main__":
    unittest.main()
