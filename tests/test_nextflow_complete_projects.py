"""complete_projects.py: projects with prideQC metrics for all of their raw data."""
from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "nextflow/bin/complete_projects.py"
SPEC = importlib.util.spec_from_file_location("prideqc_nf_complete", SCRIPT)
CP = importlib.util.module_from_spec(SPEC)
sys.modules["prideqc_nf_complete"] = CP  # pool workers find the module functions by this name
SPEC.loader.exec_module(CP)

INV_FIELDS = ["accession", "file_name", "file_category", "file_format", "file_size_bytes", "checksum",
              "run_stem", "decision", "reason", "memory_gb", "representation_evidence"]


def _sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


class CompleteProjectsTests(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.dir = Path(tmp.name)
        self.root = self.dir / "results"
        self.rows = []

    def inv(self, pxd, name, decision, reason, fmt="thermo_raw", size=2_000_000_000):
        stem = name.lower()
        for suffix in (".wiff.scan", ".wiff", ".mzml.gz", ".mzml", ".raw", ".d.zip", ".mgf", ".zip", ".mzxml"):
            if stem.endswith(suffix):
                stem = stem[: -len(suffix)]
                break
        self.rows.append(dict(accession=pxd, file_name=name, file_category="RAW", file_format=fmt,
                              file_size_bytes=size, checksum="", run_stem=stem, decision=decision, reason=reason,
                              memory_gb=16, representation_evidence=""))

    def receipt(self, tag, pxd, name, *, task=1, status="success", created="2026-11-05T00:00:00"):
        unit = self.root / tag / pxd / "_run_outputs/ready/units" / f"archive-{task:09d}-unit-00000-x"
        (unit / "qc").mkdir(parents=True)
        receipt = {"accession": pxd, "task_id": task, "source_member": name, "source_archive": name,
                   "status": status, "failure_stage": None if status == "success" else "prideqc_analysis"}
        if status == "success":
            mzqc = json.dumps({"mzQC": {"version": "1.0.0", "creationDate": created, "runQualities": [{
                "metadata": {"analysisSoftware": [{"name": "prideqc", "version": "0.2.0"}]},
                "qualityMetrics": []}]}}).encode()
            obo = b"format-version: 1.2\n"
            (unit / "qc" / f"{name}.mzQC").write_bytes(mzqc)
            (unit / "qc" / f"{name}.obo").write_bytes(obo)
            receipt.update(mzqc_file=f"{name}.mzQC", mzqc_sha256=_sha(mzqc), obo_file=f"{name}.obo",
                           obo_sha256=_sha(obo))
        (unit / "result.json").write_text(json.dumps(receipt))

    def run_it(self, archives=None):
        inv = self.dir / "inventory.tsv"
        with inv.open("w", newline="") as h:
            w = csv.DictWriter(h, fieldnames=INV_FIELDS, delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(self.rows)
        out, exc = self.dir / "complete.tsv", self.dir / "excluded.tsv"
        summary = CP.run(inv, [self.root], out, exc, workers=2, archives=archives)
        def read(p):
            return {r["pride_accession"]: r for r in csv.DictReader(p.open(), delimiter="\t")}

        return summary, read(out), read(exc)

    def test_complete_and_excluded_projects(self):
        tag = "prideqc-prod-raw-b001-x"
        # A: two Thermo raw files with QC, an mzML of the same run and a result file -> complete
        self.inv("PXD000001", "a1.raw", "QUEUE_NATIVE_RAW", "VENDOR_RAW_PREFERRED")
        self.inv("PXD000001", "a2.raw", "QUEUE_NATIVE_RAW", "VENDOR_RAW_PREFERRED")
        self.inv("PXD000001", "a1.mzML", "SKIP", "SKIP_MZML_VENDOR_RAW_PREFERRED", "mzml")
        self.inv("PXD000001", "a.mzid", "SKIP", "SKIP_NOT_A_SPECTRUM_FILE", "result")
        self.inv("PXD000001", "a.mgf", "SKIP", "SKIP_PEAK_LIST_NOT_A_RUN_REPRESENTATION", "mgf")
        self.receipt(tag, "PXD000001", "a1.raw", task=1)
        self.receipt(tag, "PXD000001", "a2.raw", task=2)
        # B: WIFF + mzML of the same run, mzML has QC -> complete (mzML only)
        self.inv("PXD000002", "b1.wiff", "SKIP", "SKIP_SCIEX_WIFF_NO_QUALIFIED_READER", "sciex_wiff")
        self.inv("PXD000002", "b1.wiff.scan", "SKIP", "SKIP_SCIEX_WIFF_SCAN_COMPANION", "sciex_wiff_scan")
        self.inv("PXD000002", "b1.mzML", "QUEUE_MZML_NO_VENDOR", "VENDOR_FILE_HAS_NO_SUPPORTED_READER", "mzml")
        self.receipt(tag, "PXD000002", "b1.mzML")
        # C: one raw failed -> excluded failed_run
        self.inv("PXD000003", "c1.raw", "QUEUE_NATIVE_RAW", "VENDOR_RAW_PREFERRED")
        self.inv("PXD000003", "c2.raw", "QUEUE_NATIVE_RAW", "VENDOR_RAW_PREFERRED")
        self.receipt(tag, "PXD000003", "c1.raw", task=1)
        self.receipt(tag, "PXD000003", "c2.raw", task=2, status="failed")
        # D: raw with QC + a Bruker run without mzML -> excluded no_reader_bruker
        self.inv("PXD000004", "d1.raw", "QUEUE_NATIVE_RAW", "VENDOR_RAW_PREFERRED")
        self.inv("PXD000004", "d2.d.zip", "SKIP", "SKIP_BRUKER_D_NO_QUALIFIED_READER", "bruker_d_archive")
        self.receipt(tag, "PXD000004", "d1.raw")
        # E: one raw not processed yet -> excluded not_processed_yet
        self.inv("PXD000005", "e1.raw", "QUEUE_NATIVE_RAW", "VENDOR_RAW_PREFERRED")
        self.inv("PXD000005", "e2.raw", "QUEUE_NATIVE_RAW", "VENDOR_RAW_PREFERRED")
        self.receipt(tag, "PXD000005", "e1.raw")
        # F: failure, then success in a later batch -> complete
        self.inv("PXD000006", "f1.raw", "QUEUE_NATIVE_RAW", "VENDOR_RAW_PREFERRED")
        self.receipt(tag, "PXD000006", "f1.raw", status="failed")
        self.receipt("prideqc-prod-raw-b002-x", "PXD000006", "f1.raw")
        # G: only result files -> no QC input
        self.inv("PXD000007", "g.mzid", "SKIP", "SKIP_NOT_A_SPECTRUM_FILE", "result")
        # H: raw + mzML-only run (no vendor file) -> complete, both kinds
        self.inv("PXD000008", "h1.raw", "QUEUE_NATIVE_RAW", "VENDOR_RAW_PREFERRED")
        self.inv("PXD000008", "h2.mzML", "QUEUE_MZML_NO_VENDOR", "NO_VENDOR_FILE_WITH_SAME_STEM", "mzml")
        self.receipt(tag, "PXD000008", "h1.raw", task=1)
        self.receipt(tag, "PXD000008", "h2.mzML", task=2)
        # I: QC from a test run tag only -> not counted
        self.inv("PXD000009", "i1.raw", "QUEUE_NATIVE_RAW", "VENDOR_RAW_PREFERRED")
        self.receipt("mass-error-v12-smoke", "PXD000009", "i1.raw")

        summary, done, excluded = self.run_it()
        self.assertEqual(sorted(done), ["PXD000001", "PXD000002", "PXD000006", "PXD000008"])
        self.assertEqual((done["PXD000001"]["runs_with_qc"], done["PXD000001"]["representation"]), ("2", "vendor_raw"))
        self.assertEqual(done["PXD000001"]["covered_other_files"], "1")
        self.assertEqual(done["PXD000002"]["representation"], "mzml")
        self.assertEqual(done["PXD000002"]["mzml_runs_with_unsupported_vendor"], "1")
        self.assertEqual(done["PXD000002"]["covered_other_files"], "2")
        self.assertEqual(done["PXD000008"]["representation"], "vendor_raw+mzml")
        self.assertEqual(done["PXD000001"]["source_gb"], "4.0")
        # PXD000007 has no raw data at all: it is in neither list
        self.assertEqual({k: v["reason"] for k, v in excluded.items()},
                         {"PXD000003": "failed_run", "PXD000004": "no_reader_bruker",
                          "PXD000005": "not_processed_yet", "PXD000009": "not_processed_yet"})
        self.assertEqual(summary["complete_projects"], 4)
        self.assertEqual(summary["complete_runs"], 6)
        self.assertEqual(summary["projects_with_qc_inputs"], 8)

    def test_archive_members_are_qc_inputs(self):
        tag = "prideqc-prod-archive-b001-x"
        z = "SKIP_ARCHIVE_CONTENTS_NOT_INSPECTED"
        # J: zip with two runs (both QC) + results zip -> complete, 2 runs from archives
        self.inv("PXD000020", "runs.zip", "SKIP", z, "zip_archive")
        self.inv("PXD000020", "results.zip", "SKIP", z, "zip_archive")
        # K: zip with one run QC and one Bruker member -> no reader
        self.inv("PXD000021", "k.zip", "SKIP", z, "zip_archive")
        # L: only a 7z (not inspected) -> no QC input yet
        self.inv("PXD000022", "l.7z", "SKIP", z, "7z_archive")
        # M: zip download failed (one archive receipt) -> failed_run
        self.inv("PXD000023", "m.zip", "SKIP", z, "zip_archive")
        for pxd, name, task in (("PXD000020", "a.raw", 1), ("PXD000020", "b.mzML", 2), ("PXD000021", "k1.raw", 3)):
            self.receipt(tag, pxd, name, task=task)
        self.receipt(tag, "PXD000023", "m.zip", task=4, status="failed")
        folder = self.dir / "arch"
        folder.mkdir()
        (folder / "archives.tsv").write_text(
            "accession\tfile_name\tdecision\n"
            "PXD000020\truns.zip\tQUEUE_ARCHIVE_RUNS\nPXD000020\tresults.zip\tSKIP_ARCHIVE_NO_RUNS\n"
            "PXD000021\tk.zip\tQUEUE_ARCHIVE_RUNS\nPXD000022\tl.7z\tSKIP_ARCHIVE_FORMAT_NOT_INSPECTED\n"
            "PXD000023\tm.zip\tQUEUE_ARCHIVE_RUNS\n")
        (folder / "members.tsv").write_text(
            "accession\tarchive_file\tmember\tmember_file\tmember_size_bytes\tmember_format\tkind\n"
            "PXD000020\truns.zip\tx/a.raw\ta.raw\t100\traw\trun\n"
            "PXD000020\truns.zip\tb.mzML\tb.mzML\t50\tmzml\trun\n"
            "PXD000021\tk.zip\tk1.raw\tk1.raw\t10\traw\trun\n"
            "PXD000021\tk.zip\tk2.d/\tk2.d\t\t\tunsupported\n"
            "PXD000023\tm.zip\tm1.raw\tm1.raw\t10\traw\trun\n")
        summary, done, excluded = self.run_it(archives=folder)
        self.assertEqual(sorted(done), ["PXD000020"])
        self.assertEqual((done["PXD000020"]["runs_with_qc"], done["PXD000020"]["runs_from_archives"],
                          done["PXD000020"]["representation"]), ("2", "2", "vendor_raw+mzml"))
        self.assertEqual({k: v["reason"] for k, v in excluded.items()},
                         {"PXD000021": "no_reader_archive_member", "PXD000022": "no_qc_input",
                          "PXD000023": "failed_run"})
        # without --archives every zip blocks
        _, done2, _ = self.run_it()
        self.assertEqual(done2, {})

    def test_listing_not_verified_blocks(self):
        self.inv("PXD000010", "x.raw", "SKIP", "SKIP_PROJECT_LISTING_NOT_VERIFIED")
        self.inv("PXD000010", "y.raw", "QUEUE_NATIVE_RAW", "VENDOR_RAW_PREFERRED")
        self.receipt("prideqc-prod-raw-b001-x", "PXD000010", "y.raw")
        _, done, excluded = self.run_it()
        self.assertEqual(done, {})
        self.assertEqual(excluded["PXD000010"]["reason"], "listing_not_verified")


if __name__ == "__main__":
    unittest.main()
