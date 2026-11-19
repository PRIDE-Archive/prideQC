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

    def run_it(self):
        inv = self.dir / "inventory.tsv"
        with inv.open("w", newline="") as h:
            w = csv.DictWriter(h, fieldnames=INV_FIELDS, delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(self.rows)
        out, exc = self.dir / "complete.tsv", self.dir / "excluded.tsv"
        summary = CP.run(inv, [self.root], out, exc, workers=2)
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

    def test_listing_not_verified_blocks(self):
        self.inv("PXD000010", "x.raw", "SKIP", "SKIP_PROJECT_LISTING_NOT_VERIFIED")
        self.inv("PXD000010", "y.raw", "QUEUE_NATIVE_RAW", "VENDOR_RAW_PREFERRED")
        self.receipt("prideqc-prod-raw-b001-x", "PXD000010", "y.raw")
        _, done, excluded = self.run_it()
        self.assertEqual(done, {})
        self.assertEqual(excluded["PXD000010"]["reason"], "listing_not_verified")


if __name__ == "__main__":
    unittest.main()
