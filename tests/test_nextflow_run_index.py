"""Run-level export/index over Nextflow results: integrity, dedup, no DB writes."""
from __future__ import annotations

import hashlib
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "nextflow/bin/export_run_index.py"
SPEC = importlib.util.spec_from_file_location("prideqc_nf_run_index", SCRIPT)
IDX = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(IDX)


def _sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


class RunIndexTests(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.root = Path(tmp.name) / "results"

    def _run(self, tag, pxd, name, *, task=1, status="success", created="2026-11-05T00:00:00",
             source_sha=None, ms1=10, tamper=False):
        acc = self.root / tag / pxd
        unit = acc / "_run_outputs/ready/units" / f"archive-{task:09d}-unit-00000-x"
        qc = unit / "qc"
        qc.mkdir(parents=True)
        receipt = {"schema_version": 2, "accession": pxd, "task_id": task, "unit_id": unit.name,
                   "source_archive": name, "source_member": name, "status": status,
                   "failure_stage": None if status == "success" else "prideqc_analysis",
                   "error_type": None if status == "success" else "PrideqcNonzeroExit"}
        if status == "success":
            mzqc = json.dumps({"mzQC": {"version": "1.0.0", "creationDate": created, "runQualities": [{
                "metadata": {"label": name, "analysisSoftware": [
                    {"name": "prideqc", "version": "0.2.0"},
                    {"name": "software", "value": "pyOpenMS", "version": "3.7"}]},
                "qualityMetrics": [{"accession": "MS:4000059", "name": "number of MS1 spectra", "value": ms1},
                                   {"accession": "MS:4000060", "name": "number of MS2 spectra", "value": 20}]}]}}).encode()
            obo = b"format-version: 1.2\n"
            (qc / f"{name}.mzQC").write_bytes(mzqc)
            (qc / f"{name}.obo").write_bytes(obo)
            (qc / f"{name}.summary.json").write_text("{}")
            receipt.update(mzqc_file=f"{name}.mzQC", mzqc_sha256=_sha(mzqc),
                           obo_file=f"{name}.obo", obo_sha256=_sha(obo))
            if tamper:
                (qc / f"{name}.mzQC").write_bytes(mzqc + b" ")
        (unit / "result.json").write_text(json.dumps(receipt))
        if source_sha:
            status_dir = acc / "_archive_status"
            status_dir.mkdir(exist_ok=True)
            (status_dir / f"{task}.prepare-status.json").write_text(json.dumps(
                {"source_evidence": {"sha256": source_sha, "size_bytes": 123}}))

    def _index(self, out):
        return [json.loads(line) for line in (out / "index.jsonl").read_text().splitlines()]

    def test_canonical_selection_layout_and_dedup(self):
        self._run("old-tag", "PXD000001", "a.raw", created="2026-10-01T00:00:00")
        self._run("new-tag", "PXD000001", "a.raw", created="2026-11-01T00:00:00", source_sha="s" * 64)
        self._run("new-tag", "PXD000001", "b.raw", task=2, source_sha="t" * 64)
        self._run("neg-tag", "PXD000002", "bad.raw", status="failed")
        out = self.root.parent / "export"
        manifest = IDX.export([self.root], out)

        self.assertEqual(manifest["canonical_runs"], 2)
        self.assertEqual(manifest["roles"], {"canonical": 2, "failed": 1, "superseded": 1})
        self.assertFalse(manifest["mongo_write_performed"])
        self.assertTrue((out / "PXD000001/a.raw.mzQC").is_file())
        self.assertTrue((out / "PXD000001/b.raw.obo").is_file())
        rows = self._index(out)
        canonical_a = [r for r in rows if r["source_file"] == "a.raw" and r["role"] == "canonical"]
        self.assertEqual(canonical_a[0]["run_tag"], "new-tag")  # source evidence wins
        superseded = [r for r in rows if r["role"] == "superseded"][0]
        self.assertEqual(superseded["superseded_by"], canonical_a[0]["mzqc_sha256"])
        self.assertTrue(all("_mzqc" not in r for r in rows))
        runs = [json.loads(line) for line in (out / "runs.ndjson").read_text().splitlines()]
        self.assertEqual({r["source_file"] for r in runs}, {"a.raw", "b.raw"})
        self.assertEqual(runs[0]["metrics"][0]["value"], 10)
        self.assertFalse(runs[0]["scientific_approval"])
        neg = json.loads((out / "PXD000002/accession.json").read_text())
        self.assertEqual(neg["runs"][0]["role"], "failed")
        self.assertFalse(neg["complete_study"])
        # Existing export is never overwritten.
        with self.assertRaisesRegex(ValueError, "Refusing existing"):
            IDX.export([self.root], out)

    def test_tampered_mzqc_is_rejected_not_exported(self):
        self._run("tag", "PXD000003", "c.raw", tamper=True)
        out = self.root.parent / "export"
        manifest = IDX.export([self.root], out)
        self.assertEqual(manifest["canonical_runs"], 0)
        self.assertEqual(self._index(out)[0]["reject_reason"], "qc_sha256_mismatch")
        self.assertFalse((out / "PXD000003").exists())

    def test_failure_superseded_by_later_success_and_byte_duplicates_flagged(self):
        self._run("t1", "PXD000004", "d.raw", status="failed")
        self._run("t2", "PXD000004", "d.raw", source_sha="u" * 64)
        self._run("t2", "PXD000005", "copy_of_d.raw", source_sha="u" * 64)
        out = self.root.parent / "export"
        manifest = IDX.export([self.root], out)
        rows = self._index(out)
        self.assertIn("superseded_failure", {r["role"] for r in rows})
        self.assertEqual(manifest["duplicate_source_bytes_groups"], 1)
        dup = [r for r in rows if r.get("source_file") == "copy_of_d.raw"][0]
        self.assertEqual(dup["duplicate_source_bytes"], ["PXD000004/d.raw"])
        self.assertEqual(manifest["canonical_runs"], 2)  # flagged, not silently merged


if __name__ == "__main__":
    unittest.main()
