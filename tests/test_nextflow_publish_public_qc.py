"""publish_public_qc.py: public folder per PRIDE project, stats, integrity, incremental writes."""
from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "nextflow/bin/publish_public_qc.py"
SPEC = importlib.util.spec_from_file_location("prideqc_nf_publish", SCRIPT)
PUB = importlib.util.module_from_spec(SPEC)
sys.modules["prideqc_nf_publish"] = PUB  # pool workers find the module functions by this name
SPEC.loader.exec_module(PUB)

MASS_SHIFT = {"accession": PUB.QS.MASS_SHIFTS, "name": "mass shifts", "value": [
    {"confidence": "high-support", "classification": "putative-ptm", "delta_mass_da": 15.9949,
     "unimod_candidates": [{"name": "Oxidation", "unimod_accession": "UNIMOD:35"}]}]}


def _sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


class PublishTests(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.root = Path(tmp.name) / "results"
        self.dest = Path(tmp.name) / "pride_qc"
        self.dest.mkdir()

    def _run(self, tag, pxd, name, *, task=1, status="success", created="2026-11-05T00:00:00", ms1=1000,
             duration=3600.0, tamper=False, extra=()):
        unit = self.root / tag / pxd / "_run_outputs/ready/units" / f"archive-{task:09d}-unit-00000-x"
        qc = unit / "qc"
        qc.mkdir(parents=True)
        receipt = {"accession": pxd, "task_id": task, "unit_id": unit.name, "source_archive": name,
                   "source_member": name, "status": status,
                   "failure_stage": None if status == "success" else "archive_preparation",
                   "download_failure_kind": None if status == "success" else "timeout"}
        if status == "success":
            metrics = [{"accession": "MS:4000059", "name": "number of MS1 spectra", "value": ms1},
                       {"accession": "MS:4000060", "name": "number of MS2 spectra", "value": 20000},
                       {"accession": "MS:4000067", "name": "run duration", "value": duration},
                       {"accession": "QCPRIDE:BIG", "name": "big evidence", "value": {"rows": [1, 2, 3]}},
                       *extra]
            mzqc = json.dumps({"mzQC": {"version": "1.0.0", "creationDate": created, "runQualities": [{
                "metadata": {"analysisSoftware": [{"name": "prideqc", "version": "0.2.0"}]},
                "qualityMetrics": metrics}]}}).encode()
            obo = b"format-version: 1.2\n"
            (qc / f"{name}.mzQC").write_bytes(mzqc + (b" " if tamper else b""))
            (qc / f"{name}.obo").write_bytes(obo)
            receipt.update(mzqc_file=f"{name}.mzQC", mzqc_sha256=_sha(mzqc), obo_file=f"{name}.obo",
                           obo_sha256=_sha(obo))
        (unit / "result.json").write_text(json.dumps(receipt))

    def _publish(self, **kw):
        return PUB.publish([self.root], self.dest, workers=2, **kw)

    def test_layout_project_json_and_stats(self):
        self._run("prideqc-prod-raw-b001-x", "PXD000001", "a.raw", task=1, ms1=1000, extra=[MASS_SHIFT])
        self._run("prideqc-prod-raw-b001-x", "PXD000001", "b.raw", task=2, ms1=3000)
        self._run("prideqc-prod-raw-b001-x", "PXD000001", "c.raw", task=3, status="failed")
        self._run("prideqc-prod-raw-b002-x", "PXD000002", "x.raw", task=1)
        self._run("prideqc-prod-raw-b002-x", "PXD000003", "y.raw", task=1, status="failed")
        manifest = self._publish()
        self.assertEqual((manifest["projects"], manifest["runs"], manifest["projectsWithoutQc"]), (2, 3, 1))
        self.assertEqual(manifest["failedFilesWithoutQc"], 1)
        p1 = self.dest / "PXD000001"
        self.assertEqual(sorted(f.name for f in p1.iterdir()),
                         ["a.raw.mzQC", "a.raw.obo", "b.raw.mzQC", "b.raw.obo", "project.json"])
        self.assertFalse((self.dest / "PXD000003").exists(), "a project without QC has no folder")
        doc = json.loads((p1 / "project.json").read_text())
        self.assertEqual(doc["summary"]["qcRuns"], 2)
        self.assertEqual(doc["summary"]["medianMs1Spectra"], 2000)
        self.assertEqual(doc["summary"]["medianRunDurationMin"], 60)
        self.assertEqual([r["sourceFile"] for r in doc["runs"]], ["a.raw", "b.raw"])
        self.assertEqual(doc["failedRuns"], [{"sourceFile": "c.raw", "failure_stage": "archive_preparation",
                                              "download_failure_kind": "timeout"}])
        run = doc["runs"][0]
        self.assertEqual(run["mzqcSha256"], _sha((p1 / "a.raw.mzQC").read_bytes()))
        self.assertEqual([m["accession"] for m in run["metrics"]], ["MS:4000059", "MS:4000060", "MS:4000067"],
                         "objects without a tolerance or cv_term are not in the short form")
        self.assertEqual(run["putativeModifications"],
                         [{"name": "Oxidation", "unimod": "UNIMOD:35", "deltaMassDa": 15.9949}])
        stats = json.loads((self.dest / "stats/archive-stats.json").read_text())
        self.assertEqual((stats["runs"], stats["projects"]), (3, 2))
        self.assertEqual(stats["putativeModifications"]["items"][0]["key"], "Oxidation")
        lines = (self.dest / "stats/projects.ndjson").read_text().splitlines()
        self.assertEqual([json.loads(x)["projectAccession"] for x in lines], ["PXD000001", "PXD000002"])
        self.assertTrue((self.dest / "README.md").is_file())
        self.assertEqual(json.loads((self.dest / "manifest.json").read_text())["runs"], 3)
        self.assertNotIn(str(self.root), (self.dest / "manifest.json").read_text(), "no internal path in public files")
        self.assertEqual([f.name for f in self.dest.rglob(".*tmp*")], [], "no temporary files left")

    def test_stats_equal_the_runs_ndjson_path(self):
        self._run("prideqc-prod-raw-b001-x", "PXD000001", "a.raw", task=1, ms1=1000, extra=[MASS_SHIFT])
        self._run("prideqc-prod-raw-b001-x", "PXD000002", "b.raw", task=1, ms1=7000, duration=7200.0)
        self._publish()
        acc = PUB.QS.Accumulator()
        for pxd in ("PXD000001", "PXD000002"):
            for f in sorted((self.dest / pxd).glob("*.mzQC")):
                q = json.loads(f.read_text())["mzQC"]
                acc.add({"pxd_accession": pxd, "metrics": PUB._metrics(q)})
        expected = acc.finish()
        got = json.loads((self.dest / "stats/archive-stats.json").read_text())
        self.assertEqual({k: got[k] for k in expected}, json.loads(json.dumps(expected)))

    def test_bad_checksum_and_test_tags_are_not_published(self):
        self._run("prideqc-prod-raw-b001-x", "PXD000001", "a.raw", task=1, tamper=True)
        self._run("mass-error-v12-smoke", "PXD000002", "b.raw", task=1)
        manifest = self._publish()
        self.assertEqual(manifest["projects"], 0)
        self.assertEqual(manifest["roles"].get("rejected"), 1)
        self.assertFalse((self.dest / "PXD000001").exists())
        self.assertFalse((self.dest / "PXD000002").exists())

    def test_newest_success_is_canonical_and_retry_supersedes_failure(self):
        self._run("prideqc-prod-raw-b001-x", "PXD000001", "a.raw", task=1, status="failed")
        self._run("prideqc-prod-raw-b002-x", "PXD000001", "a.raw", task=1, created="2026-11-06T00:00:00", ms1=5)
        self._run("prideqc-prod-raw-b003-x", "PXD000001", "a.raw", task=1, created="2026-11-07T00:00:00", ms1=9)
        self._publish()
        doc = json.loads((self.dest / "PXD000001/project.json").read_text())
        self.assertEqual(len(doc["runs"]), 1)
        self.assertEqual(doc["runs"][0]["runTag"], "prideqc-prod-raw-b003-x")
        self.assertEqual(doc["failedRuns"], [])

    def test_repeated_failures_give_one_entry(self):
        self._run("prideqc-prod-raw-b001-x", "PXD000001", "a.raw", task=1)
        self._run("prideqc-prod-raw-b001-x", "PXD000001", "z.raw", task=2, status="failed")
        self._run("prideqc-prod-raw-b002-x", "PXD000001", "z.raw", task=1, status="failed")
        manifest = self._publish()
        doc = json.loads((self.dest / "PXD000001/project.json").read_text())
        self.assertEqual([r["sourceFile"] for r in doc["failedRuns"]], ["z.raw"])
        self.assertEqual(manifest["failedFilesWithoutQc"], 1)

    def test_second_run_writes_only_new_files(self):
        self._run("prideqc-prod-raw-b001-x", "PXD000001", "a.raw", task=1)
        first = self._publish()
        self.assertEqual(first["mzqcFilesWritten"], 1)
        mtime = (self.dest / "PXD000001/a.raw.mzQC").stat().st_mtime_ns
        self._run("prideqc-prod-raw-b002-x", "PXD000001", "b.raw", task=1)
        second = self._publish()
        self.assertEqual(second["mzqcFilesWritten"], 1)
        self.assertEqual((self.dest / "PXD000001/a.raw.mzQC").stat().st_mtime_ns, mtime)
        self.assertEqual(json.loads((self.dest / "PXD000001/project.json").read_text())["summary"]["qcRuns"], 2)

    def test_dry_run_writes_nothing(self):
        self._run("prideqc-prod-raw-b001-x", "PXD000001", "a.raw", task=1)
        manifest = self._publish(dry_run=True)
        self.assertEqual((manifest["projects"], manifest["dryRun"]), (1, True))
        self.assertEqual(list(self.dest.iterdir()), [])

    def test_missing_destination_is_refused(self):
        with self.assertRaises(ValueError):
            PUB.publish([self.root], self.dest / "absent", workers=1)

    def test_compact_value(self):
        self.assertEqual(PUB.compact_value(3), 3)
        self.assertEqual(PUB.compact_value([1, 2.5]), [1, 2.5])
        self.assertIsNone(PUB.compact_value(True))
        self.assertIsNone(PUB.compact_value(["a"]))
        self.assertIsNone(PUB.compact_value({"rows": []}))
        self.assertEqual(PUB.compact_value({"suggested_tolerance": 10, "unit": "ppm", "evidence": [1] * 99}),
                         {"suggested_tolerance": 10, "unit": "ppm"})


if __name__ == "__main__":
    unittest.main()
