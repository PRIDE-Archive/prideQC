from __future__ import annotations

import importlib.util
import json
import tempfile
import unittest
import unittest.mock
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "nextflow/bin/qc_project_summary.py"
SPEC = importlib.util.spec_from_file_location("prideqc_qc_project_summary", SCRIPT)
PS = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(PS)
QS = PS.QS


def _run(pxd: str, name: str, *, rt=3600.0, ms1=3000, ms2=20000, acq="Data-dependent acquisition",
         regime="high-resolution") -> dict:
    metrics = [
        {"accession": QS.RUN_DURATION, "value": rt},
        {"accession": QS.MS1_SPECTRA, "value": ms1},
        {"accession": QS.MS2_SPECTRA, "value": ms2},
        {"accession": QS.ACQUISITION, "value": {"cv_term": {"name": acq}} if acq else {"classification": "x"}},
    ]
    if regime:
        metrics.append({"accession": QS.FRAGMENT_TOL_PPM, "value": {"suggested_tolerance": 10.0, "unit": "ppm",
                                                                   "resolution_regime": regime}})
    return {"pxd_accession": pxd, "source_file": name, "metrics": metrics}


class QcProjectSummaryTests(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.dir = Path(tmp.name)

    def _write(self, runs) -> Path:
        path = self.dir / "runs.ndjson"
        path.write_text("".join(json.dumps(r) + "\n" for r in runs), encoding="utf-8")
        return path

    def test_medians_buckets_and_sets(self):
        runs = [_run("PXD000001", "a.raw", rt=1200.0, ms1=800, ms2=5000),
                _run("PXD000001", "b.raw", rt=3000.0, ms1=1200, ms2=40000, acq="Data-independent acquisition"),
                _run("PXD000001", "c.raw", rt=4800.0, ms1=900, ms2=150000, acq=None, regime="low-resolution")]
        [s] = list(PS.summaries(self._write(runs)))
        self.assertEqual(s["qcRuns"], 3)
        self.assertEqual(s["medianRunDurationMin"], 50.0)
        self.assertEqual(s["runDurationBucket"], "30 to 60 min")
        self.assertEqual((s["medianMs1Spectra"], s["ms1SpectraBucket"]), (900, "under 1k"))
        self.assertEqual((s["medianMs2Spectra"], s["ms2SpectraBucket"]), (40000, "30k to 60k"))
        self.assertEqual(s["acquisitionMethods"], ["Data-dependent acquisition", "Data-independent acquisition"],
                         "'Not determined' is dropped when another run has a value")
        self.assertEqual(s["ms2Resolutions"], ["High resolution", "Low resolution"])

    def test_not_determined_only_when_nothing_known(self):
        [s] = list(PS.summaries(self._write([_run("PXD000002", "a.raw", acq=None, regime=None)])))
        self.assertEqual(s["acquisitionMethods"], ["Not determined"])
        self.assertEqual(s["ms2Resolutions"], ["Not determined"])

    def test_bucket_edges(self):
        self.assertEqual(PS.bucket("run_duration_min", 29.99), "under 30 min")
        self.assertEqual(PS.bucket("run_duration_min", 30.0), "30 to 60 min")
        self.assertEqual(PS.bucket("run_duration_min", 500.0), "180 min or more")
        self.assertEqual(PS.bucket("ms2_spectra", 120000), "120k or more")
        self.assertIsNone(PS.bucket("ms1_spectra", None))

    def test_labels_are_safe_filter_values(self):
        for buckets in PS.BUCKETS.values():
            for _, label in buckets:
                for bad in (",", "==", "*", "?"):
                    self.assertNotIn(bad, label)

    def test_projects_split_and_unsorted_rejected(self):
        out = list(PS.summaries(self._write([_run("PXD000001", "a.raw"), _run("PXD000002", "b.raw")])))
        self.assertEqual([s["id"] for s in out], ["PXD000001", "PXD000002"])
        with self.assertRaises(ValueError):
            list(PS.summaries(self._write([_run("PXD000001", "a.raw"), _run("PXD000002", "b.raw"),
                                           _run("PXD000001", "c.raw")])))

    def test_missing_metrics_give_no_bucket(self):
        run = {"pxd_accession": "PXD000003", "source_file": "x.raw", "metrics": []}
        [s] = list(PS.summaries(self._write([run])))
        self.assertIsNone(s["medianRunDurationMin"])
        self.assertIsNone(s["runDurationBucket"])

    def test_dry_run_writes_ndjson_without_network(self):
        path = self._write([_run("PXD000001", "a.raw"), _run("PXD000002", "b.raw")])
        out = self.dir / "projects.ndjson"
        with unittest.mock.patch.object(PS, "post") as post:
            self.assertEqual(PS.main([str(path), "--output", str(out)]), 0)
            post.assert_not_called()
        docs = [json.loads(line) for line in out.read_text(encoding="utf-8").splitlines()]
        self.assertEqual(len(docs), 2)
        self.assertFalse(docs[0]["scientificApproval"])
        self.assertEqual(len(docs[0]["sourceSha256"]), 64)

    def test_execute_batches_and_needs_key(self):
        path = self._write([_run(f"PXD{i:06d}", "a.raw") for i in range(1, 6)])
        out = self.dir / "p.ndjson"
        with self.assertRaises(SystemExit), unittest.mock.patch.dict("os.environ", {}, clear=True):
            PS.main([str(path), "--output", str(out), "--execute"])
        with unittest.mock.patch.object(PS, "post", side_effect=lambda u, h, k, b: b) as post, \
                unittest.mock.patch.dict("os.environ", {"PRIDE_MONGO_WS_API_KEY": "k"}):
            self.assertEqual(PS.main([str(path), "--output", str(out), "--execute", "--batch-size", "2"]), 0)
        self.assertEqual([len(c.args[3]) for c in post.call_args_list], [2, 2, 1])
        self.assertTrue(post.call_args_list[0].args[0].endswith("/qc-project/batch"))


if __name__ == "__main__":
    unittest.main()
