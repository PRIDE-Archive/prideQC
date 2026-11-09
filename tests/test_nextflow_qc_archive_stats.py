from __future__ import annotations

import importlib.util
import json
import tempfile
import unittest
import unittest.mock
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "nextflow/bin/qc_archive_stats.py"
SPEC = importlib.util.spec_from_file_location("prideqc_qc_archive_stats", SCRIPT)
QS = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(QS)


def _run(pxd: str, name: str, *, rt=3600.0, ms1=3000, ms2=20000, prec=6.0, frag=(15.0, "ppm"),
         acq="Data-dependent acquisition", charge=2.0, ptms=()) -> dict:
    frag_acc = QS.FRAGMENT_TOL_PPM if frag and frag[1] == "ppm" else QS.FRAGMENT_TOL_DA
    regime = "high-resolution" if frag and frag[1] == "ppm" else "low-resolution"
    metrics = [
        {"accession": QS.RUN_DURATION, "value": rt},
        {"accession": QS.MS1_SPECTRA, "value": ms1},
        {"accession": QS.MS2_SPECTRA, "value": ms2},
        {"accession": QS.MEDIAN_CHARGE, "value": charge},
        {"accession": QS.MS1_CYCLE, "value": 1.2},
        {"accession": QS.ACQUISITION, "value": {"cv_term": {"name": acq}} if acq else {"classification": "x"}},
    ]
    if prec is not None:
        metrics.append({"accession": QS.PRECURSOR_TOL, "value": {"suggested_tolerance": prec, "unit": "ppm"}})
    if frag:
        metrics.append({"accession": frag_acc, "value": {"suggested_tolerance": frag[0], "unit": frag[1],
                                                        "resolution_regime": regime}})
    if ptms is not None:
        metrics.append({"accession": QS.MASS_SHIFTS, "value": [
            {"confidence": conf, "classification": cls, "delta_mass_da": 15.995,
             "unimod_candidates": [{"name": nm, "unimod_accession": "UniMod:35"}]} for nm, cls, conf in ptms]})
    return {"pxd_accession": pxd, "source_file": name, "metrics": metrics}


OX = ("Oxidation", "putative-ptm", "high-support")


class QcArchiveStatsTests(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.dir = Path(tmp.name)

    def _stats(self, runs):
        path = self.dir / "runs.ndjson"
        path.write_text("".join(json.dumps(r) + "\n" for r in runs), encoding="utf-8")
        return QS.compute(path)

    def test_runs_and_projects_weighting(self):
        # one big project (10 DIA runs) and two small DDA projects
        runs = [_run("PXD000001", f"a{i}.raw", acq="Data-independent acquisition") for i in range(10)]
        runs += [_run("PXD000002", "b.raw"), _run("PXD000003", "c.raw")]
        s = self._stats(runs)
        self.assertEqual((s["runs"], s["projects"]), (12, 3))
        acq = {i["key"]: (i["runs"], i["projects"]) for i in s["categories"]["acquisition_method"]["items"]}
        self.assertEqual(acq["Data-independent acquisition"], (10, 1))
        self.assertEqual(acq["Data-dependent acquisition"], (2, 2))

    def test_histogram_bins_and_project_medians(self):
        runs = [_run("PXD000001", "a.raw", prec=2.5), _run("PXD000001", "b.raw", prec=4.5),
                _run("PXD000002", "c.raw", prec=99.0)]
        h = self._stats(runs)["histograms"]["precursor_tolerance_ppm"]
        self.assertEqual(h["runs"]["n"], 3)
        self.assertEqual(h["runs"]["counts"][2], 1)
        self.assertEqual(h["runs"]["counts"][4], 1)
        self.assertEqual(h["runs"]["counts"][-1], 1, "values above the last edge go to the open last bin")
        self.assertEqual(h["projects"]["n"], 2)
        self.assertEqual(h["projects"]["counts"][3], 1, "project median 3.5 ppm")
        self.assertEqual(sum(h["runs"]["counts"]), 3)

    def test_log_histograms_and_unit_conversion(self):
        s = self._stats([_run("PXD000001", "a.raw", ms2=10 ** 4.1, rt=6000.0)])
        ms2 = s["histograms"]["ms2_spectra"]
        self.assertTrue(ms2["log10"])
        self.assertEqual(ms2["edges"][ms2["runs"]["counts"].index(1)], 4.0)
        self.assertAlmostEqual(s["histograms"]["run_duration_min"]["runs"]["quantiles"]["median"], 100.0)

    def test_categories_resolution_and_charge(self):
        runs = [_run("PXD000001", "a.raw", frag=(0.5, "Da"), charge=3.0), _run("PXD000001", "b.raw", frag=None,
                                                                            charge=None, acq=None)]
        cats = self._stats(runs)["categories"]
        res = {i["key"]: i["runs"] for i in cats["ms2_resolution"]["items"]}
        self.assertEqual(res, {"Low resolution": 1, "Not determined": 1})
        ch = {i["key"]: i["runs"] for i in cats["median_precursor_charge"]["items"]}
        self.assertEqual(ch, {"3+": 1, "Not determined": 1})
        acq = {i["key"]: i["runs"] for i in cats["acquisition_method"]["items"]}
        self.assertEqual(acq["Not determined"], 1)

    def test_putative_modifications_filter_and_counts(self):
        runs = [
            _run("PXD000001", "a.raw", ptms=[OX, ("Oxidation", "putative-ptm", "high-support")]),
            _run("PXD000001", "b.raw", ptms=[OX]),
            _run("PXD000002", "c.raw", ptms=[OX, ("Dimer", "isotope-like", "high-support"),
                                             ("Weak", "putative-ptm", "moderate-support")]),
            _run("PXD000003", "d.raw", ptms=None),
        ]
        pm = self._stats(runs)["putativeModifications"]
        self.assertEqual([i["key"] for i in pm["items"]], ["Oxidation"],
                         "isotope-like and non-high-support families are left out")
        self.assertEqual((pm["items"][0]["runs"], pm["items"][0]["projects"]), (3, 2))
        self.assertEqual(pm["items"][0]["deltaMassDa"], 15.995)
        self.assertEqual((pm["runsWithEvidence"], pm["projectsWithEvidence"]), (3, 2))

    def test_unsorted_input_is_rejected(self):
        runs = [_run("PXD000001", "a.raw"), _run("PXD000002", "b.raw"), _run("PXD000001", "c.raw")]
        with self.assertRaises(ValueError):
            self._stats(runs)

    def test_dry_run_writes_output_and_no_network(self):
        path = self.dir / "runs.ndjson"
        path.write_text(json.dumps(_run("PXD000001", "a.raw")) + "\n", encoding="utf-8")
        out = self.dir / "stats.json"
        with unittest.mock.patch.object(QS, "post") as post:
            self.assertEqual(QS.main([str(path), "--output", str(out)]), 0)
            post.assert_not_called()
        doc = json.loads(out.read_text())
        self.assertTrue(doc["id"].startswith("prideqc-archive-stats-"))
        self.assertIs(doc["scientificApproval"], False)
        self.assertEqual(len(doc["sourceSha256"]), 64)

    def test_execute_requires_key(self):
        path = self.dir / "runs.ndjson"
        path.write_text(json.dumps(_run("PXD000001", "a.raw")) + "\n", encoding="utf-8")
        with unittest.mock.patch.dict("os.environ", {}, clear=True), self.assertRaises(SystemExit):
            QS.main([str(path), "--output", str(self.dir / "s.json"), "--execute"])


if __name__ == "__main__":
    unittest.main()
