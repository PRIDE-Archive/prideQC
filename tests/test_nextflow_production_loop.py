"""Unattended production loop: retry rules, guards, controller layout, stop conditions."""
from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import shutil
import tempfile
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
BIN = REPO / "nextflow/bin"


def load(name, file):
    spec = importlib.util.spec_from_file_location(name, BIN / file)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


L = load("prideqc_nf_loop", "production_loop.py")
PLAN = load("prideqc_nf_plan_loop", "plan_manifest.py")


def receipt(root, tag, acc, name, unit, **data):
    d = root / "results" / tag / acc / "_run_outputs/ready/units" / unit
    d.mkdir(parents=True)
    (d / "result.json").write_text(json.dumps(dict({"accession": acc, "source_member": name}, **data)))


class ProductionLoopTests(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.root = Path(tmp.name)
        (self.root / "results").mkdir()
        (self.root / "work/controllers").mkdir(parents=True)
        src = self.root / "src"
        (src / "nextflow/bin").mkdir(parents=True)
        shutil.copy(BIN / "queue_to_manifest.py", src / "nextflow/bin/queue_to_manifest.py")
        (src / ".source-sha").write_text("40b8b9a9aa43e21d90fc5d178d5a883c4999b05c\n")
        self.sif = self.root / "x.sif"
        self.sif.write_bytes(b"sif")
        rows = [{"accession": f"PXD00000{a}", "file_name": f"r{a}_{i}.raw", "file_format": "thermo_raw",
                 "file_size_bytes": str(10 ** 9 * (1 + i % 3)), "checksum": "", "memory_gb": "32",
                 "decision": "QUEUE_NATIVE_RAW"} for a in range(1, 4) for i in range(4)]
        rows.append({"accession": "PXD000009", "file_name": "big.raw", "file_format": "thermo_raw",
                     "file_size_bytes": str(9 * 10 ** 9), "checksum": "", "memory_gb": "128",
                     "decision": "QUEUE_NATIVE_RAW"})
        rows.append({"accession": "PXD000009", "file_name": "m.mzML", "file_format": "mzml",
                     "file_size_bytes": "10", "checksum": "", "memory_gb": "16", "decision": "QUEUE_MZML_NO_VENDOR"})
        self.queue = self.root / "queue.tsv"
        with self.queue.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
            w.writeheader()
            w.writerows(rows)

    def args(self, **kw):
        a = L.parser().parse_args([
            "--root", str(self.root), "--loop-id", "t1", "--source", str(self.root / "src"),
            "--sif", str(self.sif), "--sif-sha256", hashlib.sha256(b"sif").hexdigest(),
            "--queue", str(self.queue), "--batches", "2", "--first-batch", "3",
            "--max-files", "5", "--min-files", "1", "--min-free-gb", "0", "--min-hours-left", "0"])
        for k, v in kw.items():
            setattr(a, k, v)
        return a

    def fake_runner(self, success=True, rc=0, completed=True):
        def run(ctrl, env):
            prov = json.loads((ctrl / "screen-provenance.json").read_text())
            (ctrl / "launch-attempt.txt").write_text("x")
            with (ctrl / "file-manifest.tsv").open() as fh:
                for r in csv.DictReader(fh, delimiter="\t"):
                    receipt(self.root, prov["tag"], r["pxd_accession"], r["archive_file"], f"u{r['task_id']}",
                            status="success" if success else "failed",
                            failure_stage=None if success else "prideqc_analysis")
            with (ctrl / "file-manifest.tsv").open() as fh:
                n = sum(1 for _ in fh) - 1
            if completed:
                (ctrl / "nextflow.log").write_text(
                    f"Workflow completed > WorkflowStats[succeededCount={n}; failedCount=0; ...]")
            (ctrl / "controller.exitcode").write_text(f"{rc}\n")
            return rc
        return run

    def test_classify_retry_rules(self):
        rows = [{"a": "P1", "f": "ok", "s": "success"}, {"a": "P1", "f": "ok", "s": "failed", "st": "archive_preparation"},
                {"a": "P1", "f": "dl", "s": "failed", "st": "archive_preparation"},
                {"a": "P1", "f": "dl3", "s": "failed", "st": "archive_preparation"},
                {"a": "P1", "f": "dl3", "s": "failed", "st": "archive_preparation"},
                {"a": "P1", "f": "dl3", "s": "failed", "st": "archive_preparation"},
                {"a": "P1", "f": "oom", "s": "failed", "st": "prideqc_analysis", "x": -9, "e": "PrideqcNonzeroExit"},
                {"a": "P1", "f": "oom2", "s": "failed", "st": "prideqc_analysis", "x": -9, "e": "PrideqcNonzeroExit"},
                {"a": "P1", "f": "oom2", "s": "failed", "st": "prideqc_analysis", "x": -9, "e": "PrideqcNonzeroExit"},
                {"a": "P1", "f": "newoom", "s": "failed", "st": "prideqc_analysis", "x": -9,
                 "e": "AnalysisKilledResourceLimit"},
                {"a": "P1", "f": "reader", "s": "failed", "st": "prideqc_analysis", "e": "ThermoRawOpenError"}]
        skip, retry = L.classify(rows)
        self.assertEqual(retry, {("P1", "dl"), ("P1", "oom")})
        self.assertEqual(skip, {("P1", f) for f in ("ok", "dl3", "oom2", "newoom", "reader")})

    def test_dry_run_cuts_without_controller(self):
        busy = self.root / "work/controllers/live"
        busy.mkdir()
        (busy / "launch-attempt.txt").write_text("x")                      # dry run may run beside a live batch
        (busy / "file-manifest.tsv").write_text("task_id\tpxd_accession\tarchive_file\n1\tPXD000002\tr2_0.raw\n"
                                              "2\tPXD000003\tr3_1.raw\n")
        receipt(self.root, "old", "PXD000001", "r1_0.raw", "u1", status="success")
        receipt(self.root, "old", "PXD000002", "R2_0.RAW", "u2", status="failed", failure_stage="archive_preparation")
        result = L.Loop(self.args(dry_run=True)).run()
        self.assertEqual(result, "DRY_RUN")
        self.assertEqual([p.name for p in (self.root / "work/controllers").iterdir()], ["live"])
        cut = [json.loads(line) for line in (self.root / "work/loops/t1/loop.jsonl").read_text().splitlines()
               if json.loads(line)["event"] == "CUT"][0]
        self.assertEqual((cut["retry_keys"], cut["retry_in_batch"], cut["in_flight"]), (0, 0, 2))
        self.assertEqual(cut["files"], 5)
        with self.assertRaises(ValueError):   # no new controller: the dry run left nothing to collide with
            L.write_controller(self.root, "live", [], source=self.root / "src", sif=self.sif, sif_sha256="",
                               rules={}, nxf_opts="", max_forks=1, loop_id="t1")

    def test_two_batches_disjoint_and_controller_layout(self):
        result = L.Loop(self.args(), runner=self.fake_runner(), sleep=lambda s: None).run()
        self.assertEqual(result, "DONE")
        ctrls = sorted((self.root / "work/controllers").iterdir())
        self.assertEqual([c.name[:22] for c in ctrls], ["prideqc-prod-raw-b003-", "prideqc-prod-raw-b004-"])
        names = []
        for c in ctrls:
            with (c / "file-manifest.tsv").open() as fh:
                rows = list(csv.DictReader(fh, delimiter="\t"))
            names += [(r["pxd_accession"], r["archive_file"]) for r in rows]
            self.assertTrue(all(int(r["file_size_bytes"]) <= 6e9 and r["archive_file"].endswith(".raw") for r in rows))
            params = json.loads((c / "params.json").read_text())
            self.assertFalse(params["allow_complete_study"])
            runner = (c / "run-nextflow.sh").read_text()
            self.assertIn("export NXF_OPTS='-Xms1g -Xmx4g'", runner)
            self.assertIn("maxForks = 120", (c / "resources.config").read_text())
            prov = json.loads((c / "screen-provenance.json").read_text())
            self.assertEqual(prov["controller_sha256"]["run-nextflow.sh"],
                             hashlib.sha256(runner.encode()).hexdigest())
            plan = PLAN.prepare(c / "file-manifest.tsv", c / "report-manifest.tsv", self.root / f"p-{c.name}", set())
            self.assertEqual(plan["tasks"], len(rows))                          # pinned planner accepts it
        self.assertEqual(len(names), len(set(names)))                           # no file in two batches
        self.assertEqual(len(names), 10)

    def test_stops_on_low_success_and_never_relaunches_existing(self):
        result = L.Loop(self.args(), runner=self.fake_runner(success=False), sleep=lambda s: None).run()
        self.assertEqual(result, "STOP_LOW_SUCCESS")
        self.assertEqual(len(list((self.root / "work/controllers").iterdir())), 1)
        with self.assertRaises(ValueError):                                    # same tag again: refused
            L.Loop(self.args(batches=1), runner=self.fake_runner(), sleep=lambda s: None).run()

    def test_completed_workflow_with_controller_kill_still_passes(self):
        result = L.Loop(self.args(batches=1), runner=self.fake_runner(rc=137), sleep=lambda s: None).run()
        self.assertEqual(result, "DONE")
        result = L.Loop(self.args(batches=1, first_batch=9, loop_id="t2"),
                        runner=self.fake_runner(rc=1, completed=False), sleep=lambda s: None).run()
        self.assertEqual(result, "STOP_CONTROLLER_FAILED")

    def test_excluded_accession_is_never_cut(self):
        result = L.Loop(self.args(batches=1, dry_run=True, exclude_accession=["pxd000001"])).run()
        self.assertEqual(result, "DRY_RUN")
        cut = [json.loads(line) for line in (self.root / "work/loops/t1/loop.jsonl").read_text().splitlines()
               if json.loads(line)["event"] == "CUT"][0]
        self.assertEqual(cut["files"], 5)
        L.Loop(self.args(batches=1, loop_id="t3", exclude_accession=["PXD000001", "PXD000002"]),
               runner=self.fake_runner(), sleep=lambda s: None).run()
        ctrl = next((self.root / "work/controllers").iterdir())
        with (ctrl / "file-manifest.tsv").open() as fh:
            accs = {r["pxd_accession"] for r in csv.DictReader(fh, delimiter="\t")}
        self.assertEqual(accs, {"PXD000003"})
        prov = json.loads((ctrl / "screen-provenance.json").read_text())
        self.assertEqual(prov["rules"]["exclude_accession"], ["PXD000001", "PXD000002"])

    def test_download_failures_do_not_count_as_analysis_failures_until_outage_limit(self):
        def runner(fail_first):
            def run(ctrl, env):
                prov = json.loads((ctrl / "screen-provenance.json").read_text())
                (ctrl / "launch-attempt.txt").write_text("x")
                with (ctrl / "file-manifest.tsv").open() as fh:
                    tasks = list(csv.DictReader(fh, delimiter="\t"))
                for i, r in enumerate(tasks):
                    dl = i < fail_first
                    receipt(self.root, prov["tag"], r["pxd_accession"], r["archive_file"], f"u{r['task_id']}",
                            status="failed" if dl else "success",
                            failure_stage="archive_preparation" if dl else None,
                            **({"download_failure_kind": "timeout"} if dl and i == 0 else {}))
                (ctrl / "nextflow.log").write_text(
                    f"Workflow completed > WorkflowStats[succeededCount={len(tasks)}; failedCount=0; ...]")
                (ctrl / "controller.exitcode").write_text("0\n")
                return 0
            return run
        # 1 of 5 downloads fails (20%): below the 25% limit, analysis 4/4 = 100% -> PASS.
        result = L.Loop(self.args(batches=1), runner=runner(1), sleep=lambda s: None).run()
        self.assertEqual(result, "DONE")
        last = [json.loads(line) for line in (self.root / "work/loops/t1/loop.jsonl").read_text().splitlines()
                if json.loads(line)["event"] == "RESULT"][-1]
        self.assertEqual((last["success_rate"], last["analysis_success_rate"], last["download_failures"]),
                         (0.8, 1.0, 1))
        self.assertEqual(last["download_failure_kinds"], {"timeout": 1})
        # 2 of 5 (40%): above the limit -> STOP_DOWNLOAD_OUTAGE.
        result = L.Loop(self.args(batches=1, first_batch=8, loop_id="t5"), runner=runner(2),
                        sleep=lambda s: None).run()
        self.assertEqual(result, "STOP_DOWNLOAD_OUTAGE")
        last = [json.loads(line) for line in (self.root / "work/loops/t5/loop.jsonl").read_text().splitlines()
                if json.loads(line)["event"] == "RESULT"][-1]
        self.assertEqual(last["download_failure_kinds"], {"timeout": 1, "unrecorded": 1})

    def test_mzml_lane_only_accession_runs_beside_raw_loop(self):
        with self.queue.open("a", newline="") as fh:
            for i in range(3):
                fh.write(f"PXD082223\tm{i}.mzML\tmzml\t230000\t\t16\tQUEUE_MZML_NO_VENDOR\n")
        raw = self.root / "work/controllers/prideqc-prod-raw-b005-x"
        raw.mkdir()
        (raw / "launch-attempt.txt").write_text("x")
        # Without --concurrent-with the mzML loop waits for (and here stops on) the RAW controller.
        self.assertTrue(L.Loop(self.args(lane="QUEUE_MZML_NO_VENDOR", only_accession=["PXD082223"]),
                               runner=self.fake_runner(), sleep=lambda s: None).run().startswith("ACTIVE"))
        result = L.Loop(self.args(batches=1, loop_id="m1", lane="QUEUE_MZML_NO_VENDOR",
                                  only_accession=["PXD082223"], concurrent_with=["prideqc-prod-raw-"],
                                  tag_prefix="prideqc-prod-mzml"),
                        runner=self.fake_runner(), sleep=lambda s: None).run()
        self.assertEqual(result, "DONE")
        ctrl = next(c for c in (self.root / "work/controllers").iterdir() if c.name.startswith("prideqc-prod-mzml"))
        with (ctrl / "file-manifest.tsv").open() as fh:
            names = {r["archive_file"] for r in csv.DictReader(fh, delimiter="\t")}
        self.assertEqual(names, {"m0.mzML", "m1.mzML", "m2.mzML"})     # mzML only, PXD082223 only
        prov = json.loads((ctrl / "screen-provenance.json").read_text())
        self.assertEqual((prov["rules"]["lane"], prov["rules"]["only_accession"]),
                         ("QUEUE_MZML_NO_VENDOR", ["PXD082223"]))

    def test_killed_runner_still_gets_exitcode(self):
        def killed(ctrl, env):
            (ctrl / "launch-attempt.txt").write_text("x")
            return -9
        result = L.Loop(self.args(), runner=killed, sleep=lambda s: None).run()
        self.assertEqual(result, "STOP_CONTROLLER_FAILED")
        ctrl = next((self.root / "work/controllers").iterdir())
        self.assertEqual((ctrl / "controller.exitcode").read_text(), "-9\n")
        self.assertEqual(L.active_controllers(self.root), [])

    def test_guards(self):
        busy = self.root / "work/controllers/live"
        busy.mkdir()
        (busy / "launch-attempt.txt").write_text("x")
        self.assertTrue(L.Loop(self.args(), sleep=lambda s: None, runner=self.fake_runner()).run()
                        .startswith("ACTIVE_CONTROLLER"))
        (busy / "controller.exitcode").write_text("0\n")
        (self.root / "work/loops/t1/STOP").write_text("")
        self.assertEqual(L.Loop(self.args(), sleep=lambda s: None, runner=self.fake_runner()).run(), "STOP_FILE")
        (self.root / "work/loops/t1/STOP").unlink()
        self.assertEqual(L.Loop(self.args(min_hours_left=8.0), hours_left=lambda: 2.0,
                                runner=self.fake_runner()).run(), "SLURM_TIME_LOW")
        self.assertEqual(L.Loop(self.args(sif_sha256="0" * 64), runner=self.fake_runner()).run(), "SIF_HASH_MISMATCH")
        self.assertEqual(L.parse_slurm_left("6-23:59:10"), 6 * 24 + 23 + 59 / 60 + 10 / 3600)
        self.assertEqual(L.parse_slurm_left("45:00"), 0.75)

    def test_receipt_cache_only_for_finished_tags(self):
        receipt(self.root, "done", "PXD000001", "a.raw", "u1", status="success")
        receipt(self.root, "live", "PXD000001", "b.raw", "u1", status="success")
        live = self.root / "work/controllers/live"
        live.mkdir()
        rows = L.scan_results(self.root, self.root / "cache")
        self.assertEqual(len(rows), 2)
        self.assertEqual(sorted(p.name for p in (self.root / "cache").iterdir()), ["done.receipts.json"])


if __name__ == "__main__":
    unittest.main()
