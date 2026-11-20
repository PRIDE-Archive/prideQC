from __future__ import annotations

import io
import json
import sys
import tempfile
import unittest
import zipfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "nextflow" / "bin"))
import archive_members as AM  # noqa: E402
import inspect_archives as IA  # noqa: E402


def make_zip(entries: dict) -> bytes:
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
        for name, data in entries.items():
            if name.endswith("/"):
                zf.writestr(zipfile.ZipInfo(name), b"")
            else:
                zf.writestr(name, data)
    return buf.getvalue()


class FakeResponse:
    def __init__(self, data: bytes, start: int, total: int, status: int = 206):
        self.data, self.status = data, status
        self.headers = {"Content-Range": f"bytes {start}-{start + len(data) - 1}/{total}"}

    def read(self):
        return self.data

    def __enter__(self):
        return self

    def __exit__(self, *a):
        return False


def range_opener(blob: bytes, log: list):
    def opener(request, timeout=None):
        spec = request.headers["Range"].split("=", 1)[1]
        log.append(spec)
        if spec.startswith("-"):
            start = max(0, len(blob) - int(spec[1:]))
            end = len(blob) - 1
        else:
            a, b = spec.split("-")
            start, end = int(a), min(int(b), len(blob) - 1)
        return FakeResponse(blob[start:end + 1], start, len(blob))
    return opener


class PlanMembersTest(unittest.TestCase):
    def plan(self, entries):
        return AM.list_zip(io.BytesIO(make_zip(entries)))

    def test_runs_unsupported_and_ignored(self):
        plan = self.plan({"a/run1.raw": b"x" * 10, "run2.mzML.gz": b"y", "res/peptides.txt": b"t",
                          "b.d/": b"", "b.d/analysis.tdf": b"z", "b.d/analysis.tdf_bin": b"z",
                          "w.raw/": b"", "w.raw/_FUNC001.DAT": b"q", "c.wiff": b"w", "__MACOSX/._a": b""})
        self.assertEqual([r["member"] for r in plan["runs"]], ["a/run1.raw", "run2.mzML.gz"])
        self.assertEqual(plan["runs"][0]["size"], 10)
        self.assertEqual(plan["unsupported"], ["b.d/", "c.wiff", "w.raw/"])
        self.assertEqual(plan["ignored"], 2)

    def test_vendor_first_over_mzml(self):
        plan = self.plan({"x/run1.raw": b"r", "y/run1.mzML": b"m", "run2.mzML": b"m"})
        self.assertEqual([r["member"] for r in plan["runs"]], ["run2.mzML", "x/run1.raw"])
        self.assertEqual(plan["skipped_mzml_with_raw"], 1)

    def test_unsafe_and_duplicate_members_fail_closed(self):
        for entries in ({"../evil.raw": b"x"}, {"/abs.raw": b"x"}, {"a/r.raw": b"x", "b/R.RAW": b"y"},
                        {"bad name.raw": b"x"}):
            with self.subTest(entries=list(entries)):
                with self.assertRaises(AM.ArchiveMemberError):
                    self.plan(entries)

    def test_not_a_zip(self):
        with self.assertRaises(AM.ArchiveMemberError):
            AM.list_zip(io.BytesIO(b"not a zip at all" * 100))


class RemoteListingTest(unittest.TestCase):
    def test_one_suffix_request_lists_small_directory(self):
        blob = make_zip({f"run{i}.raw": b"x" * 5000 for i in range(20)})
        log: list = []
        plan = AM.list_zip("https://example.org/a.zip", opener=range_opener(blob, log))
        self.assertEqual(len(plan["runs"]), 20)
        self.assertEqual(log, ["-1048576"])

    def test_large_file_reads_directory_from_tail_only(self):
        big = make_zip({"run1.raw": bytes(range(256)) * 20000, "notes.txt": b"n"})
        log: list = []
        file = AM.HttpRangeFile("https://example.org/a.zip", opener=range_opener(big, log), tail=4096)
        self.assertEqual(file.size, len(big))
        with zipfile.ZipFile(io.BufferedReader(file, 1 << 16)) as zf:
            self.assertEqual([i.filename for i in zf.infolist()], ["run1.raw", "notes.txt"])
        self.assertEqual(log[0], "-4096")

    def test_no_range_support_fails_closed(self):
        def opener(request, timeout=None):
            return FakeResponse(b"x", 0, 1, status=200)
        with self.assertRaises(AM.ArchiveMemberError):
            AM.list_zip("https://example.org/a.zip", opener=opener)


class ExtractMemberTest(unittest.TestCase):
    def test_extract_one_member_with_hash(self):
        import hashlib
        with tempfile.TemporaryDirectory() as tmp:
            archive = Path(tmp) / "a.zip"
            archive.write_bytes(make_zip({"d/run1.raw": b"abc" * 1000, "run2.raw": b"z"}))
            out = AM.extract_member(archive, "d/run1.raw", Path(tmp) / "work")
            self.assertEqual(out["path"].name, "run1.raw")
            self.assertEqual(out["size"], 3000)
            self.assertEqual(out["sha256"], hashlib.sha256(b"abc" * 1000).hexdigest())
            self.assertEqual(sorted(p.name for p in (Path(tmp) / "work").iterdir()), ["run1.raw"])
            with self.assertRaises(AM.ArchiveMemberError):
                AM.extract_member(archive, "d/run1.raw", Path(tmp) / "work2", max_bytes=100)

    def test_extract_refuses_non_run_member(self):
        with tempfile.TemporaryDirectory() as tmp:
            archive = Path(tmp) / "a.zip"
            archive.write_bytes(make_zip({"notes.txt": b"n"}))
            with self.assertRaises(AM.ArchiveMemberError):
                AM.extract_member(archive, "notes.txt", Path(tmp) / "w")


class InspectArchivesTest(unittest.TestCase):
    def write_inventory(self, path: Path) -> None:
        rows = [
            ("PXD1", "runs.zip", "RAW", "1000", "SKIP_ARCHIVE_CONTENTS_NOT_INSPECTED"),
            ("PXD1", "results.zip", "SEARCH", "50", "SKIP_ARCHIVE_CONTENTS_NOT_INSPECTED"),
            ("PXD1", "bruker.zip", "RAW", "50", "SKIP_ARCHIVE_CONTENTS_NOT_INSPECTED"),
            ("PXD1", "other.7z", "RAW", "70", "SKIP_ARCHIVE_CONTENTS_NOT_INSPECTED"),
            ("PXD1", "broken.zip", "RAW", "70", "SKIP_ARCHIVE_CONTENTS_NOT_INSPECTED"),
            ("PXD1", "direct.raw", "RAW", "10", "QUEUE_NATIVE_RAW"),
            ("PXD2", "nourl.zip", "RAW", "9", "SKIP_ARCHIVE_CONTENTS_NOT_INSPECTED"),
        ]
        with path.open("w", encoding="utf-8") as fh:
            fh.write("accession\tfile_name\tfile_category\tfile_size_bytes\treason\n")
            for r in rows:
                fh.write("\t".join(r) + "\n")

    def test_build_outputs(self):
        zips = {"runs.zip": {"a/r1.raw": b"x" * 2000, "r2.raw": b"y" * 5, "b.d/x.tdf": b"t"},
                "results.zip": {"peptides.txt": b"p"}, "bruker.zip": {"s.d/analysis.tdf": b"t"}}

        def fetch(url):
            if "PXD1" in url and url.endswith("page=0"):
                return [{"fileName": n, "publicFileLocations": [
                    {"name": "FTP Protocol", "value": f"ftp://ftp.pride.ebi.ac.uk/pride/data/archive/x/PXD1/{n}"}]}
                    for n in ["runs.zip", "results.zip", "bruker.zip", "other.7z", "broken.zip"]]
            return []

        def lister(url):
            name = url.rsplit("/", 1)[1]
            self.assertTrue(url.startswith("https://ftp.pride.ebi.ac.uk/"))
            if name == "broken.zip":
                raise AM.ArchiveMemberError("Not a valid ZIP archive: x")
            return AM.list_zip(io.BytesIO(make_zip(zips[name])))

        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            self.write_inventory(tmp / "inventory.tsv")
            summary = IA.build(tmp / "inventory.tsv", tmp / "out", cache=tmp / "cache", workers=2,
                               fetch=fetch, lister=lister)
            self.assertEqual(summary["by_decision"], {"QUEUE_ARCHIVE_RUNS": 1, "SKIP_ARCHIVE_FORMAT_NOT_INSPECTED": 1,
                                                      "SKIP_ARCHIVE_NO_RUNS": 1, "SKIP_ARCHIVE_UNREADABLE": 2,
                                                      "SKIP_ARCHIVE_UNSUPPORTED_RUNS": 1})
            self.assertEqual(summary["queue_runs"], 2)
            import csv
            queue = list(csv.DictReader(io.StringIO((tmp / "out/queue.tsv").read_text()), delimiter="\t"))
            self.assertEqual([(q["file_name"], q["memory_gb"], q["file_format"]) for q in queue],
                             [("runs.zip", "16", "zip_archive")])
            members = list(csv.DictReader(io.StringIO((tmp / "out/members.tsv").read_text()), delimiter="\t"))
            self.assertEqual([(m["archive_file"], m["member"], m["member_file"], m["kind"]) for m in members],
                             [("bruker.zip", "s.d/", "s.d", "unsupported"),
                              ("runs.zip", "a/r1.raw", "r1.raw", "run"), ("runs.zip", "r2.raw", "r2.raw", "run"),
                              ("runs.zip", "b.d/", "b.d", "unsupported")])
            archives = {r["file_name"]: r for r in csv.DictReader(io.StringIO((tmp / "out/archives.tsv").read_text()),
                                                                    delimiter="\t")}
            self.assertEqual(archives["runs.zip"]["reason"], "RUNS_AND_UNSUPPORTED_MEMBERS")
            self.assertEqual(archives["nourl.zip"]["reason"], "no_public_url")
            self.assertEqual(archives["broken.zip"]["reason"], "not_a_zip")
            # the cache keeps inspected archives only; a rerun does not list them again
            calls: list = []
            IA.build(tmp / "inventory.tsv", tmp / "out2", cache=tmp / "cache", workers=1, fetch=fetch,
                     lister=lambda url: calls.append(url) or lister(url))
            self.assertEqual([c.rsplit("/", 1)[1] for c in calls], ["broken.zip"])
            self.assertEqual((tmp / "out2/members.tsv").read_text(), (tmp / "out/members.tsv").read_text())
            with self.assertRaises(ValueError):
                IA.build(tmp / "inventory.tsv", tmp / "out", fetch=fetch, lister=lister)
            self.assertEqual(json.loads((tmp / "out/summary.json").read_text())["queue_archives"], 1)

    def test_memory_tier(self):
        self.assertEqual([IA.memory_tier_gb(s) for s in (0, 2_000_000_000, 5_000_000_000, 9_000_000_000)],
                         [16, 32, 64, 128])


class ArchiveLaneTest(unittest.TestCase):
    """source_qc archive lane: one receipt per run member, RAW never leaves scratch."""

    def setUp(self):
        import importlib
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.root = Path(tmp.name)
        self.analyze_unit = importlib.import_module("analyze_unit")
        self.source_qc = importlib.import_module("source_qc")
        self.finalize = importlib.import_module("finalize_accession").finalize
        self.archive = self.root / "runs.zip"
        self.archive.write_bytes(make_zip({"set/a.raw": b"A" * 300, "set/b.raw": b"B" * 200,
                                           "c.mzML": b"C", "set/a.mzML": b"dup of a", "peptides.txt": b"t",
                                           "x.d/analysis.tdf": b"bruker"}))
        self.task = self.root / "task.json"
        self.task.write_text(json.dumps({"accession": "PXD000123", "task_id": 7, "archive_file": "runs.zip",
                                         "expected_runs": 1, "source_manifest": "m.tsv",
                                         "source_kind": "archive_members"}))

    def fake_qc(self, fail: str | None = None, kill: str | None = None, seen: list | None = None):
        from types import SimpleNamespace

        def run(cmd, **kwargs):
            source = Path(cmd[2])
            if seen is not None:
                seen.append((source.name, sorted(p.name for p in source.parent.iterdir()), kwargs.get("timeout")))
            if source.name == kill:
                return SimpleNamespace(returncode=-9, stdout="", stderr="")
            if source.name == fail:
                return SimpleNamespace(returncode=1, stdout="", stderr="reader error")
            out = Path(cmd[cmd.index("-o") + 1])
            out.mkdir()
            (out / "manifest.json").write_text(json.dumps({"files": [{"status": "success"}]}))
            (out / f"{source.name}.mzQC").write_text(json.dumps({"mzQC": {"runQualities": [{"qualityMetrics": [
                {"accession": "MS:4000059", "value": 3}, {"accession": "MS:4000060", "value": 4}]}]}}))
            (out / f"{source.name}.obo").write_text("format-version: 1.2\n")
            return SimpleNamespace(returncode=0, stdout="", stderr="")
        return run

    def run_lane(self, **kw):
        from unittest.mock import patch
        with patch.object(self.analyze_unit.subprocess, "run", side_effect=self.fake_qc(**kw)), \
                patch.object(self.source_qc.Path, "cwd", return_value=self.root):
            return self.source_qc.run(self.task, "ftp", self.root / "node-local-raw", self.root / "ready",
                                      source=self.archive, mass_error=True)

    def receipts(self):
        return sorted((json.loads(p.read_text()) for p in (self.root / "ready/units").glob("*/result.json")),
                      key=lambda r: r["member_index"])

    def test_one_receipt_per_member_one_member_on_scratch(self):
        import hashlib
        seen: list = []
        status = self.run_lane(seen=seen)
        # one member on scratch at a time; vendor .raw first, mzML of the same run skipped
        self.assertEqual([(name, files) for name, files, _ in seen],
                         [("c.mzML", ["c.mzML", "unit.json"]), ("a.raw", ["a.raw", "unit.json"]),
                          ("b.raw", ["b.raw", "unit.json"])])
        self.assertTrue(all(t == self.analyze_unit.MIN_ANALYSIS_TIMEOUT_S for _, _, t in seen))
        self.assertEqual((status["emitted_units"], status["members_success"]), (3, 3))
        self.assertEqual(status["member_plan"]["unsupported"], ["x.d/"])
        self.assertEqual(status["member_plan"]["mzml_skipped_for_raw"], 1)
        rows = self.receipts()
        self.assertEqual([(r["source_archive"], r["source_member"], r["status"]) for r in rows],
                         [("runs.zip", "c.mzML", "success"), ("runs.zip", "set/a.raw", "success"),
                          ("runs.zip", "set/b.raw", "success")])
        self.assertEqual(rows[1]["source_member_sha256"], hashlib.sha256(b"A" * 300).hexdigest())
        self.assertEqual(rows[1]["source_member_size_bytes"], 300)
        self.assertEqual(rows[1]["source_archive_sha256"], hashlib.sha256(self.archive.read_bytes()).hexdigest())
        self.assertEqual(rows[0]["requested_estimators"]["mass_error"], True)
        self.assertEqual({r["source_kind"] for r in rows}, {"archive_member"})
        self.assertFalse((self.root / "node-local-raw").exists())
        staged = [p.name.casefold() for p in (self.root / "ready").rglob("*")]
        self.assertFalse(any(n.endswith((".raw", ".mzml", ".zip")) for n in staged))
        self.assertEqual(len(list((self.root / "ready").rglob("*.mzQC"))), 3)
        # the accession finalizer counts members (units), and one archive task
        summary = self.finalize("PXD000123", 3, self.root / "ready/units", self.root / "final")
        self.assertEqual((summary["successful_runs"], summary["failed_runs"]), (3, 0))

    def test_failed_member_does_not_hide_the_others(self):
        status = self.run_lane(fail="a.raw")
        rows = self.receipts()
        self.assertEqual([r["status"] for r in rows], ["success", "failed", "success"])
        self.assertEqual(rows[1]["error_type"], "PrideqcNonzeroExit")
        self.assertEqual(status["members_success"], 2)
        self.assertEqual(self.source_qc.retry_exit_code(status, 1, 2), 0)

    def test_killed_member_asks_for_memory_retry(self):
        status = self.run_lane(kill="b.raw")
        self.assertEqual(status["members_killed"], 1)
        self.assertEqual(self.source_qc.retry_exit_code(status, 1, 2), 137)
        self.assertEqual(self.source_qc.retry_exit_code(status, 2, 2), 0)

    def test_worker_kill_keeps_done_members_and_fails_closed_for_the_rest(self):
        import signal
        out = self.root / "ready"
        units = out / "units"
        plan = {"base": {"schema_version": 2, "accession": "PXD000123", "task_id": 7, "source_archive": "runs.zip",
                         "expected_file_tasks": 1, "source_kind": "archive_member"},
                "members": [{"unit_id": f"u{i}", "member": f"m{i}.raw", "member_index": i} for i in range(3)]}
        (units / "u0").mkdir(parents=True)
        (units / "u0/result.json").write_text(json.dumps({"status": "success", "unit_id": "u0"}))
        (units / "u1").mkdir()
        (units / "u1/partial.log").write_text("x")
        (out / "member-plan.json").write_text(json.dumps(plan))
        status = self.source_qc.killed_receipt(self.task, self.root / "node-local-raw", out, signal.SIGKILL)
        got = {p.parent.name: json.loads(p.read_text()) for p in units.glob("*/result.json")}
        self.assertEqual(got["u0"]["status"], "success")
        self.assertEqual((got["u1"]["error_type"], got["u1"]["source_member"]), ("AnalysisKilledResourceLimit", "m1.raw"))
        self.assertEqual(got["u2"]["error_type"], "NotAttemptedAfterWorkerKill")
        self.assertFalse((units / "u1/partial.log").exists())
        self.assertEqual((status["emitted_units"], status["members_killed"]), (3, 1))

    def test_archive_without_runs_gives_one_failure_receipt(self):
        self.archive.write_bytes(make_zip({"peptides.txt": b"t"}))
        status = self.run_lane()
        rows = [json.loads(p.read_text()) for p in (self.root / "ready/units").glob("*/result.json")]
        self.assertEqual([(r["status"], r["failure_stage"], r["error_type"]) for r in rows],
                         [("failed", "archive_listing", "ArchiveMemberError")])
        self.assertEqual(status["emitted_units"], 1)
        self.assertFalse((self.root / "node-local-raw").exists())


class PlanManifestArchiveTest(unittest.TestCase):
    def test_archive_task_resources(self):
        import importlib
        pm = importlib.import_module("plan_manifest")
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            (tmp / "files.tsv").write_text(
                "task_id\tpxd_accession\tarchive_file\tsdrf_path\tneeds_file_map\tfile_size_bytes\tsource_kind"
                "\tlargest_member_bytes\tmember_bytes\n"
                "1\tPXD000123\truns.zip\t\t0\t9000000000\tarchive_members\t2500000000\t20000000000\n"
                "2\tPXD000123\tone.raw\t\t0\t500\t\t\t\n")
            (tmp / "reports.tsv").write_text("report_id\tpxd_accession\texpected_file_tasks\n1\tPXD000123\t2\n")
            pm.prepare(tmp / "files.tsv", tmp / "reports.tsv", tmp / "plan", set())
            zip_task = json.loads((tmp / "plan/tasks/PXD000123__000000001.json").read_text())
            raw_task = json.loads((tmp / "plan/tasks/PXD000123__000000002.json").read_text())
            self.assertEqual((zip_task["source_kind"], zip_task["memory_gb"], zip_task["time_hours"]),
                             ("archive_members", 32, 11))
            self.assertNotIn("time_hours", raw_task)
            self.assertEqual(pm.archive_time_hours(1), 4)
            self.assertEqual(pm.archive_time_hours(10**15), 72)
            bad = (tmp / "files.tsv").read_text().replace("2500000000", "")
            (tmp / "bad.tsv").write_text(bad)
            with self.assertRaises(ValueError):
                pm.prepare(tmp / "bad.tsv", tmp / "reports.tsv", tmp / "plan2", set())


if __name__ == "__main__":
    unittest.main()
