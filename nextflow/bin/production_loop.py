#!/usr/bin/env python3
"""Run N production RAW batches in sequence, unattended, inside one Slurm job.

One operator approval covers N batches. For each batch:
  1. guards: STOP file, no other active controller, free disk, Slurm time left;
  2. scan receipts under <root>/results (finished run tags are cached);
  3. cut a batch from the inventory queue with queue_to_manifest.select rules;
  4. write a new immutable controller (refuses an existing one);
  5. run Nextflow in the foreground (run-nextflow.sh writes controller.exitcode);
  6. check: the workflow completed, every task has a receipt, success rate >= --min-success.
The loop stops at the first failed check. It never deletes data. No Mongo or SDRF writes.
allow_complete_study is always false. Python 3.9 compatible (Codon host python3).

Retry rules (per (accession, file), from the receipts):
  - any success: done, never queued again;
  - download/preparation failures only: retried while there are fewer than 3 such receipts;
  - legacy OOM (exit -9 + PrideqcNonzeroExit, engine before 6431b84): retried once;
  - every other failure (reader, AnalysisKilledResourceLimit after the 2x retry): not retried.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import importlib.util
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
import time
from collections import defaultdict
from collections.abc import Callable
from pathlib import Path
from typing import Any

Key = tuple[str, str]
# Codon host python3 is 3.9: datetime.UTC (3.11+) does not exist there.
UTC = dt.timezone.utc  # noqa: UP017
MAX_DOWNLOAD_RECEIPTS = 3
COMPLETED = re.compile(r"Workflow completed > WorkflowStats\[succeededCount=(\d+); failedCount=(\d+)")
FILE_NAME = re.compile(r"[A-Za-z0-9_.+\-]+\.raw\Z", flags=re.IGNORECASE)
ACCESSION = re.compile(r"PXD\d{6,}\Z")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def now() -> str:
    return dt.datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ")


# ---------------------------------------------------------------- receipts
def scan_tag(tag_dir: Path) -> list[dict]:
    """One compact row per receipt of one run tag. Never raises on a bad receipt."""
    rows = []
    for receipt in sorted(tag_dir.glob("*/_run_outputs/ready/units/*/result.json")):
        try:
            if receipt.is_symlink() or receipt.stat().st_size > 2_000_000:
                continue
            d = json.loads(receipt.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            continue
        if not isinstance(d, dict) or not isinstance(d.get("accession"), str):
            continue
        name = Path(str(d.get("source_member") or d.get("source_archive") or "")).name
        if not name:
            continue
        rows.append({"a": d["accession"].upper(), "f": name.casefold(), "s": d.get("status"),
                     "st": d.get("failure_stage"), "x": d.get("exit_code"), "e": d.get("error_type")})
    return rows


def tag_finished(root: Path, tag: str) -> bool:
    ctrl = root / "work/controllers" / tag
    return (ctrl / "controller.exitcode").is_file() or not ctrl.is_dir()


def scan_results(root: Path, cache: Path) -> list[dict]:
    """Receipts of every run tag. A finished tag is scanned once and cached in `cache`."""
    cache.mkdir(parents=True, exist_ok=True)
    rows: list[dict] = []
    for tag_dir in sorted(p for p in (root / "results").iterdir() if p.is_dir() and not p.is_symlink()):
        cached = cache / f"{tag_dir.name}.receipts.json"
        if cached.is_file():
            rows.extend(json.loads(cached.read_text(encoding="utf-8")))
            continue
        found = scan_tag(tag_dir)
        if tag_finished(root, tag_dir.name):
            tmp = cached.with_suffix(".partial")
            tmp.write_text(json.dumps(found), encoding="utf-8")
            tmp.rename(cached)
        rows.extend(found)
    return rows


def is_legacy_oom(r: dict) -> bool:
    return r.get("st") == "prideqc_analysis" and r.get("x") == -9 and r.get("e") == "PrideqcNonzeroExit"


def classify(rows: list[dict]) -> tuple[set[Key], set[Key]]:
    """Return (skip, retry_first). Keys are (ACCESSION, file casefold)."""
    by: dict[Key, list[dict]] = defaultdict(list)
    for r in rows:
        by[(r["a"], r["f"])].append(r)
    skip: set[Key] = set()
    retry: set[Key] = set()
    for key, items in by.items():
        if any(r.get("s") == "success" for r in items):
            skip.add(key)
        elif all(r.get("st") == "archive_preparation" for r in items) and len(items) < MAX_DOWNLOAD_RECEIPTS:
            retry.add(key)
        elif len(items) == 1 and is_legacy_oom(items[0]):
            retry.add(key)
        else:
            skip.add(key)
    return skip, retry


# ---------------------------------------------------------------- guards
def active_controllers(root: Path, exclude: str = "") -> list[str]:
    """Controllers with a launch marker and no exit code (a Nextflow run may be live)."""
    ctrls = root / "work/controllers"
    if not ctrls.is_dir():
        return []
    return sorted(c.name for c in ctrls.iterdir()
                  if c.is_dir() and c.name != exclude and (c / "launch-attempt.txt").is_file()
                  and not (c / "controller.exitcode").is_file())


def in_flight_keys(root: Path) -> set[Key]:
    """Files listed in the manifest of an active controller (receipts may not exist yet)."""
    keys: set[Key] = set()
    for tag in active_controllers(root):
        manifest = root / "work/controllers" / tag / "file-manifest.tsv"
        if manifest.is_file():
            with manifest.open(encoding="utf-8", newline="") as fh:
                for r in csv.DictReader(fh, delimiter="\t"):
                    keys.add((r["pxd_accession"].upper(), r["archive_file"].casefold()))
    return keys


def slurm_hours_left(job_id: str | None) -> float | None:
    if not job_id:
        return None
    try:
        out = subprocess.run(["squeue", "-h", "-j", job_id, "-o", "%L"], capture_output=True,
                             text=True, timeout=60).stdout.strip()
    except (OSError, subprocess.SubprocessError):
        return None
    return parse_slurm_left(out)


def parse_slurm_left(text: str) -> float | None:
    m = re.fullmatch(r"(?:(\d+)-)?(?:(\d+):)?(\d+):(\d+)", text.strip())
    if not m:
        return None
    d, h, mi, s = (int(g) if g else 0 for g in m.groups())
    return d * 24 + h + mi / 60 + s / 3600


def free_gb(path: Path) -> float:
    return shutil.disk_usage(path).free / 1e9


# ---------------------------------------------------------------- controller
def load_queue_module(source: Path) -> Any:
    spec = importlib.util.spec_from_file_location("prideqc_loop_q", source / "nextflow/bin/queue_to_manifest.py")
    if spec is None or spec.loader is None:
        raise ValueError("QUEUE_MODULE_NOT_LOADABLE")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def memory_tier_gb(size: int) -> int:
    # Mirror of nextflow/bin/plan_manifest.py; the pipeline's own PRIDEQC_MANIFEST step re-derives it.
    for bound, gb in ((1_000_000_000, 16), (3_000_000_000, 32), (6_000_000_000, 64)):
        if size < bound:
            return gb
    return 128


def write_tsv(path: Path, fields: list[str], rows: list[dict]) -> None:
    with path.open("x", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def resources_config(max_forks: int) -> str:
    return ("// Production loop controller resource override.\n"
            "process {\n    withName: PRIDEQC_SOURCE_QC {\n"
            f"        maxForks = {max_forks}\n    }}\n}}\n"
            f"executor {{\n    queueSize = {max_forks + 8}\n    submitRateLimit = '5/1sec'\n}}\n")


def write_controller(root: Path, tag: str, picked: list[dict], *, source: Path, sif: Path,
                     sif_sha256: str, rules: dict, nxf_opts: str, max_forks: int, loop_id: str) -> Path:
    ctrl_root = root / "work/controllers"
    ctrl, work, out = ctrl_root / tag, root / "work" / tag, root / "results" / tag
    if any(p.exists() or p.is_symlink() for p in (ctrl, work, out)):
        raise ValueError(f"EXISTING_CONTROLLER_OUTPUT_OR_WORK_{tag}")
    for r in picked:
        if not ACCESSION.fullmatch(r["accession"]) or not FILE_NAME.fullmatch(r["file_name"]):
            raise ValueError(f"UNSAFE_ROW_{r['accession']}_{r['file_name']}")
    ctrl.mkdir(mode=0o750)
    tasks = [{"task_id": str(i), "pxd_accession": r["accession"], "archive_file": r["file_name"],
              "sdrf_path": "", "needs_file_map": "0", "file_size_bytes": str(r["file_size_bytes"])}
             for i, r in enumerate(picked, 1)]
    per: dict[str, int] = defaultdict(int)
    for r in picked:
        per[r["accession"]] += 1
    reports = [{"report_id": str(i), "pxd_accession": a, "expected_file_tasks": str(n)}
               for i, (a, n) in enumerate(sorted(per.items()), 1)]
    write_tsv(ctrl / "file-manifest.tsv", ["task_id", "pxd_accession", "archive_file", "sdrf_path",
                                           "needs_file_map", "file_size_bytes"], tasks)
    write_tsv(ctrl / "report-manifest.tsv", ["report_id", "pxd_accession", "expected_file_tasks"], reports)
    (ctrl / "resources.config").write_text(resources_config(max_forks), encoding="utf-8")
    params = {"input": str(ctrl / "file-manifest.tsv"), "report_manifest": str(ctrl / "report-manifest.tsv"),
              "outdir": str(out), "sif": str(sif), "protocol": "fire", "checksum_check": True,
              "allow_complete_study": False, "estimate_mass_error": True, "estimate_mass_shifts": True,
              "tmt_reporter_evidence": True, "tmt_profile_reporter_evidence": True}
    (ctrl / "params.json").write_text(json.dumps(params, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    q = shlex.quote
    invoke = (f"nextflow -log {q(str(ctrl / 'nextflow.log'))} run ./main.nf -profile codon_slurm"
              f" -c {q(str(ctrl / 'resources.config'))} -params-file {q(str(ctrl / 'params.json'))}"
              f" -work-dir {q(str(work))} -ansi-log false")
    runner = ("#!/usr/bin/env bash\numask 027\n"
              f"export NXF_OPTS={q(nxf_opts)}\n"
              f"cd {q(str(source / 'nextflow'))} || {{ printf '93\\n' > {q(str(ctrl / 'controller.exitcode'))}; exit 93; }}\n"
              f"{invoke} > {q(str(ctrl / 'controller.log'))} 2>&1\n"
              "status=$?\n"
              f"printf '%s\\n' \"$status\" > {q(str(ctrl / 'controller.exitcode'))}\n"
              'exit "$status"\n')
    (ctrl / "run-nextflow.sh").write_text(runner, encoding="utf-8")
    (ctrl / "run-nextflow.sh").chmod(0o750)
    tiers: dict[str, int] = defaultdict(int)
    for r in picked:
        tiers[str(memory_tier_gb(int(r["file_size_bytes"])))] += 1
    names = ("file-manifest.tsv", "report-manifest.tsv", "resources.config", "params.json", "run-nextflow.sh")
    prov = {"scope": "PRODUCTION_LOOP_ENGINEERING_QUALIFIED_LANE", "scientific_approval": False,
            "loop_id": loop_id, "tag": tag, "source": str(source),
            "source_revision": (source / ".source-sha").read_text().strip(),
            "sif": str(sif), "sif_sha256": sif_sha256, "rules": rules, "files": len(picked),
            "accessions": len(per), "total_gb": round(sum(int(r["file_size_bytes"]) for r in picked) / 1e9, 2),
            "memory_gb_tiers": dict(sorted(tiers.items())), "allow_complete_study": False,
            "mongo_write": False, "sdrf_write": False,
            "controller_sha256": {n: sha256(ctrl / n) for n in names}, "created_utc": now()}
    (ctrl / "screen-provenance.json").write_text(json.dumps(prov, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return ctrl


def run_controller(ctrl: Path, env: dict[str, str]) -> int:
    """Create the one-shot marker, then run Nextflow in the foreground."""
    fd = os.open(str(ctrl / "launch-attempt.txt"), os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o640)
    with os.fdopen(fd, "w") as fh:
        fh.write(f"launched_utc={now()}\nhost={os.uname().nodename}\nslurm_job={os.environ.get('SLURM_JOB_ID', '')}\n")
    return subprocess.run(["bash", str(ctrl / "run-nextflow.sh")], env=env).returncode


def evaluate(root: Path, tag: str, expected: int, min_success: float) -> dict:
    ctrl = root / "work/controllers" / tag
    try:
        rc: int | None = int((ctrl / "controller.exitcode").read_text().strip())
    except (OSError, ValueError):
        rc = None
    log = (ctrl / "nextflow.log").read_text(encoding="utf-8", errors="replace") if (ctrl / "nextflow.log").is_file() else ""
    m = COMPLETED.search(log)
    completed = m is not None and int(m.group(2)) == 0
    rows = scan_tag(root / "results" / tag) if (root / "results" / tag).is_dir() else []
    keys = {(r["a"], r["f"]) for r in rows}
    ok_keys = {(r["a"], r["f"]) for r in rows if r.get("s") == "success"}
    rate = len(ok_keys) / expected if expected else 0.0
    verdict = "PASS"
    if rc != 0 and not completed:
        verdict = "STOP_CONTROLLER_FAILED"
    elif len(keys) < expected:
        verdict = "STOP_MISSING_RECEIPTS"
    elif rate < min_success:
        verdict = "STOP_LOW_SUCCESS"
    return {"tag": tag, "exitcode": rc, "workflow_completed_no_failed_tasks": completed,
            "receipts": len(keys), "success": len(ok_keys), "expected": expected,
            "success_rate": round(rate, 4), "work_left": (root / "work" / tag).exists(), "verdict": verdict}


# ---------------------------------------------------------------- loop
class Loop:
    def __init__(self, args: argparse.Namespace, *, runner: Callable[[Path, dict[str, str]], int] = run_controller,
                 sleep: Callable[[float], None] = time.sleep,
                 hours_left: Callable[[], float | None] | None = None):
        self.a = args
        self.root = Path(args.root)
        self.dir = self.root / "work/loops" / args.loop_id
        self.runner = runner
        self.sleep = sleep
        self.hours_left = hours_left or (lambda: slurm_hours_left(os.environ.get("SLURM_JOB_ID")))

    def log(self, event: str, **kw: Any) -> None:
        row = dict(kw, event=event, utc=now())
        with (self.dir / "loop.jsonl").open("a", encoding="utf-8") as fh:
            fh.write(json.dumps(row, sort_keys=True) + "\n")
        print(json.dumps(row, sort_keys=True), flush=True)

    def tag(self, n: int) -> str:
        rev = (Path(self.a.source) / ".source-sha").read_text().strip()[:7]
        return f"{self.a.tag_prefix}-b{n:03d}-{rev}-{dt.datetime.now(UTC):%Y%m%d}-v1"

    def guard(self) -> str | None:
        if (self.dir / "STOP").exists():
            return "STOP_FILE"
        if free_gb(self.root) < self.a.min_free_gb:
            return "LOW_DISK"
        left = self.hours_left()
        if left is not None and left < self.a.min_hours_left:
            return "SLURM_TIME_LOW"
        return None

    def wait_for_idle(self) -> str | None:
        deadline = time.time() + self.a.wait_active_hours * 3600
        while True:
            busy = active_controllers(self.root)
            if not busy:
                return None
            if time.time() >= deadline:
                return "ACTIVE_CONTROLLER " + ",".join(busy)
            self.log("WAIT_ACTIVE_CONTROLLER", active=busy)
            self.sleep(600)

    def cut(self, Q: Any, rows: list[dict]) -> tuple[list[dict], dict]:
        receipts = scan_results(self.root, self.dir / "receipt-cache")
        skip, retry = classify(receipts)
        busy = in_flight_keys(self.root)
        skip |= busy
        retry -= busy
        picked = Q.select(rows, max_files=self.a.max_files, max_bytes=int(self.a.max_gb * 1e9), large_files=0,
                          per_accession=self.a.per_accession, skip=skip,
                          max_file_bytes=int(self.a.max_file_gb * 1e9), max_accessions=self.a.max_accessions,
                          first=retry, lane="QUEUE_NATIVE_RAW")
        for r in picked:
            r["accession"] = r["accession"].upper()
        info = {"receipts": len(receipts), "skip_keys": len(skip), "retry_keys": len(retry), "in_flight": len(busy),
                "retry_in_batch": sum(1 for r in picked if (r["accession"], r["file_name"].casefold()) in retry),
                "files": len(picked), "total_gb": round(sum(int(r["file_size_bytes"]) for r in picked) / 1e9, 2)}
        return picked, info

    def run(self) -> str:
        self.dir.mkdir(parents=True, exist_ok=True)
        src, sif = Path(self.a.source), Path(self.a.sif)
        if sha256(sif) != self.a.sif_sha256:
            self.log("STOP", reason="SIF_HASH_MISMATCH")
            return "SIF_HASH_MISMATCH"
        Q = load_queue_module(src)
        with Path(self.a.queue).open(encoding="utf-8", newline="") as fh:
            keep = ("accession", "file_name", "file_size_bytes", "memory_gb", "decision")
            rows = [{k: r.get(k, "") for k in keep} for r in csv.DictReader(fh, delimiter="\t")
                    if r.get("decision") == "QUEUE_NATIVE_RAW"]
        for r in rows:
            r["accession"] = r["accession"].upper()
        self.log("START", loop_id=self.a.loop_id, batches=self.a.batches, first_batch=self.a.first_batch,
                 queue=self.a.queue, queue_sha256=sha256(Path(self.a.queue)), queue_rows=len(rows),
                 source=str(src), sif=str(sif), dry_run=self.a.dry_run)
        rules = {k: getattr(self.a, k) for k in ("max_files", "max_gb", "max_file_gb", "per_accession",
                                                  "max_accessions", "min_success", "max_forks")}
        for i in range(self.a.batches):
            n = self.a.first_batch + i
            reason = self.guard() or (None if self.a.dry_run else self.wait_for_idle())
            if reason:
                self.log("STOP", reason=reason, batch=n)
                return reason
            if self.a.dry_run and active_controllers(self.root):
                self.log("DRY_RUN_NOTE_ACTIVE_CONTROLLER", active=active_controllers(self.root))
            picked, info = self.cut(Q, rows)
            tag = self.tag(n)
            self.log("CUT", batch=n, tag=tag, **info)
            if len(picked) < self.a.min_files:
                self.log("STOP", reason="QUEUE_EXHAUSTED_OR_SMALL_BATCH", batch=n, files=len(picked))
                return "QUEUE_EXHAUSTED_OR_SMALL_BATCH"
            if self.a.dry_run:
                self.log("DRY_RUN_STOP", batch=n)
                return "DRY_RUN"
            ctrl = write_controller(self.root, tag, picked, source=src, sif=sif, sif_sha256=self.a.sif_sha256,
                                    rules=rules, nxf_opts=self.a.nxf_opts, max_forks=self.a.max_forks,
                                    loop_id=self.a.loop_id)
            self.log("LAUNCH", batch=n, tag=tag, controller=str(ctrl))
            start = time.time()
            rc = self.runner(ctrl, dict(os.environ))
            if not (ctrl / "controller.exitcode").is_file():   # runner shell killed: record it, do not hang
                (ctrl / "controller.exitcode").write_text(f"{rc}\n", encoding="utf-8")
            result = evaluate(self.root, tag, len(picked), self.a.min_success)
            result["hours"] = round((time.time() - start) / 3600, 2)
            result["gb_per_min"] = round(info["total_gb"] / max((time.time() - start) / 60, 1e-9), 1)
            self.log("RESULT", batch=n, **result)
            if result["verdict"] != "PASS":
                self.log("STOP", reason=result["verdict"], batch=n)
                return result["verdict"]
        self.log("DONE", batches=self.a.batches)
        return "DONE"


def parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", required=True)
    ap.add_argument("--loop-id", required=True)
    ap.add_argument("--source", required=True, help="Pinned source dir with .source-sha")
    ap.add_argument("--sif", required=True)
    ap.add_argument("--sif-sha256", required=True)
    ap.add_argument("--queue", required=True, help="Inventory queue.tsv")
    ap.add_argument("--batches", type=int, required=True)
    ap.add_argument("--first-batch", type=int, required=True)
    ap.add_argument("--tag-prefix", default="prideqc-prod-raw")
    ap.add_argument("--max-files", type=int, default=10000)
    ap.add_argument("--max-gb", type=float, default=10000.0)
    ap.add_argument("--max-file-gb", type=float, default=6.0)
    ap.add_argument("--per-accession", type=int, default=200)
    ap.add_argument("--max-accessions", type=int, default=300)
    ap.add_argument("--min-files", type=int, default=1000)
    ap.add_argument("--min-success", type=float, default=0.95)
    ap.add_argument("--max-forks", type=int, default=120)
    ap.add_argument("--nxf-opts", default="-Xms1g -Xmx4g")
    ap.add_argument("--min-free-gb", type=float, default=500.0)
    ap.add_argument("--min-hours-left", type=float, default=8.0)
    ap.add_argument("--wait-active-hours", type=float, default=0.0,
                    help="Wait up to N hours for another active controller to finish")
    ap.add_argument("--dry-run", action="store_true",
                    help="Scan and cut the first batch; write only the loop log and receipt cache; launch nothing")
    return ap


def main() -> int:
    args = parser().parse_args()
    if not 1 <= args.batches <= 40:
        print(json.dumps({"loop": "BLOCKED", "reason": "batches must be 1..40"}))
        return 2
    try:
        result = Loop(args).run()
    except (OSError, ValueError, subprocess.SubprocessError) as err:
        print(json.dumps({"loop": "BLOCKED", "reason": f"{type(err).__name__}: {err}"}))
        return 2
    return 0 if result in ("DONE", "DRY_RUN") else 3


if __name__ == "__main__":
    sys.exit(main())
