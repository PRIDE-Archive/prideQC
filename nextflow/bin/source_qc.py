#!/usr/bin/env python3
"""Pilot-only: download/expand/analyze ONE logical run in one scratch task.

Only the verified receipt, mzQC, OBO and bounded logs may leave node-local
scratch. A multi-run bundle is explicitly unsupported in this pilot lane.
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import signal
import subprocess
import sys
from pathlib import Path

from analyze_unit import analyze
from prepare_archive import prepare

ARCHIVE_MEMBERS = "archive_members"
MEMBER_PLAN = "member-plan.json"


def run(task_file: Path, protocol: str, scratch_dir: Path, outdir: Path,
        *, source: Path | None = None, checksum: bool = True,
        mass_error: bool = False, mass_shifts: bool = False,
        tmt_reporter: bool = False, tmt_profile: bool = False) -> dict:
    if scratch_dir.exists() or (outdir.exists() and any(outdir.iterdir())):
        raise ValueError("Refusing existing scratch or nonempty output directory")
    if scratch_dir.resolve().parent != Path.cwd().resolve():
        raise ValueError("Source scratch directory must be directly under task-local cwd")
    outdir.mkdir(parents=True, exist_ok=True)
    ready_units = outdir / "units"
    ready_units.mkdir()
    task = json.loads(task_file.read_text(encoding="utf-8"))
    if task.get("source_kind") == ARCHIVE_MEMBERS:
        return run_archive(task, protocol, scratch_dir, outdir, source=source, checksum=checksum,
                           mass_error=mass_error, mass_shifts=mass_shifts, tmt_reporter=tmt_reporter,
                           tmt_profile=tmt_profile)
    try:
        prepared = scratch_dir / "prepared"
        status = prepare(task_file, protocol, prepared, source=source, checksum=checksum)
        units = sorted((prepared / "units").iterdir())
        if len(units) != 1 or status["status"] != "success":
            if status["status"] == "success":
                status["status"] = "failed"
                status["error_type"] = "MultipleLogicalRunsUnsupported"
                status["discovered_logical_runs"] = len(units)
                status["emitted_units"] = 1
            # For genuine discovery failures the existing failure sentinel is reused.
            # For multi-run archives, create one synthetic failure instead of silently
            # dropping units or processing a misleading subset of the archive.
            if len(units) == 1 and json.loads((units[0] / "unit.json").read_text()).get("discovery_status") == "failed":
                unit = units[0]
            else:
                task = json.loads(task_file.read_text(encoding="utf-8"))
                unit = scratch_dir / "pilot-failure"
                unit.mkdir()
                (unit / "unit.json").write_text(json.dumps({
                    "schema_version": 2, "accession": task["accession"],
                    "task_id": task["task_id"], "source_archive": task["archive_file"],
                    "expected_file_tasks": task["expected_runs"],
                    "unit_id": f"archive-{int(task['task_id']):09d}-pilot-failure",
                    "discovery_status": "failed", "failure_stage": "archive_preparation",
                    "error_type": status["error_type"] or "InvalidUnitCount",
                }, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        else:
            unit = units[0]
        unit_id = json.loads((unit / "unit.json").read_text(encoding="utf-8"))["unit_id"]
        result = analyze(unit, ready_units / unit_id, mass_error=mass_error,
                         mass_shifts=mass_shifts, tmt_reporter=tmt_reporter,
                         tmt_profile=tmt_profile)
        # Record the requested scientific estimators and the provenance of a
        # curated primary representation even if the actual reader/QC failed.
        # These are intent/provenance markers, NOT evidence that metrics passed.
        task = json.loads(task_file.read_text(encoding="utf-8"))
        result["requested_estimators"] = {"mass_error": bool(mass_error),
                                          "mass_shifts": bool(mass_shifts),
                                          "tmt_reporter_evidence": bool(tmt_reporter),
                                          "tmt_profile_reporter_evidence": bool(tmt_profile)}
        if task.get("representation_policy"):
            result["source_representation"] = {
                "logical_run_id": task.get("logical_run_id"),
                "policy": task["representation_policy"],
                "inventory_sha256": task.get("representation_inventory_sha256"),
            }
        (ready_units / unit_id / "result.json").write_text(
            json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        (outdir / "prepare-status.json").write_text(
            json.dumps(status, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        return result
    finally:
        # PREPARE's download, extracted RAW, and all RAW-containing unit paths
        # never become declared Nextflow outputs, even if analysis fails.
        shutil.rmtree(scratch_dir, ignore_errors=True)


def run_archive(task: dict, protocol: str, scratch_dir: Path, outdir: Path, *, source: Path | None,
                checksum: bool, mass_error: bool, mass_shifts: bool, tmt_reporter: bool,
                tmt_profile: bool) -> dict:
    """Archive lane: one receipt per run member (archive_member_qc.py)."""
    from archive_member_qc import run_members

    requested = {"mass_error": bool(mass_error), "mass_shifts": bool(mass_shifts),
                 "tmt_reporter_evidence": bool(tmt_reporter), "tmt_profile_reporter_evidence": bool(tmt_profile)}

    def analyzer(unit: Path, output: Path, **kw: object) -> dict:
        result = analyze(unit, output, mass_error=mass_error, mass_shifts=mass_shifts,
                         tmt_reporter=tmt_reporter, tmt_profile=tmt_profile, **kw)  # type: ignore[arg-type]
        result["requested_estimators"] = requested
        return result

    try:
        status = run_members(task, protocol, scratch_dir, outdir / "units", source=source, checksum=checksum,
                             analyzer=analyzer, plan_file=outdir / MEMBER_PLAN)
    finally:
        shutil.rmtree(scratch_dir, ignore_errors=True)
    (outdir / "prepare-status.json").write_text(json.dumps(status, indent=2, sort_keys=True) + "\n",
                                                 encoding="utf-8")
    return status


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--task", required=True, type=Path)
    ap.add_argument("--protocol", required=True, choices=("ftp", "aspera", "s3", "globus", "fire"))
    ap.add_argument("--scratch-dir", required=True, type=Path)
    ap.add_argument("--outdir", required=True, type=Path)
    ap.add_argument("--local-source", type=Path, help="Only for offline fixtures")
    ap.add_argument("--no-checksum-check", action="store_true")
    ap.add_argument("--estimate-mass-error", action="store_true")
    ap.add_argument("--estimate-mass-shifts", action="store_true")
    ap.add_argument("--tmt-reporter-evidence", action="store_true")
    ap.add_argument("--tmt-profile-reporter-evidence", action="store_true")
    ap.add_argument("--attempt", type=int, default=1, help="Nextflow task.attempt")
    ap.add_argument("--max-attempts", type=int, default=1, help="Nextflow maxRetries + 1")
    ap.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    args = ap.parse_args()
    if not args.worker:
        raise SystemExit(supervise(args, sys.argv[1:]))
    # Make the worker (and the prideqc child it starts) the first OOM victim, never the
    # small supervisor. Raising the score needs no privilege; best effort.
    try:
        Path("/proc/self/oom_score_adj").write_text("1000\n")
    except OSError:
        pass
    result = run(args.task, args.protocol, args.scratch_dir, args.outdir,
                 source=args.local_source, checksum=not args.no_checksum_check,
                 mass_error=args.estimate_mass_error,
                 mass_shifts=args.estimate_mass_shifts,
                 tmt_reporter=args.tmt_reporter_evidence,
                 tmt_profile=args.tmt_profile_reporter_evidence)
    print(json.dumps(result, sort_keys=True))
    raise SystemExit(retry_exit_code(result, args.attempt, args.max_attempts))


def retry_exit_code(result: dict, attempt: int, max_attempts: int) -> int:
    """137 asks Nextflow to retry with more memory; only for a resource kill before the last attempt.
    On the last attempt the fail-closed receipt stands and the task exits 0.
    Archive lane: the summary has no error_type of its own; a killed member counts."""
    if result.get("error_type") == "AnalysisKilledResourceLimit" and attempt < max_attempts:
        return 137
    if result.get("members_killed") and attempt < max_attempts:
        return 137
    return 0


def killed_receipt(task_file: Path, scratch_dir: Path, outdir: Path, signal_number: int) -> dict:
    """Fail-closed receipt for a worker that a signal killed (no receipt of its own).

    Batch 041 (2026-11-15): PXD072437/2zup_RCC12F3.raw (0.8 GB) exceeded 16 GB, then 32 GB.
    On attempt 2 the Slurm memory cgroup killed the worker python too, so no receipt was
    written, the task exited 137 on its last attempt and Nextflow aborted the batch.
    """
    task = json.loads(task_file.read_text(encoding="utf-8"))
    # The worker could not run its own cleanup: RAW must not stay in task scratch.
    shutil.rmtree(scratch_dir, ignore_errors=True)
    if task.get("source_kind") == ARCHIVE_MEMBERS and (outdir / MEMBER_PLAN).is_file():
        return killed_archive_receipts(task, outdir, signal_number)
    shutil.rmtree(outdir, ignore_errors=True)
    unit_id = f"archive-{int(task['task_id']):09d}-worker-killed"
    unit = outdir / "units" / unit_id
    unit.mkdir(parents=True)
    name = signal.Signals(signal_number).name if signal_number in signal.valid_signals() else str(signal_number)
    result = {"schema_version": 2, "accession": task["accession"], "task_id": task["task_id"],
              "source_archive": task["archive_file"], "source_member": None, "unit_id": unit_id,
              "expected_file_tasks": task["expected_runs"], "status": "failed",
              "failure_stage": "prideqc_analysis", "error_type": "AnalysisKilledResourceLimit",
              "resource_failure_kind": "task_killed_by_signal", "exit_code": -signal_number,
              "worker_signal": name}
    (unit / "result.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    status = {"schema_version": 2, "accession": task["accession"], "task_id": task["task_id"],
              "source_archive": task["archive_file"], "status": "unknown_worker_killed",
              "error_type": "WorkerKilled", "emitted_units": 1,
              "source_manifest": task.get("source_manifest"),
              "source_manifest_task_count": task["expected_runs"], "source_evidence": {}}
    (outdir / "prepare-status.json").write_text(json.dumps(status, indent=2, sort_keys=True) + "\n",
                                                 encoding="utf-8")
    return result


def killed_archive_receipts(task: dict, outdir: Path, signal_number: int) -> dict:
    """Archive lane: keep the receipts of finished members; fail closed for the others.

    The member in progress (result folder without result.json) gets a resource-kill receipt.
    A member that did not start gets NotAttemptedAfterWorkerKill. So every planned member
    has exactly one receipt and no member is lost.
    """
    plan = json.loads((outdir / MEMBER_PLAN).read_text(encoding="utf-8"))
    units = outdir / "units"
    name = signal.Signals(signal_number).name if signal_number in signal.valid_signals() else str(signal_number)
    killed = 0
    for entry in plan["members"]:
        unit = units / entry["unit_id"]
        if (unit / "result.json").is_file():
            continue
        started = unit.is_dir()
        shutil.rmtree(unit, ignore_errors=True)
        unit.mkdir(parents=True)
        result = dict(plan["base"], unit_id=entry["unit_id"], source_member=entry["member"],
                      member_index=entry["member_index"], member_count=len(plan["members"]), status="failed",
                      failure_stage="prideqc_analysis" if started else "archive_member_extraction",
                      error_type="AnalysisKilledResourceLimit" if started else "NotAttemptedAfterWorkerKill",
                      worker_signal=name)
        if started:
            result.update(resource_failure_kind="task_killed_by_signal", exit_code=-signal_number)
            killed += 1
        (unit / "result.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    status = dict(plan["base"], status="unknown_worker_killed", error_type="WorkerKilled",
                  emitted_units=len(plan["members"]), members_killed=killed,
                  source_manifest=task.get("source_manifest"), source_manifest_task_count=task["expected_runs"])
    (outdir / "prepare-status.json").write_text(json.dumps(status, indent=2, sort_keys=True) + "\n",
                                                 encoding="utf-8")
    return status


def supervise(args: argparse.Namespace, argv: list[str],
              runner: object = subprocess.run) -> int:
    """Run the work in a child process; this small parent survives a memory kill of the child.

    A signal kill before the last attempt exits 137 (Nextflow retries with 2x memory).
    On the last attempt the parent writes a fail-closed receipt and exits 0, so one bad
    file never aborts the batch.
    """
    cmd = [sys.executable, os.path.abspath(__file__), *argv, "--worker"]
    proc = runner(cmd, check=False)  # type: ignore[operator]
    code = proc.returncode
    if code >= 0:
        return code
    if args.attempt < args.max_attempts:
        shutil.rmtree(args.scratch_dir, ignore_errors=True)
        return 137
    result = killed_receipt(args.task, args.scratch_dir, args.outdir, -code)
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    main()
