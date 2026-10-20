#!/usr/bin/env python3
"""Pilot-only: download/expand/analyze ONE logical run in one scratch task.

Only the verified receipt, mzQC, OBO and bounded logs may leave node-local
scratch. A multi-run bundle is explicitly unsupported in this pilot lane.
"""
from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

from analyze_unit import analyze
from prepare_archive import prepare


def run(task_file: Path, protocol: str, scratch_dir: Path, outdir: Path,
        *, source: Path | None = None, checksum: bool = True,
        mass_error: bool = False, mass_shifts: bool = False) -> dict:
    if scratch_dir.exists() or (outdir.exists() and any(outdir.iterdir())):
        raise ValueError("Refusing existing scratch or nonempty output directory")
    if scratch_dir.resolve().parent != Path.cwd().resolve():
        raise ValueError("Source scratch directory must be directly under task-local cwd")
    outdir.mkdir(parents=True, exist_ok=True)
    ready_units = outdir / "units"
    ready_units.mkdir()
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
                         mass_shifts=mass_shifts)
        # Record the requested scientific estimators and the provenance of a
        # curated primary representation even if the actual reader/QC failed.
        # These are intent/provenance markers, NOT evidence that metrics passed.
        task = json.loads(task_file.read_text(encoding="utf-8"))
        result["requested_estimators"] = {"mass_error": bool(mass_error),
                                          "mass_shifts": bool(mass_shifts)}
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
    args = ap.parse_args()
    print(json.dumps(run(args.task, args.protocol, args.scratch_dir, args.outdir,
                         source=args.local_source, checksum=not args.no_checksum_check,
                         mass_error=args.estimate_mass_error,
                         mass_shifts=args.estimate_mass_shifts), sort_keys=True))


if __name__ == "__main__":
    main()
