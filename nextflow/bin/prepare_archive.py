#!/usr/bin/env python3
"""Acquire one manifest archive using pridepy, then emit logical run folders.

Scientific/discovery failures emit a fail-closed synthetic unit rather than
silently reducing the expected study run count. Runtime termination still fails
and is managed by Nextflow's retry policy.
"""
from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

from archive_units import ArchiveError, expand


def prepare(task_file: Path, protocol: str, outdir: Path,
            *, source: Path | None = None, checksum: bool = True) -> dict:
    task = json.loads(task_file.read_text(encoding="utf-8"))
    pxd = task["accession"]
    task_id = int(task["task_id"])
    archive_name = task["archive_file"]
    if outdir.exists() and any(outdir.iterdir()):
        raise ArchiveError("Refusing nonempty preparation output")
    outdir.mkdir(parents=True, exist_ok=True)
    units = outdir / "units"
    units.mkdir()
    transfer = outdir / "download"
    status = "success"
    error_type = None
    try:
        if source is None:
            from prideqc.pride import DownloadOptions, PrideRepository

            if protocol == "fire":
                from pridepy.pridepy import PROTOCOL_CHOICES
                choices = getattr(PROTOCOL_CHOICES, "choices", PROTOCOL_CHOICES)
                if "fire" not in choices:
                    raise RuntimeError("FIRE capability missing in installed pridepy")
            paths = PrideRepository().download_files(
                pxd, [archive_name], transfer,
                options=DownloadOptions(protocol=protocol, checksum_check=checksum),
            )
            archive = paths[0]
        else:
            archive = source.resolve(strict=True)
        # Use a temporary work area outside units; only final unit folders stage to workers.
        expanded = outdir / "extracted"
        inventory = expand(archive, expanded)
        for folder in sorted(expanded.glob("unit-*")):
            target = units / f"archive-{task_id:09d}-{folder.name}"
            shutil.move(str(folder), str(target))
            info_file = target / "unit.json"
            info = json.loads(info_file.read_text())
            info.update({"schema_version": 2, "accession": pxd,
                         "unit_id": target.name,
                         "source_archive": archive_name, "task_id": task_id,
                         "expected_file_tasks": task["expected_runs"],
                         "discovery_status": "success"})
            info_file.write_text(json.dumps(info, indent=2, sort_keys=True) + "\n")
        shutil.rmtree(expanded, ignore_errors=True)
        count = inventory["logical_run_count"]
        if count < 1:
            raise ArchiveError("Archive contained zero logical runs")
    except Exception as exc:
        # Never serialize upstream download exceptions; they may contain private URLs.
        status = "failed"
        error_type = type(exc).__name__
        # Drop any partially materialized units: exactly one failure sentinel remains.
        shutil.rmtree(units, ignore_errors=True)
        units.mkdir()
        folder = units / f"archive-{task_id:09d}-failure"
        folder.mkdir()
        (folder / "unit.json").write_text(json.dumps({
            "schema_version": 2, "accession": pxd, "task_id": task_id,
            "source_archive": archive_name, "expected_file_tasks": task["expected_runs"],
            "discovery_status": "failed", "failure_stage": "archive_preparation",
            "error_type": error_type,
            "unit_id": f"archive-{task_id:09d}-failure",
        }, indent=2, sort_keys=True) + "\n")
        count = 1
    summary = {"schema_version": 2, "accession": pxd, "task_id": task_id,
               "source_archive": archive_name, "status": status,
               "error_type": error_type, "emitted_units": count,
               "source_manifest": task.get("source_manifest"),
               "source_manifest_task_count": task["expected_runs"]}
    (outdir / "prepare-status.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    shutil.rmtree(transfer, ignore_errors=True)
    return summary


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--task", type=Path, required=True)
    ap.add_argument("--protocol", choices=("ftp", "aspera", "s3", "globus", "fire"), required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--local-source", type=Path)
    ap.add_argument("--no-checksum-check", action="store_true")
    args = ap.parse_args()
    result = prepare(args.task, args.protocol, args.outdir,
                     source=args.local_source, checksum=not args.no_checksum_check)
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
