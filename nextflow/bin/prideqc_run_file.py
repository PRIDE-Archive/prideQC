#!/usr/bin/env python3
"""One-file execution adapter. Never publishes a failed or partial mzQC."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
from pathlib import Path


def check_fire_supported() -> None:
    try:
        from pridepy.pridepy import PROTOCOL_CHOICES

        choices = getattr(PROTOCOL_CHOICES, "choices", [])
    except Exception as exc:
        raise RuntimeError(f"Cannot inspect installed pridepy FIRE capability: {exc}") from exc
    if "fire" not in choices:
        raise RuntimeError(
            "Installed pridepy does not support --protocol fire. Pin/build the FIRE-capable pridepy dev revision before running Codon."
        )


def run(
    task_file: Path,
    protocol: str,
    output: Path,
    download: Path,
    checksum: bool,
    mass_error: bool,
    mass_shifts: bool,
    tmt_reporter: bool = False,
    tmt_profile: bool = False,
) -> None:
    task = json.loads(task_file.read_text(encoding="utf-8"))
    pxd, name = task["accession"], task["archive_file"]
    if protocol not in {"ftp", "aspera", "s3", "globus", "fire"}:
        raise ValueError(f"Unsupported protocol: {protocol}")
    if protocol == "fire":
        check_fire_supported()
    output.mkdir(parents=True, exist_ok=True)
    download.mkdir(parents=True, exist_ok=True)
    if any(output.iterdir()) or any(download.iterdir()):
        raise ValueError("Refusing nonempty task output/download directory")
    cmd = [
        "prideqc",
        "analyze",
        "--accession",
        pxd,
        "--file",
        name,
        "--protocol",
        protocol,
        "--download-dir",
        str(download),
        "--workers",
        "1",
        "--no-progress",
        "-o",
        str(output),
    ]
    if not checksum:
        cmd.append("--no-checksum-check")
    if mass_error:
        cmd.append("--estimate-mass-error")
    if mass_shifts:
        cmd.append("--estimate-mass-shifts")
    if tmt_reporter:
        cmd.append("--tmt-reporter-evidence")
    if tmt_profile:
        cmd.append("--tmt-profile-reporter-evidence")
    subprocess.run(cmd, check=True)
    manifest_file = output / "manifest.json"
    if not manifest_file.is_file():
        raise RuntimeError("Missing prideQC per-run manifest")
    manifest = json.loads(manifest_file.read_text(encoding="utf-8"))
    records = manifest.get("files", [])
    if len(records) != 1 or records[0].get("status") != "success" or records[0].get("error"):
        raise RuntimeError("Per-run manifest does not report precisely one successful file")
    mzqcs = list(output.rglob("*.mzQC"))
    if len(mzqcs) != 1 or mzqcs[0].name != f"{name}.mzQC":
        raise RuntimeError(f"Expected one {name}.mzQC, found {mzqcs}")
    obo = mzqcs[0].with_suffix(".obo")
    if not obo.is_file():
        raise RuntimeError("Missing OBO sidecar")
    (output / "nextflow-run.json").write_text(
        json.dumps(
            {
                "task_id": task["task_id"],
                "accession": pxd,
                "archive_file": name,
                "expected_runs": task["expected_runs"],
                "requested_protocol": protocol,
                "mzqc_sha256": hashlib.sha256(mzqcs[0].read_bytes()).hexdigest(),
                "obo_sha256": hashlib.sha256(obo.read_bytes()).hexdigest(),
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    # Downloads are not published and are removed only after a verified successful run.
    if download.resolve().parent != Path.cwd().resolve():
        raise RuntimeError("Refusing to clean a download directory outside this task workdir")
    shutil.rmtree(download)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--task", type=Path, required=True)
    ap.add_argument("--protocol", required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--download", type=Path, required=True)
    ap.add_argument("--no-checksum-check", action="store_true")
    ap.add_argument("--estimate-mass-error", action="store_true")
    ap.add_argument("--estimate-mass-shifts", action="store_true")
    ap.add_argument("--tmt-reporter-evidence", action="store_true")
    ap.add_argument("--tmt-profile-reporter-evidence", action="store_true")
    a = ap.parse_args()
    run(
        a.task,
        a.protocol,
        a.output,
        a.download,
        not a.no_checksum_check,
        a.estimate_mass_error,
        a.estimate_mass_shifts,
        a.tmt_reporter_evidence,
        a.tmt_profile_reporter_evidence,
    )


if __name__ == "__main__":
    main()
