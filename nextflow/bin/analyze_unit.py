#!/usr/bin/env python3
"""Analyze one logical MS run; always write an explicit unit success/failure receipt."""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from pathlib import Path


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def analyze(unit_dir: Path, output: Path, *, mass_error: bool = False,
            mass_shifts: bool = False) -> dict:
    unit = json.loads((unit_dir / "unit.json").read_text(encoding="utf-8"))
    if output.exists() and any(output.iterdir()):
        raise ValueError("Refusing nonempty unit result directory")
    output.mkdir(parents=True, exist_ok=True)
    result = {"schema_version": 2, "accession": unit["accession"],
              "task_id": unit["task_id"], "source_archive": unit["source_archive"],
              "unit_id": unit["unit_id"], "source_member": unit.get("source_member"),
              "expected_file_tasks": unit["expected_file_tasks"], "status": "failed",
              "failure_stage": None, "error_type": None}
    try:
        if unit["discovery_status"] != "success":
            result["failure_stage"] = unit.get("failure_stage", "archive_preparation")
            result["error_type"] = unit.get("error_type", "DiscoveryFailed")
        else:
            source = unit_dir / unit["input_basename"]
            if source.is_symlink() or not source.exists():
                raise ValueError("Missing or symlinked logical run input")
            qc = output / "qc"
            cmd = ["prideqc", "analyze", str(source), "-o", str(qc),
                   "--project-accession", unit["accession"], "--no-progress"]
            if mass_error:
                cmd.append("--estimate-mass-error")
            if mass_shifts:
                cmd.append("--estimate-mass-shifts")
            proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
            # Keep bounded local logs; do not include full command, credentials, or URLs.
            (output / "prideqc.stdout.log").write_text(proc.stdout[-16000:], encoding="utf-8")
            (output / "prideqc.stderr.log").write_text(proc.stderr[-16000:], encoding="utf-8")
            if proc.returncode:
                result["failure_stage"] = "prideqc_analysis"
                result["error_type"] = "PrideqcNonzeroExit"
                result["exit_code"] = proc.returncode
                # Failure triage only: a .raw suffix does not prove Thermo
                # content, and a reader-open error does not prove corruption.
                if "Thermo RAW reader: Could not open RAW file" in proc.stderr:
                    result["reader_failure_kind"] = "thermo_raw_open_error_unclassified"
            else:
                manifest = json.loads((qc / "manifest.json").read_text())
                records = manifest.get("files", [])
                files = list(qc.rglob("*.mzQC"))
                if (len(records) != 1 or records[0].get("status") != "success"
                        or records[0].get("error") or len(files) != 1):
                    raise ValueError("QC output integrity or per-run status check failed")
                mzqc = files[0]
                obo = mzqc.with_suffix(".obo")
                if not obo.is_file():
                    raise ValueError("Missing per-run OBO")
                result.update({"status": "success", "failure_stage": None,
                               "mzqc_file": mzqc.name, "mzqc_sha256": digest(mzqc),
                               "obo_file": obo.name, "obo_sha256": digest(obo)})
    except (OSError, ValueError, KeyError, TypeError) as exc:
        result["failure_stage"] = result["failure_stage"] or "qc_validation"
        result["error_type"] = type(exc).__name__
    finally:
        (output / "result.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--unit", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--estimate-mass-error", action="store_true")
    ap.add_argument("--estimate-mass-shifts", action="store_true")
    args = ap.parse_args()
    print(json.dumps(analyze(args.unit, args.output,
                             mass_error=args.estimate_mass_error,
                             mass_shifts=args.estimate_mass_shifts), sort_keys=True))


if __name__ == "__main__":
    main()
