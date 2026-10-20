#!/usr/bin/env python3
"""Analyze one logical MS run; always write an explicit unit success/failure receipt."""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
from pathlib import Path


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class EmptySpectrumRun(ValueError):
    """The reader returned no MS spectra; a QC receipt must not claim success."""


def require_observed_spectra(path: Path) -> int:
    """Check authoritative mzQC scan counts, including MS3+ acquisitions.

    Uses the PSI-MS MS1/MS2 accessions and PRIDEQC local names for higher
    levels. Do not infer scan counts from an mzQC metric count or file size.
    """
    document = json.loads(path.read_text(encoding="utf-8"))
    run_qualities = document.get("mzQC", {}).get("runQualities")
    if not isinstance(run_qualities, list) or len(run_qualities) != 1:
        raise ValueError("Expected exactly one mzQC runQuality")
    quality_metrics = run_qualities[0].get("qualityMetrics")
    if not isinstance(quality_metrics, list):
        raise ValueError("Missing mzQC qualityMetrics")
    counts: dict[int, int] = {}
    for item in quality_metrics:
        if not isinstance(item, dict):
            raise ValueError("Malformed mzQC qualityMetric")
        accession, name = item.get("accession"), item.get("name")
        level = {"MS:4000059": 1, "MS:4000060": 2}.get(accession)
        match = re.fullmatch(r"NumberOfSpectra_MS([3-9][0-9]*)", name or "") if isinstance(name, str) else None
        if level is None and match:
            level = int(match.group(1))
        if level is None:
            continue
        value = item.get("value")
        if (level in counts or type(value) is not int or value < 0):
            raise ValueError(f"Invalid or duplicate MS{level} spectrum count")
        counts[level] = value
    if 1 not in counts or 2 not in counts:
        raise ValueError("Missing mzQC MS1/MS2 spectrum counts")
    total = sum(counts.values())
    if total == 0:
        raise EmptySpectrumRun("No MS spectra; zero-scan mzQC is diagnostic, not a successful QC run")
    return total


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
                # A successful CLI exit can still contain zero read spectra.
                # Fail closed so the accession finalizer cannot count an empty
                # acquisition as a scientifically usable successful run.
                require_observed_spectra(mzqc)
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
