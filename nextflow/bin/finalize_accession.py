#!/usr/bin/env python3
"""Write accession success/failure audit and merge ONLY fully successful logical runs.

A pipeline completing with failed unit receipts is NOT an accession QC success.
Consumers must gate on accession-status.json:status == complete.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
import subprocess
from collections import Counter
from pathlib import Path


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def finalize(accession: str, expected_units: int, root: Path, outdir: Path, *, verified_complete_inventory: bool = False) -> dict:
    if expected_units < 1:
        raise ValueError("expected_units must be positive")
    if outdir.exists() and any(outdir.iterdir()):
        raise ValueError("Refusing nonempty accession output")
    outdir.mkdir(parents=True, exist_ok=True)
    audit = outdir / "_audit"
    audit.mkdir()
    receipts = sorted(root.rglob("result.json"))
    results: list[tuple[dict, Path]] = []
    for receipt in receipts:
        item = json.loads(receipt.read_text(encoding="utf-8"))
        if item.get("accession") != accession:
            raise ValueError("Cross-accession receipt encountered")
        results.append((item, receipt.parent))
    ids = [item["unit_id"] for item, _ in results]
    tasks = {item["task_id"] for item, _ in results}
    status_counts = Counter(item["status"] for item, _ in results)
    expected_archives = {item.get("expected_file_tasks") for item, _ in results}
    archive_count_ok = (len(expected_archives) == 1 and
                        len(tasks) == next(iter(expected_archives)))
    complete = (archive_count_ok and len(results) == expected_units and len(ids) == len(set(ids))
                and status_counts.get("success", 0) == expected_units
                and status_counts.get("failed", 0) == 0)
    failures = []
    if not archive_count_ok:
        failures.append("source_archive_task_count_mismatch")
    if len(results) != expected_units:
        failures.append("missing_or_extra_unit_receipts")
    if len(ids) != len(set(ids)):
        failures.append("duplicate_unit_ids")
    if not complete and not failures:
        failures.append("one_or_more_units_failed")
    for item, _ in results:
        if item["status"] != "success":
            failures.append(f"unit_failed:{item['unit_id']}")
    with (audit / "file-status.tsv").open("w", encoding="utf-8", newline="") as fh:
        columns = ("accession", "task_id", "source_archive", "source_member",
                   "unit_id", "status", "failure_stage", "error_type", "exit_code",
                   "mzqc_sha256")
        writer = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for item, _ in results:
            writer.writerow({key: item.get(key, "") for key in columns})
    if complete and verified_complete_inventory:
        source = outdir / "_merge_inputs"
        source.mkdir()
        try:
            seen: set[str] = set()
            for item, folder in results:
                qc = folder / "qc"
                name = item["mzqc_file"]
                obo_name = item["obo_file"]
                if name.casefold() in seen or obo_name != Path(name).with_suffix(".obo").name:
                    raise ValueError("Ambiguous per-run output basename")
                seen.add(name.casefold())
                for basename, expected in ((name, item["mzqc_sha256"]),
                                            (obo_name, item["obo_sha256"])):
                    path = next(iter(qc.rglob(basename)), None)
                    if path is None or digest(path) != expected:
                        raise ValueError("Missing or changed per-run QC output")
                    shutil.copy2(path, source / basename)
            study = outdir / "_study"
            study.mkdir()
            cmd = ["prideqc", "mzqc", "merge", "--input-dir", str(source),
                   "--output", str(study / f"{accession}.mzQC"),
                   "--accession", accession, "--expected-runs", str(expected_units)]
            proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
            if proc.returncode:
                raise RuntimeError(f"Study merge exited {proc.returncode}")
            merged_manifest = json.loads((study / f"{accession}.manifest.json").read_text())
            if merged_manifest.get("status") != "complete" or merged_manifest.get("merged_runs") != expected_units:
                raise ValueError("Study merge count or status mismatch")
        except (OSError, ValueError, RuntimeError, KeyError) as exc:
            failures.append(f"merge_failed:{type(exc).__name__}")
            complete = False
            shutil.rmtree(outdir / "_study", ignore_errors=True)
        finally:
            shutil.rmtree(source, ignore_errors=True)
    state = "incomplete" if not complete else ("complete" if verified_complete_inventory else "partial")
    summary = {"schema_version": 2, "accession": accession,
               "status": state,
               "inventory_scope": "operator_verified_complete" if verified_complete_inventory else "pilot_subset_unverified",
               "study_merge_performed": state == "complete",
               "expected_logical_runs": expected_units,
               "observed_receipts": len(results),
               "successful_runs": status_counts.get("success", 0),
               "failed_runs": status_counts.get("failed", 0),
               "source_archive_tasks_observed": len(tasks),
               "expected_source_archive_tasks": next(iter(expected_archives)) if len(expected_archives) == 1 else None,
               "failures": sorted(failures), "mongo_write_performed": False}
    (audit / "accession-status.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    return summary


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--accession", required=True)
    ap.add_argument("--expected-units", type=int, required=True)
    ap.add_argument("--root", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--verified-complete-inventory", action="store_true")
    args = ap.parse_args()
    print(json.dumps(finalize(args.accession, args.expected_units, args.root, args.outdir,
                              verified_complete_inventory=args.verified_complete_inventory), sort_keys=True))


if __name__ == "__main__":
    main()
