#!/usr/bin/env python3
"""Fail-closed, deterministic conversion of frozen Slurm manifests to Nextflow tasks."""

from __future__ import annotations

import argparse
import csv
import json
import re
from collections import defaultdict
from pathlib import Path

ACCESSION = re.compile(r"PXD[0-9]{6,}$")
# An explicit, narrow v1 name policy. No paths, traversal, shell syntax or whitespace.
FILENAME = re.compile(r"[A-Za-z0-9_.+\-]+$")
SUPPORTED = (".raw", ".mzml", ".mzml.gz", ".zip", ".tar.gz", ".tgz", ".d", ".raw.gz")
MIXED_RUN_ID = re.compile(r"[A-Za-z0-9_.+\-]{1,128}$")
HEX_SHA256 = re.compile(r"[0-9a-f]{64}$")


def read_tsv(path: Path, required: set[str]) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path}: missing columns: {sorted(missing)}")
        rows = list(reader)
    if not rows:
        raise ValueError(f"Empty manifest: {path}")
    return rows


def prepare(file_manifest: Path, report_manifest: Path, output: Path, accessions: set[str]) -> dict:
    files = read_tsv(
        file_manifest, {"task_id", "pxd_accession", "archive_file", "sdrf_path", "needs_file_map"}
    )
    reports = read_tsv(report_manifest, {"report_id", "pxd_accession", "expected_file_tasks"})
    wanted = {x.upper() for x in accessions}
    selected = [r for r in files if not wanted or r["pxd_accession"].upper() in wanted]
    if not selected:
        raise ValueError("No selected RAW task rows")
    by_accession: dict[str, list[dict[str, str]]] = defaultdict(list)
    task_ids: set[str] = set()
    names_by_pxd: dict[str, set[str]] = defaultdict(set)
    representations_by_pxd: dict[str, set[str]] = defaultdict(set)
    run_ids_by_pxd: dict[str, set[str]] = defaultdict(set)
    stems_by_pxd: dict[str, set[str]] = defaultdict(set)
    for row in selected:
        pxd = row["pxd_accession"].upper()
        task_id = row["task_id"].strip()
        name = row["archive_file"]
        if not ACCESSION.fullmatch(pxd):
            raise ValueError(f"Invalid accession: {pxd}")
        if (
            not task_id.isascii()
            or not task_id.isdecimal()
            or int(task_id) < 1
            or task_id in task_ids
        ):
            raise ValueError(f"Duplicate/invalid task_id: {task_id}")
        task_ids.add(task_id)
        if not FILENAME.fullmatch(name) or name in {".", ".."}:
            raise ValueError(f"Unsafe filename: {name!r}; v1 requires simple basenames")
        if not name.lower().endswith(SUPPORTED):
            raise ValueError(f"Unsupported input format: {name}; do not silently skip files")
        if name.casefold() in names_by_pxd[pxd]:
            raise ValueError(f"Duplicate filename for {pxd}: {name}")
        names_by_pxd[pxd].add(name.casefold())
        lower = name.casefold()
        suffix = next((s for s in (".mzml.gz", ".raw.gz", ".mzml", ".raw") if lower.endswith(s)), None)
        if suffix:
            kind = "mzml" if suffix.startswith(".mzml") else "native_raw"
            representations_by_pxd[pxd].add(kind)
            stem = lower[:-len(suffix)]
            if stem in stems_by_pxd[pxd]:
                raise ValueError(f"Duplicate or ambiguous RAW/mzML representations: {pxd}/{stem}")
            stems_by_pxd[pxd].add(stem)
        # v1 is intentionally metadata-only. SDRF mappings need a separate, tested lane.
        if (
            any(
                row.get(key, "").strip() for key in ("sdrf_path", "sdrf_basename", "sdrf_data_file")
            )
            or row.get("needs_file_map", "0").strip() != "0"
        ):
            raise ValueError(f"Task {task_id}: SDRF mapping not supported in v1 Nextflow pilot")
        by_accession[pxd].append(row)
    # Different original vendor files and converted mzMLs in the same accession
    # can represent the same acquisition. Require provenance from a reviewed
    # representation-selection inventory before allowing a mixed-format plan.
    for pxd, kinds in representations_by_pxd.items():
        if kinds == {"mzml", "native_raw"}:
            inventories = set()
            for row in by_accession[pxd]:
                run_id = row.get("logical_run_id", "").strip()
                proof = row.get("representation_inventory_sha256", "").strip()
                if (row.get("representation_policy", "") != "curated_single_primary_v1"
                        or not MIXED_RUN_ID.fullmatch(run_id)
                        or not HEX_SHA256.fullmatch(proof)):
                    raise ValueError(f"Mixed vendor/mzML project needs a curated manifest: {pxd}")
                if run_id in run_ids_by_pxd[pxd]:
                    raise ValueError(f"Duplicate logical MS acquisition: {pxd}/{run_id}")
                run_ids_by_pxd[pxd].add(run_id)
                inventories.add(proof)
            if len(inventories) != 1:
                raise ValueError(f"Inconsistent representation-inventory provenance: {pxd}")
    expected: dict[str, int] = {}
    report_ids: set[str] = set()
    for row in reports:
        pxd = row["pxd_accession"].upper()
        if pxd not in by_accession:
            continue
        count = row["expected_file_tasks"].strip()
        rid = row["report_id"].strip()
        if (
            pxd in expected
            or rid in report_ids
            or not count.isascii()
            or not count.isdecimal()
            or int(count) < 1
        ):
            raise ValueError(f"Invalid/duplicate report row: {row}")
        expected[pxd] = int(count)
        report_ids.add(rid)
    if set(expected) != set(by_accession):
        raise ValueError(
            f"Missing accession report rows: {sorted(set(by_accession) - set(expected))}"
        )
    for pxd, rows in by_accession.items():
        if len(rows) != expected[pxd]:
            raise ValueError(f"Count mismatch for {pxd}: {len(rows)} vs {expected[pxd]}")
    output.mkdir(parents=True, exist_ok=True)
    tasks = output / "tasks"
    tasks.mkdir(exist_ok=True)
    if any(tasks.iterdir()):
        raise ValueError(f"Output task directory is not empty: {tasks}")
    for pxd in sorted(by_accession):
        for row in sorted(by_accession[pxd], key=lambda x: int(x["task_id"])):
            entry = {
                "schema_version": 1,
                "task_id": int(row["task_id"]),
                "accession": pxd,
                "archive_file": row["archive_file"],
                "expected_runs": expected[pxd],  # number of source archive tasks, not logical MS runs
                "source_manifest": file_manifest.name,
                "logical_run_id": row.get("logical_run_id", "").strip(),
                "representation_policy": row.get("representation_policy", "").strip(),
                "representation_inventory_sha256": row.get("representation_inventory_sha256", "").strip(),
            }
            (tasks / f"{pxd}__{int(row['task_id']):09d}.json").write_text(
                json.dumps(entry, sort_keys=True, indent=2) + "\n", encoding="utf-8"
            )
    summary = {
        "schema_version": 1,
        "expected_by_accession": expected,
        "tasks": sum(expected.values()),
        "selected_accessions": sorted(expected),
    }
    (output / "plan.json").write_text(json.dumps(summary, sort_keys=True, indent=2) + "\n")
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--report-manifest", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--accession", action="append", default=[])
    args = parser.parse_args()
    print(
        json.dumps(
            prepare(args.manifest, args.report_manifest, args.outdir, set(args.accession)),
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
