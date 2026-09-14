#!/usr/bin/env python3
"""Resolve curated ground-truth RAW files to task IDs in an existing manifest.

The generated TSV always uses Unix LF line endings so its final ``task_id``
column is safe to consume directly from POSIX shell tools and Slurm array
specifications.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

OUTPUT_FIELDS = ["gt_id", "pxd_accession", "raw_file", "task_id"]
REQUIRED_GT_FIELDS = set(OUTPUT_FIELDS[:3])
REQUIRED_MANIFEST_FIELDS = {"task_id", "pxd_accession", "archive_file"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gt", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def _read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = list(reader.fieldnames or [])
        rows = [dict(row) for row in reader]
    return fieldnames, rows


def resolve_tasks(gt_path: Path, manifest_path: Path) -> list[dict[str, str]]:
    gt_fields, gt_rows = _read_tsv(gt_path)
    manifest_fields, manifest_rows = _read_tsv(manifest_path)

    missing_gt = REQUIRED_GT_FIELDS.difference(gt_fields)
    if missing_gt:
        raise ValueError(
            "Ground-truth TSV is missing required columns: "
            + ", ".join(sorted(missing_gt))
        )

    missing_manifest = REQUIRED_MANIFEST_FIELDS.difference(manifest_fields)
    if missing_manifest:
        raise ValueError(
            "Manifest TSV is missing required columns: "
            + ", ".join(sorted(missing_manifest))
        )

    index: dict[tuple[str, str], list[dict[str, str]]] = {}
    for row in manifest_rows:
        key = (row["pxd_accession"].strip(), row["archive_file"].strip())
        index.setdefault(key, []).append(row)

    resolved: list[dict[str, str]] = []
    for gt in gt_rows:
        gt_id = gt["gt_id"].strip()
        pxd = gt["pxd_accession"].strip()
        raw_file = gt["raw_file"].strip()
        matches = index.get((pxd, raw_file), [])
        if len(matches) != 1:
            raise ValueError(
                f"{gt_id}: expected exactly one manifest match for "
                f"({pxd!r}, {raw_file!r}), found {len(matches)}"
            )
        resolved.append(
            {
                "gt_id": gt_id,
                "pxd_accession": pxd,
                "raw_file": raw_file,
                "task_id": matches[0]["task_id"].strip(),
            }
        )

    task_ids = [row["task_id"] for row in resolved]
    if len(task_ids) != len(set(task_ids)):
        raise ValueError(f"Ground-truth rows resolve to duplicate physical task IDs: {task_ids}")

    return resolved


def write_task_subset(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=OUTPUT_FIELDS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def slurm_array(rows: list[dict[str, str]]) -> str:
    return ",".join(row["task_id"] for row in rows)


def main() -> int:
    args = parse_args()
    rows = resolve_tasks(args.gt, args.manifest)
    write_task_subset(args.output, rows)
    print(f"resolved={len(rows)}")
    print(f"output={args.output}")
    print(f"task_array={slurm_array(rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
