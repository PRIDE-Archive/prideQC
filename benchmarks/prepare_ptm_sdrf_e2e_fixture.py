#!/usr/bin/env python3
"""Prepare a controlled SDRF applicator E2E fixture from a real SDRF.

This is an evaluation helper only. It:
1. selects exactly one v4 review row by PXD + UniMod accession;
2. marks only that row accepted;
3. removes any existing repeated modification-parameter column containing that
   accession from a COPY of the source SDRF;
4. writes hashes and metadata.

The source SDRF is never modified.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from pathlib import Path

MOD = "comment[modification parameters]"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def accessions(value: str) -> set[str]:
    return {
        x.strip().upper()
        for x in re.findall(r"(?:^|;)\s*AC=([^;]+)", value, flags=re.IGNORECASE)
    }


def read_sdrf(path: Path):
    preamble = []
    with path.open(encoding="utf-8-sig", newline="") as handle:
        first = handle.readline()
        while first.startswith("#") or first in ("\n", "\r\n"):
            preamble.append(first)
            first = handle.readline()
        if not first:
            raise SystemExit("source SDRF is empty")
        columns = next(csv.reader([first], delimiter="\t"))
        rows = list(csv.reader(handle, delimiter="\t"))
    if any(len(row) != len(columns) for row in rows):
        raise SystemExit("source SDRF contains non-rectangular rows")
    return preamble, columns, rows


def write_sdrf(path: Path, preamble, columns, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        handle.writelines(preamble)
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(columns)
        writer.writerows(rows)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--source-sdrf", required=True, type=Path)
    p.add_argument("--review", required=True, type=Path)
    p.add_argument("--pxd-accession", required=True)
    p.add_argument("--unimod-accession", required=True)
    p.add_argument("--output-dir", required=True, type=Path)
    args = p.parse_args()

    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    target_pxd = args.pxd_accession.strip().upper()
    target_unimod = args.unimod_accession.strip().upper()

    with args.review.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        review_rows = list(reader)
    matches = [
        row for row in review_rows
        if (row.get("pxd_accession") or "").strip().upper() == target_pxd
        and (row.get("unimod_accession") or "").strip().upper() == target_unimod
    ]
    if len(matches) != 1:
        raise SystemExit(
            f"expected exactly one review row for {target_pxd} {target_unimod}; "
            f"found {len(matches)}"
        )
    selected = dict(matches[0])
    if (selected.get("sdrf_status") or "").strip() != "review-required":
        raise SystemExit(
            "selected v4 row is not review-required: "
            + (selected.get("sdrf_status") or "")
        )
    selected["sdrf_status"] = "accepted"

    accepted_path = out / "accepted-review.tsv"
    with accepted_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerow(selected)

    preamble, columns, rows = read_sdrf(args.source_sdrf)
    remove_indices = []
    for index, column in enumerate(columns):
        if column.strip().casefold() != MOD.casefold():
            continue
        if any(target_unimod in accessions(row[index]) for row in rows):
            remove_indices.append(index)
    if not remove_indices:
        raise SystemExit(
            f"source SDRF has no repeated modification column containing {target_unimod}; "
            "controlled removal fixture cannot be constructed"
        )

    keep = [i for i in range(len(columns)) if i not in set(remove_indices)]
    base_columns = [columns[i] for i in keep]
    base_rows = [[row[i] for i in keep] for row in rows]
    base_path = out / "base-without-target.sdrf.tsv"
    write_sdrf(base_path, preamble, base_columns, base_rows)

    metadata = {
        "purpose": "controlled-real-SDRF-applicator-mechanics-test",
        "source_sdrf": str(args.source_sdrf.resolve()),
        "source_sdrf_sha256": sha256(args.source_sdrf.resolve()),
        "pxd_accession": target_pxd,
        "unimod_accession": target_unimod,
        "removed_modification_columns": len(remove_indices),
        "removed_column_indices_zero_based": remove_indices,
        "base_sdrf": str(base_path),
        "base_sdrf_sha256": sha256(base_path),
        "accepted_review": str(accepted_path),
        "accepted_review_sha256": sha256(accepted_path),
        "warning": (
            "Evaluation fixture only: the target modification was deliberately removed "
            "from a copy of a trusted SDRF to test guarded re-addition."
        ),
    }
    metadata_path = out / "fixture-metadata.json"
    metadata_path.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    print(f"removed_modification_columns={len(remove_indices)}")
    print(f"base_sdrf={base_path}")
    print(f"accepted_review={accepted_path}")
    print(f"metadata={metadata_path}")


if __name__ == "__main__":
    main()
