#!/usr/bin/env python3
"""Apply explicitly accepted prideQC PTM review rows to an SDRF proposal copy.

This tool is intentionally human-gated. The probability calibrator writes
``sdrf_status=review-required``. A reviewer must explicitly change a row to
``accepted`` before this script will add it to ``comment[modification parameters]``.

Only NT and AC are written. RAW mass-shift evidence cannot establish search
configuration (fixed/variable), target amino acid, position, or localization.
"""

from __future__ import annotations

import argparse
import csv
import re
from dataclasses import dataclass
from pathlib import Path

MODIFICATION_COLUMN = "comment[modification parameters]"


@dataclass(frozen=True, slots=True)
class AcceptedModification:
    pxd_accession: str
    unimod_accession: str
    unimod_name: str
    sdrf_value: str


def _accessions(value: str) -> set[str]:
    return {
        match.strip().upper()
        for match in re.findall(r"(?:^|;)\s*AC=([^;]+)", value, flags=re.IGNORECASE)
    }


def read_review(path: Path) -> list[AcceptedModification]:
    accepted: list[AcceptedModification] = []
    with path.open(newline="", encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if (row.get("sdrf_status") or "").strip().casefold() != "accepted":
                continue
            accession = (row.get("unimod_accession") or "").strip()
            name = (row.get("unimod_name") or "").strip()
            value = (row.get("sdrf_value") or "").strip()
            pxd = (row.get("pxd_accession") or "").strip()
            if not accession or not name or not value:
                raise ValueError(
                    "Accepted review rows require UniMod accession, name, and SDRF value."
                )
            if accession.upper() not in _accessions(value):
                raise ValueError(
                    f"Accepted SDRF value does not contain its UniMod accession: {value!r}"
                )
            accepted.append(AcceptedModification(pxd, accession, name, value))
    return accepted


def read_sdrf(path: Path) -> tuple[list[str], list[str], list[list[str]]]:
    preamble: list[str] = []
    with path.open(encoding="utf-8-sig", newline="") as handle:
        first = handle.readline()
        while first.startswith("#") or first in ("\n", "\r\n"):
            preamble.append(first)
            first = handle.readline()
        if not first:
            raise ValueError("SDRF is empty.")
        columns = next(csv.reader([first], delimiter="\t"))
        rows = list(csv.reader(handle, delimiter="\t"))
    if any(len(row) != len(columns) for row in rows):
        raise ValueError("SDRF contains a row with a different number of fields than its header.")
    return preamble, columns, rows


def apply_review(
    columns: list[str],
    rows: list[list[str]],
    accepted: list[AcceptedModification],
) -> tuple[list[str], list[list[str]], list[dict[str, str]]]:
    new_columns = list(columns)
    new_rows = [list(row) for row in rows]
    modification_indices = [
        index
        for index, name in enumerate(columns)
        if name.strip().casefold() == MODIFICATION_COLUMN.casefold()
    ]
    existing_accessions: set[str] = set()
    for row in rows:
        for index in modification_indices:
            existing_accessions.update(_accessions(row[index]))

    audit: list[dict[str, str]] = []
    seen: set[str] = set()
    for item in accepted:
        accession = item.unimod_accession.upper()
        if accession in seen:
            audit.append(
                {
                    "pxd_accession": item.pxd_accession,
                    "unimod_accession": item.unimod_accession,
                    "unimod_name": item.unimod_name,
                    "sdrf_value": item.sdrf_value,
                    "action": "skipped-duplicate-review",
                }
            )
            continue
        seen.add(accession)
        if accession in existing_accessions:
            audit.append(
                {
                    "pxd_accession": item.pxd_accession,
                    "unimod_accession": item.unimod_accession,
                    "unimod_name": item.unimod_name,
                    "sdrf_value": item.sdrf_value,
                    "action": "skipped-existing-accession",
                }
            )
            continue
        new_columns.append(MODIFICATION_COLUMN)
        for row in new_rows:
            row.append(item.sdrf_value)
        existing_accessions.add(accession)
        audit.append(
            {
                "pxd_accession": item.pxd_accession,
                "unimod_accession": item.unimod_accession,
                "unimod_name": item.unimod_name,
                "sdrf_value": item.sdrf_value,
                "action": "added",
            }
        )
    return new_columns, new_rows, audit


def write_sdrf(
    path: Path,
    preamble: list[str],
    columns: list[str],
    rows: list[list[str]],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        handle.writelines(preamble)
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(columns)
        writer.writerows(rows)


def write_audit(path: Path, rows: list[dict[str, str]]) -> None:
    fields = [
        "pxd_accession",
        "unimod_accession",
        "unimod_name",
        "sdrf_value",
        "action",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Apply only explicitly accepted prideQC PTM review rows to a new SDRF proposal file."
        )
    )
    parser.add_argument("--sdrf", required=True, type=Path)
    parser.add_argument("--review", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--audit", required=True, type=Path)
    args = parser.parse_args()

    source = args.sdrf.resolve()
    output = args.output.resolve()
    if source == output:
        raise SystemExit("Refusing to overwrite the source SDRF; choose a distinct --output path.")

    accepted = read_review(args.review)
    if not accepted:
        raise SystemExit("No sdrf_status=accepted rows were found in the review TSV.")
    preamble, columns, rows = read_sdrf(args.sdrf)
    proposed_columns, proposed_rows, audit = apply_review(columns, rows, accepted)
    write_sdrf(args.output, preamble, proposed_columns, proposed_rows)
    write_audit(args.audit, audit)

    added = sum(item["action"] == "added" for item in audit)
    print(f"accepted_review_rows={len(accepted)}")
    print(f"added_modification_columns={added}")
    print(f"output={args.output.resolve()}")
    print(f"audit={args.audit.resolve()}")
    print("validation_required=1")


if __name__ == "__main__":
    main()
