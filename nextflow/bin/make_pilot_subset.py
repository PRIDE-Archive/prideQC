#!/usr/bin/env python3
"""Derive an explicitly partial accession pilot manifest; never overwrite sources."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--file-manifest", type=Path, required=True)
    p.add_argument("--accession", required=True)
    p.add_argument("--take", type=int, default=2)
    p.add_argument("--outdir", type=Path, required=True)
    a = p.parse_args()
    if a.take < 1:
        raise ValueError("--take must be positive")
    if a.outdir.exists() and any(a.outdir.iterdir()):
        raise ValueError("Refusing nonempty output directory")
    with a.file_manifest.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames
        if not fields or not {"task_id", "pxd_accession", "archive_file"} <= set(fields):
            raise ValueError("Invalid original file manifest headers")
        rows = [row for row in reader if row["pxd_accession"].upper() == a.accession.upper()]
    rows = sorted(rows, key=lambda r: int(r["task_id"]))[: a.take]
    if len(rows) != a.take:
        raise ValueError("Not enough selected input rows")
    a.outdir.mkdir(parents=True, exist_ok=True)
    f = a.outdir / "pilot.tasks.tsv"
    with f.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    (a.outdir / "pilot.reports.tsv").write_text(
        f"report_id\tpxd_accession\texpected_file_tasks\n1\t{a.accession.upper()}\t{a.take}\n",
        encoding="utf-8",
    )
    (a.outdir / "PILOT_SCOPE.txt").write_text(
        f"PILOT_ONLY accession={a.accession.upper()} selected_tasks={a.take}; NOT A FULL-ACCESSION COMPLETENESS CLAIM\n"
    )
    print(f"PILOT_ONLY={a.accession.upper()} selected_tasks={a.take} output={a.outdir}")


if __name__ == "__main__":
    main()
