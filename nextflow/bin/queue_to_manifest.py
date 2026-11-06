#!/usr/bin/env python3
"""Cut one bounded batch from a pride_inventory.py queue.tsv into plan_manifest inputs.

Writes a NEW directory (refuses an existing one):
  file-manifest.tsv    task_id, pxd_accession, archive_file, sdrf_path, needs_file_map, file_size_bytes
  report-manifest.tsv  report_id, pxd_accession, expected_file_tasks (= files in THIS batch)
  batch.json           selection rules, counts, bytes, memory tiers, input SHA256

Rules: only QUEUE_NATIVE_RAW rows (the mzML lane is not qualified); rows already
canonical or failed in an export index are excluded; one batch never claims a
complete study. Selection is deterministic: the largest files first up to
--large-files, then a round-robin over accessions (smallest first) up to the limits.
Python 3.9 compatible.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path

DECISION = "QUEUE_NATIVE_RAW"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def excluded(index: Path | None) -> set[tuple[str, str]]:
    keys: set[tuple[str, str]] = set()
    if index:
        for line in index.read_text(encoding="utf-8").splitlines():
            row = json.loads(line)
            if row.get("role") in ("canonical", "failed", "superseded_failure"):
                keys.add((str(row.get("accession")), str(row.get("source_file", "")).casefold()))
    return keys


def select(rows: list[dict], *, max_files: int, max_bytes: int, large_files: int,
           per_accession: int, skip: set[tuple[str, str]]) -> list[dict]:
    pool = [r for r in rows if r.get("decision") == DECISION
            and (r["accession"], r["file_name"].casefold()) not in skip
            and str(r.get("file_size_bytes", "")).isdigit() and int(r["file_size_bytes"]) > 0]
    for r in pool:
        r["_size"] = int(r["file_size_bytes"])
    picked: list[dict] = []
    seen: set[tuple[str, str]] = set()
    used = 0
    count: dict[str, int] = defaultdict(int)

    def take(r: dict) -> bool:
        nonlocal used
        key = (r["accession"], r["file_name"].casefold())
        if key in seen or len(picked) >= max_files or used + r["_size"] > max_bytes:
            return False
        if count[r["accession"]] >= per_accession:
            return False
        picked.append(r)
        seen.add(key)
        used += r["_size"]
        count[r["accession"]] += 1
        return True

    for r in sorted(pool, key=lambda x: (-x["_size"], x["accession"], x["file_name"].casefold()))[:large_files]:
        take(r)
    by_acc: dict[str, list[dict]] = defaultdict(list)
    for r in sorted(pool, key=lambda x: (x["_size"], x["file_name"].casefold())):
        by_acc[r["accession"]].append(r)
    progressed = True
    while progressed and len(picked) < max_files:
        progressed = False
        for acc in sorted(by_acc):
            while by_acc[acc]:
                r = by_acc[acc].pop(0)
                if (r["accession"], r["file_name"].casefold()) in seen:
                    continue
                if take(r):
                    progressed = True
                break
    return sorted(picked, key=lambda r: (r["accession"], r["file_name"].casefold()))


def write_tsv(path: Path, fields: list[str], rows: list[dict]) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def build(queue: Path, outdir: Path, *, max_files: int, max_gb: float, large_files: int = 0,
          per_accession: int = 1_000_000, processed_index: Path | None = None) -> dict:
    if outdir.exists() or outdir.is_symlink():
        raise ValueError(f"Refusing existing output directory: {outdir}")
    with queue.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    picked = select(rows, max_files=max_files, max_bytes=int(max_gb * 1e9), large_files=large_files,
                    per_accession=per_accession, skip=excluded(processed_index))
    if not picked:
        raise ValueError("Empty batch")
    stage = outdir.with_name(f".{outdir.name}.partial")
    stage.mkdir(parents=True)
    tasks = [{"task_id": str(i), "pxd_accession": r["accession"], "archive_file": r["file_name"],
              "sdrf_path": "", "needs_file_map": "0", "file_size_bytes": r["file_size_bytes"]}
             for i, r in enumerate(picked, 1)]
    per: dict[str, int] = defaultdict(int)
    for r in picked:
        per[r["accession"]] += 1
    reports = [{"report_id": str(i), "pxd_accession": a, "expected_file_tasks": str(n)}
               for i, (a, n) in enumerate(sorted(per.items()), 1)]
    write_tsv(stage / "file-manifest.tsv", ["task_id", "pxd_accession", "archive_file", "sdrf_path",
                                            "needs_file_map", "file_size_bytes"], tasks)
    write_tsv(stage / "report-manifest.tsv", ["report_id", "pxd_accession", "expected_file_tasks"], reports)
    tiers: dict[str, int] = defaultdict(int)
    for r in picked:
        tiers[str(r.get("memory_gb", ""))] += 1
    summary = {"schema_version": 1, "queue": str(queue), "queue_sha256": sha256(queue),
               "processed_index": str(processed_index) if processed_index else None,
               "rules": {"decision": DECISION, "max_files": max_files, "max_gb": max_gb,
                         "large_files": large_files, "per_accession": per_accession},
               "files": len(picked), "accessions": len(per),
               "total_gb": round(sum(r["_size"] for r in picked) / 1e9, 2),
               "largest_gb": round(max(r["_size"] for r in picked) / 1e9, 2),
               "memory_gb_tiers": dict(sorted(tiers.items())),
               "complete_study_asserted": False,
               "file_manifest_sha256": sha256(stage / "file-manifest.tsv"),
               "report_manifest_sha256": sha256(stage / "report-manifest.tsv")}
    (stage / "batch.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    stage.rename(outdir)
    return summary


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--queue", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--max-files", type=int, required=True)
    ap.add_argument("--max-gb", type=float, required=True)
    ap.add_argument("--large-files", type=int, default=0, help="Take the N largest files first")
    ap.add_argument("--per-accession", type=int, default=1_000_000)
    ap.add_argument("--processed-index", type=Path, help="export_run_index.py index.jsonl to exclude")
    args = ap.parse_args()
    try:
        print(json.dumps(build(args.queue, args.outdir, max_files=args.max_files, max_gb=args.max_gb,
                               large_files=args.large_files, per_accession=args.per_accession,
                               processed_index=args.processed_index), indent=2, sort_keys=True))
        return 0
    except (ValueError, OSError) as err:
        print(json.dumps({"batch": "BLOCKED", "reason": f"{type(err).__name__}: {err}"}))
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
