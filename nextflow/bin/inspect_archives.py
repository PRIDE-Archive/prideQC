#!/usr/bin/env python3
"""Inspect the ZIP archives of PRIDE projects without downloading them. Read-only, public API.

For each archive file of the inventory (reason SKIP_ARCHIVE_CONTENTS_NOT_INSPECTED), read the
ZIP central directory with HTTP range reads (a few KB at the end of the file) and list the
members. The same member rules as the compute node apply (archive_members.py).

Outputs (new directory, refuses an existing one):
  archives.tsv   one row per archive: runs, run bytes, largest run, unsupported members, decision
  members.tsv    one row per run or unsupported member: the link archive -> member file, for SDRF
  queue.tsv      QUEUE_ARCHIVE_RUNS rows for production_loop.py (memory from the largest run)
  summary.json   counts and bytes
A cache directory (one JSON per archive) makes the inspection resumable.

Decisions:
  QUEUE_ARCHIVE_RUNS               at least one .raw/.mzML/.mzML.gz run member
  SKIP_ARCHIVE_NO_RUNS             no run and no unsupported member (results, peak lists, text):
                                   not raw data, it does not block a complete project
  SKIP_ARCHIVE_UNSUPPORTED_RUNS    only runs without a reader (.d or Waters .raw folder, .wiff,
                                   .raw.gz, nested archive)
  SKIP_ARCHIVE_FORMAT_NOT_INSPECTED  7z, rar, tar: the contents are not inspected (yet)
  SKIP_ARCHIVE_UNREADABLE          the listing failed (not a ZIP, no range reads, network error)
Python 3.9 compatible (Codon host python3).
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import re
import sys
import threading
import time
import urllib.error
from collections.abc import Callable
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import archive_members as AM  # noqa: E402

API = "https://www.ebi.ac.uk/pride/ws/archive/v3"
PAGE = 100
ZIP_SUFFIXES = (".zip",)
OTHER_ARCHIVES = (".7z", ".rar", ".tar", ".tar.gz", ".tgz")
MEMORY_TIERS = ((1_000_000_000, 16), (3_000_000_000, 32), (6_000_000_000, 64))
ARCHIVE_FIELDS = ["accession", "file_name", "file_size_bytes", "status", "runs", "run_bytes",
                  "largest_run_bytes", "unsupported", "ignored", "mzml_skipped_for_raw", "decision", "reason"]
MEMBER_FIELDS = ["accession", "archive_file", "member", "member_file", "member_size_bytes", "member_format",
                 "kind"]
QUEUE_FIELDS = ["accession", "file_name", "file_category", "file_format", "file_size_bytes", "checksum",
                "memory_gb", "decision", "runs", "run_bytes", "largest_run_bytes"]

Fetcher = Callable[[str], Any]


def memory_tier_gb(size: int) -> int:
    for bound, gb in MEMORY_TIERS:
        if size < bound:
            return gb
    return 128


def http_json(url: str, *, retries: int = 4, timeout: int = 60) -> Any:
    import urllib.request
    delay = 2.0
    for attempt in range(retries + 1):
        try:
            req = urllib.request.Request(url, headers={"Accept": "application/json",
                                                       "User-Agent": "prideqc-archive-inspect/1"})
            with urllib.request.urlopen(req, timeout=timeout) as response:
                return json.loads(response.read().decode("utf-8"))
        except (urllib.error.URLError, TimeoutError, ValueError, ConnectionError) as err:
            if attempt == retries or getattr(err, "code", None) in (400, 404):
                raise
            time.sleep(delay)
            delay *= 2
    raise RuntimeError("unreachable")


def public_urls(fetch: Fetcher, accession: str) -> dict:
    """file name -> https URL of the PRIDE FTP copy (from the v3 API file list)."""
    out: dict = {}
    page = 0
    while True:
        batch = fetch(f"{API}/projects/{accession}/files?pageSize={PAGE}&page={page}")
        if not isinstance(batch, list) or not batch:
            break
        for rec in batch:
            for loc in rec.get("publicFileLocations") or []:
                value = str(loc.get("value") or "")
                if value.startswith("ftp://ftp.pride.ebi.ac.uk/"):
                    out[str(rec.get("fileName"))] = "https://" + value[len("ftp://"):]
        if len(batch) < PAGE:
            break
        page += 1
    return out


def inspect_one(accession: str, name: str, size: int, url: str | None, lister: Callable = AM.list_zip) -> dict:
    """Archive row (and its member rows under '_members'); never raises."""
    row: dict[str, Any] = {"accession": accession, "file_name": name, "file_size_bytes": size, "status": "",
                           "runs": 0, "run_bytes": 0, "largest_run_bytes": 0, "unsupported": 0, "ignored": 0,
                           "mzml_skipped_for_raw": 0, "decision": "", "reason": "", "_members": []}
    lower = name.casefold()
    if lower.endswith(OTHER_ARCHIVES) or not lower.endswith(ZIP_SUFFIXES):
        row.update(status="not_inspected", decision="SKIP_ARCHIVE_FORMAT_NOT_INSPECTED",
                   reason=next((s.lstrip(".") for s in OTHER_ARCHIVES if lower.endswith(s)), "other"))
        return row
    if not url:
        row.update(status="no_public_url", decision="SKIP_ARCHIVE_UNREADABLE", reason="no_public_url")
        return row
    try:
        plan = lister(url)
    except AM.ArchiveMemberError as err:
        row.update(status="unreadable", decision="SKIP_ARCHIVE_UNREADABLE", reason=_reason(str(err)))
        return row
    except (OSError, urllib.error.URLError, ValueError, EOFError) as err:
        detail = str(getattr(err, "reason", "") or getattr(err, "code", "") or "")
        detail = re.sub(r"[^A-Za-z0-9_ .:-]", "", detail)[:60]
        row.update(status="unreadable", decision="SKIP_ARCHIVE_UNREADABLE",
                   reason=type(err).__name__ + (f":{detail}" if detail else ""))
        return row
    runs = plan["runs"]
    row.update(status="inspected", runs=len(runs), run_bytes=sum(r["size"] for r in runs),
               largest_run_bytes=max((r["size"] for r in runs), default=0),
               unsupported=len(plan["unsupported"]), ignored=plan["ignored"],
               mzml_skipped_for_raw=plan["skipped_mzml_with_raw"])
    if runs:
        row.update(decision="QUEUE_ARCHIVE_RUNS",
                   reason="RUNS_AND_UNSUPPORTED_MEMBERS" if plan["unsupported"] else "RUN_MEMBERS")
    elif plan["unsupported"]:
        row.update(decision="SKIP_ARCHIVE_UNSUPPORTED_RUNS", reason="UNSUPPORTED_MEMBERS_ONLY")
    else:
        row.update(decision="SKIP_ARCHIVE_NO_RUNS", reason="NO_RAW_DATA")
    members = [{"accession": accession, "archive_file": name, "member": r["member"], "member_file": r["basename"],
                "member_size_bytes": r["size"], "member_format": r["format"], "kind": "run"} for r in runs]
    members += [{"accession": accession, "archive_file": name, "member": m, "member_file": m.rstrip("/").split("/")[-1],
                 "member_size_bytes": "", "member_format": "", "kind": "unsupported"} for m in plan["unsupported"]]
    row["_members"] = members
    return row


def _reason(text: str) -> str:
    for key, label in (("range", "no_range_reads"), ("Not a valid ZIP", "not_a_zip"), ("Duplicate", "duplicate_member"),
                       ("Unsafe", "unsafe_member_name"), ("too many", "too_many_members"), ("size", "bad_member_size")):
        if key in text:
            return label
    return "zip_error"


def read_archives(inventory: Path, accessions: set) -> dict:
    """accession -> [(file name, size)] of the archive files that the inventory did not inspect."""
    out: dict = {}
    with inventory.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["reason"].split(";")[0] != "SKIP_ARCHIVE_CONTENTS_NOT_INSPECTED":
                continue
            if accessions and row["accession"] not in accessions:
                continue
            size = int(row["file_size_bytes"]) if row["file_size_bytes"].isdigit() else 0
            out.setdefault(row["accession"], []).append((row["file_name"], size, row.get("file_category", "")))
    return out


def _cache_file(cache: Path, accession: str, name: str) -> Path:
    import hashlib
    return cache / accession / (hashlib.sha256(name.encode()).hexdigest()[:24] + ".json")


def _cached(cache: Path | None, accession: str, name: str, size: int) -> dict | None:
    target = _cache_file(cache, accession, name) if cache else None
    if target is None or not target.is_file():
        return None
    try:
        cached = json.loads(target.read_text(encoding="utf-8"))
    except ValueError:
        return None
    return cached if cached.get("file_size_bytes") == size and cached.get("decision") else None


def _store(cache: Path | None, row: dict) -> None:
    # only final answers are cached; an unreadable archive is listed again on the next run
    if not cache or row["status"] not in ("inspected", "not_inspected"):
        return
    target = _cache_file(cache, row["accession"], row["file_name"])
    target.parent.mkdir(parents=True, exist_ok=True)
    tmp = target.with_suffix(f".json.partial-{os.getpid()}-{threading.get_ident()}")
    tmp.write_text(json.dumps({k: v for k, v in row.items() if k != "_category"}, sort_keys=True), encoding="utf-8")
    tmp.replace(target)


def project_urls(accession: str, files: list, fetch: Fetcher) -> dict:
    if not any(n.casefold().endswith(ZIP_SUFFIXES) for n, _, _ in files):
        return {}
    try:
        return public_urls(fetch, accession)
    except (OSError, urllib.error.URLError, ValueError):
        return {}


def write_tsv(path: Path, fields: list, rows: list) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def build(inventory: Path, outdir: Path, *, cache: Path | None = None, workers: int = 16,
          accessions: set | None = None, limit: int | None = None, fetch: Fetcher = http_json,
          lister: Callable = AM.list_zip, progress: Callable | None = None) -> dict:
    if outdir.exists() or outdir.is_symlink():
        raise ValueError(f"Refusing existing output directory: {outdir}")
    projects = read_archives(inventory, accessions or set())
    items = sorted(projects.items())
    if limit is not None:
        items = items[:limit]
    rows: list = []
    todo: list = []
    for acc, files in items:
        for name, size, category in files:
            hit = _cached(cache, acc, name, size)
            if hit is not None:
                rows.append(dict(hit, _category=category))
            else:
                todo.append((acc, name, size, category))
    pending: dict = {}
    for acc, name, size, category in todo:
        pending.setdefault(acc, []).append((name, size, category))
    # two stages, both parallel over all archives: the file URLs of each project, then one
    # listing per archive (a project with 1,000 archives does not run on one thread)
    with ThreadPoolExecutor(max_workers=workers) as pool:
        url_jobs = {acc: pool.submit(project_urls, acc, files, fetch) for acc, files in pending.items()}
        urls = {acc: job.result() for acc, job in url_jobs.items()}

        def one(job: tuple) -> dict:
            acc, name, size, category = job
            row = inspect_one(acc, name, size, urls[acc].get(name), lister)
            row["_category"] = category
            _store(cache, row)
            return row

        for i, row in enumerate(pool.map(one, todo), 1):
            rows.append(row)
            if progress:
                progress(i, len(todo))
    rows.sort(key=lambda r: (r["accession"], r["file_name"].casefold()))
    stage = outdir.with_name(f".{outdir.name}.partial-{os.getpid()}")
    stage.mkdir(parents=True)
    write_tsv(stage / "archives.tsv", ARCHIVE_FIELDS, rows)
    write_tsv(stage / "members.tsv", MEMBER_FIELDS, [m for r in rows for m in r["_members"]])
    queue = [{"accession": r["accession"], "file_name": r["file_name"], "file_category": r.get("_category", ""),
              "file_format": "zip_archive", "file_size_bytes": r["file_size_bytes"], "checksum": "",
              "memory_gb": memory_tier_gb(int(r["largest_run_bytes"])), "decision": "QUEUE_ARCHIVE_RUNS",
              "runs": r["runs"], "run_bytes": r["run_bytes"], "largest_run_bytes": r["largest_run_bytes"]}
             for r in rows if r["decision"] == "QUEUE_ARCHIVE_RUNS"]
    write_tsv(stage / "queue.tsv", QUEUE_FIELDS, queue)

    def count(key: str) -> dict:
        out: dict = {}
        for r in rows:
            out[r[key]] = out.get(r[key], 0) + 1
        return dict(sorted(out.items()))

    summary = {"schema_version": 1, "inventory": str(inventory), "archives": len(rows), "projects": len(items),
               "by_decision": count("decision"), "by_reason": count("reason"),
               "queue_archives": len(queue), "queue_projects": len({r["accession"] for r in queue}),
               "queue_runs": sum(int(r["runs"]) for r in queue),
               "queue_run_tb": round(sum(int(r["run_bytes"]) for r in queue) / 1e12, 2),
               "queue_archive_tb": round(sum(int(r["file_size_bytes"]) for r in queue) / 1e12, 2),
               "queue_memory_gb": {str(k): sum(1 for r in queue if r["memory_gb"] == k) for k in (16, 32, 64, 128)}}
    (stage / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    stage.rename(outdir)
    return summary


def main(argv: list | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--inventory", type=Path, required=True, help="pride_inventory.py inventory.tsv")
    ap.add_argument("--outdir", type=Path, required=True, help="new output directory")
    ap.add_argument("--cache", type=Path, help="resumable cache directory (one JSON per archive)")
    ap.add_argument("--workers", type=int, default=16)
    ap.add_argument("--accession", action="append", default=[], help="only these projects (repeatable)")
    ap.add_argument("--limit", type=int, help="only the first N projects (for a trial)")
    args = ap.parse_args(argv)

    def progress(i: int, n: int) -> None:
        if i % 500 == 0 or i == n:
            print(json.dumps({"event": "PROGRESS", "archives": i, "of": n}), flush=True)

    try:
        summary = build(args.inventory, args.outdir, cache=args.cache, workers=args.workers,
                        accessions=set(args.accession), limit=args.limit, progress=progress)
    except (ValueError, OSError) as err:
        print(json.dumps({"event": "BLOCKED", "reason": f"{type(err).__name__}: {err}"}), file=sys.stderr)
        return 2
    print(json.dumps(dict(summary, event="DONE"), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
