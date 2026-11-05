#!/usr/bin/env python3
"""Build a PRIDE-wide source-file inventory and QC work queue. Read-only, public API.

For each PRIDE project (PRIDE Archive REST v3), list every file, verify the
listing against the project file count, classify each file by format, select one
representation per logical run, and decide one of:

  QUEUE_NATIVE_RAW        Thermo .raw, qualified reader: send to the QC pipeline.
  QUEUE_MZML_NO_VENDOR    .mzML with no vendor file of the same stem; reader
                          qualification pending, so a separate lane.
  DONE                    already canonical in a run index (export index.jsonl).
  FAILED_PREVIOUSLY       a failure receipt exists and no success; do not retry blindly.
  SKIP_<reason>           not queued; the reason is explicit.

Stem matching between vendor and mzML files is a heuristic. Rows carry
`representation_evidence=stem_match_unverified`. It never asserts that a study is complete.

Outputs (new directory, refuses an existing one):
  inventory.tsv       one row per file with format, decision, reason
  queue.tsv           QUEUE_* rows only, in plan_manifest-compatible columns
  summary.json        counts, bytes, memory tiers, listing problems
  projects.tsv        per-project listing status (count verified or not)
A cache directory (one JSON per project) makes the crawl resumable.
Python 3.9 compatible.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import re
import sys
import time
import urllib.error
import urllib.request
from collections.abc import Callable
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any

API = "https://www.ebi.ac.uk/pride/ws/archive/v3"
PAGE = 100
PROJECT_PAGE = 1000
ACCESSION = re.compile(r"PXD\d{6,}\Z")
SAFE_NAME = re.compile(r"[A-Za-z0-9_.+\-]{1,255}\Z")
MEMORY_TIERS = ((1_000_000_000, 16), (3_000_000_000, 32), (6_000_000_000, 64))
# Longest suffix first. (suffix, format)
FORMATS = (
    (".wiff.scan", "sciex_wiff_scan"), (".wiff2", "sciex_wiff"), (".wiff", "sciex_wiff"),
    (".mzml.gz", "mzml"), (".mzml", "mzml"), (".mzxml.gz", "mzxml"), (".mzxml", "mzxml"),
    (".raw.zip", "archive"), (".d.zip", "bruker_d_archive"), (".d.tar.gz", "bruker_d_archive"),
    (".d.tar", "bruker_d_archive"), (".tar.gz", "archive"), (".tgz", "archive"),
    (".zip", "archive"), (".7z", "archive"), (".rar", "archive"), (".tar", "archive"),
    (".raw.gz", "thermo_raw_gz"), (".raw", "thermo_raw"), (".mgf.gz", "mgf"), (".mgf", "mgf"),
    (".baf", "bruker_baf"), (".tdf", "bruker_tdf"), (".yep", "bruker_yep"),
    (".lcd", "shimadzu_lcd"), (".dat", "dat"), (".mzid.gz", "result"), (".mzid", "result"),
    (".mztab", "result"), (".fasta", "fasta"), (".xml.gz", "result"), (".xml", "other"),
    (".txt", "other"), (".tsv", "other"), (".csv", "other"), (".pdf", "other"),
    (".xlsx", "other"), (".xls", "other"), (".msf", "result"), (".pdresult", "result"),
    (".mzqc", "other"), (".sdrf.tsv", "other"),
)
SPECTRUM_FORMATS = {"thermo_raw", "thermo_raw_gz", "sciex_wiff", "sciex_wiff_scan", "mzml", "mzxml",
                    "bruker_d_archive", "bruker_baf", "bruker_tdf", "bruker_yep", "shimadzu_lcd",
                    "archive", "mgf"}
VENDOR_FORMATS = {"thermo_raw", "thermo_raw_gz", "sciex_wiff", "bruker_d_archive", "bruker_baf",
                  "bruker_tdf", "bruker_yep", "shimadzu_lcd"}
SKIP_REASON = {
    "sciex_wiff": "SKIP_SCIEX_WIFF_NO_QUALIFIED_READER",
    "sciex_wiff_scan": "SKIP_SCIEX_WIFF_SCAN_COMPANION",
    "mzxml": "SKIP_MZXML_READER_NOT_QUALIFIED",
    "archive": "SKIP_ARCHIVE_CONTENTS_NOT_INSPECTED",
    "bruker_d_archive": "SKIP_BRUKER_D_NO_QUALIFIED_READER",
    "bruker_baf": "SKIP_BRUKER_NO_QUALIFIED_READER",
    "bruker_tdf": "SKIP_BRUKER_NO_QUALIFIED_READER",
    "bruker_yep": "SKIP_BRUKER_NO_QUALIFIED_READER",
    "shimadzu_lcd": "SKIP_SHIMADZU_NO_QUALIFIED_READER",
    "thermo_raw_gz": "SKIP_COMPRESSED_RAW_NOT_QUALIFIED",
    "mgf": "SKIP_PEAK_LIST_NOT_A_RUN_REPRESENTATION",
}
FIELDS = ["accession", "file_name", "file_category", "file_format", "file_size_bytes", "checksum",
          "run_stem", "decision", "reason", "memory_gb", "representation_evidence"]

Fetcher = Callable[[str], Any]


def http_fetch(url: str, *, retries: int = 4, timeout: int = 60) -> Any:
    delay = 2.0
    for attempt in range(retries + 1):
        try:
            request = urllib.request.Request(url, headers={"Accept": "application/json",
                                                           "User-Agent": "prideqc-inventory/1"})
            with urllib.request.urlopen(request, timeout=timeout) as response:
                return json.loads(response.read().decode("utf-8"))
        except (urllib.error.URLError, TimeoutError, ValueError, ConnectionError) as err:
            status = getattr(err, "code", None)
            if attempt == retries or status in (400, 404):
                raise
            time.sleep(delay)
            delay *= 2
    raise RuntimeError("unreachable")


def file_format(name: str) -> str:
    lower = name.casefold()
    for suffix, kind in FORMATS:
        if lower.endswith(suffix):
            return kind
    return "other"


def run_stem(name: str) -> str:
    lower = name.casefold()
    for suffix, _ in FORMATS:
        if lower.endswith(suffix):
            return lower[: -len(suffix)]
    return lower


def memory_tier_gb(size: int | None) -> int | None:
    if size is None:
        return None
    for bound, gb in MEMORY_TIERS:
        if size < bound:
            return gb
    return 128


def list_projects(fetch: Fetcher, limit: int | None = None) -> list[str]:
    found: list[str] = []
    page = 0
    while True:
        batch = fetch(f"{API}/projects?pageSize={PROJECT_PAGE}&page={page}&sortDirection=ASC&sortFields=accession")
        if not isinstance(batch, list) or not batch:
            break
        found.extend(str(p.get("accession", "")) for p in batch if isinstance(p, dict))
        if limit and len(found) >= limit:
            break
        if len(batch) < PROJECT_PAGE:
            break
        page += 1
    return sorted({a for a in found if ACCESSION.fullmatch(a)})[: limit or None]


def project_files(fetch: Fetcher, accession: str, cache: Path | None) -> dict[str, Any]:
    """Return {'accession','expected','files','verified'}; cached when verified."""
    target = cache / f"{accession}.json" if cache else None
    if target and target.is_file():
        try:
            cached = json.loads(target.read_text(encoding="utf-8"))
            if cached.get("verified"):
                return cached
        except ValueError:
            pass
    expected = fetch(f"{API}/projects/{accession}/files/count")
    expected = int(expected) if isinstance(expected, (int, str)) and str(expected).isdigit() else None
    records: list[dict] = []
    page = 0
    while True:
        batch = fetch(f"{API}/projects/{accession}/files?pageSize={PAGE}&page={page}")
        if not isinstance(batch, list) or not batch:
            break
        records.extend(r for r in batch if isinstance(r, dict))
        if len(batch) < PAGE:
            break
        page += 1
    files = []
    for r in records:
        category = r.get("fileCategory") or {}
        files.append({"file_name": str(r.get("fileName") or ""),
                      "file_size_bytes": r.get("fileSizeBytes"),
                      "checksum": str(r.get("checksum") or ""),
                      "file_category": str(category.get("value") or "") if isinstance(category, dict) else ""})
    names = [f["file_name"] for f in files]
    result = {"accession": accession, "expected": expected, "listed": len(files),
              "unique_names": len(set(names)), "files": files,
              "verified": expected is not None and expected == len(files) == len(set(names))}
    if target:
        partial = target.with_suffix(f".json.partial-{os.getpid()}")
        partial.write_text(json.dumps(result, sort_keys=True), encoding="utf-8")
        partial.replace(target)
    return result


def load_processed(index: Path | None) -> dict[tuple, str]:
    """(accession, file casefold) -> 'DONE' | 'FAILED_PREVIOUSLY' from export index.jsonl."""
    state: dict[tuple, str] = {}
    if not index:
        return state
    for line in index.read_text(encoding="utf-8").splitlines():
        row = json.loads(line)
        key = (row.get("accession"), str(row.get("source_file", "")).casefold())
        if row.get("role") == "canonical":
            state[key] = "DONE"
        elif row.get("role") == "failed" and state.get(key) != "DONE":
            state[key] = "FAILED_PREVIOUSLY"
    return state


def decide(project: dict[str, Any], processed: dict[tuple, str]) -> list[dict]:
    accession = project["accession"]
    rows = []
    for f in project["files"]:
        size = f["file_size_bytes"]
        size = int(size) if isinstance(size, int) or (isinstance(size, str) and size.isdigit()) else None
        kind = file_format(f["file_name"])
        rows.append(dict(accession=accession, file_name=f["file_name"], file_category=f["file_category"],
                         file_format=kind, file_size_bytes=size, checksum=f["checksum"],
                         run_stem=run_stem(f["file_name"]), decision="", reason="",
                         memory_gb=memory_tier_gb(size), representation_evidence=""))
    vendor_stems = {r["run_stem"] for r in rows if r["file_format"] in VENDOR_FORMATS}
    unsupported_vendor_stems = {r["run_stem"] for r in rows
                                if r["file_format"] in VENDOR_FORMATS and r["file_format"] != "thermo_raw"}
    for r in rows:
        kind = r["file_format"]
        done = processed.get((accession, r["file_name"].casefold()))
        if kind not in SPECTRUM_FORMATS:
            r["decision"], r["reason"] = "SKIP", "SKIP_NOT_A_SPECTRUM_FILE"
        elif kind not in ("thermo_raw", "mzml"):
            r["decision"], r["reason"] = "SKIP", SKIP_REASON.get(kind, "SKIP_UNCLASSIFIED")
        elif not project["verified"]:
            r["decision"], r["reason"] = "SKIP", "SKIP_PROJECT_LISTING_NOT_VERIFIED"
        elif not SAFE_NAME.fullmatch(r["file_name"]):
            r["decision"], r["reason"] = "SKIP", "SKIP_UNSAFE_FILE_NAME"
        elif not r["file_size_bytes"]:
            r["decision"], r["reason"] = "SKIP", "SKIP_SIZE_UNKNOWN_OR_ZERO"
        elif kind == "thermo_raw":
            r["decision"] = done or "QUEUE_NATIVE_RAW"
            r["reason"] = "VENDOR_RAW_PREFERRED"
        elif kind == "mzml":
            if r["run_stem"] in vendor_stems:
                r["decision"] = "SKIP"
                r["reason"] = ("SKIP_MZML_VENDOR_RAW_PREFERRED"
                               if r["run_stem"] not in unsupported_vendor_stems
                               else "SKIP_MZML_VENDOR_PRESENT_BUT_UNSUPPORTED_POLICY_PENDING")
                r["representation_evidence"] = "stem_match_unverified"
            else:
                r["decision"] = done or "QUEUE_MZML_NO_VENDOR"
                r["reason"] = "NO_VENDOR_FILE_WITH_SAME_STEM"
                r["representation_evidence"] = "stem_match_unverified"
        else:
            r["decision"], r["reason"] = "SKIP", "SKIP_UNCLASSIFIED"
        if r["decision"] in ("DONE", "FAILED_PREVIOUSLY"):
            r["reason"] = f"{r['reason']};{r['decision']}"
    return rows


def write_tsv(path: Path, fields: list[str], rows: list[dict]) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def build(accessions: list[str], fetch: Fetcher, outdir: Path, *, cache: Path | None = None,
          processed_index: Path | None = None, sleep: float = 0.0, workers: int = 1,
          progress: Callable[[int, int], None] | None = None) -> dict[str, Any]:
    if outdir.exists() or outdir.is_symlink():
        raise ValueError(f"Refusing existing output directory: {outdir}")
    if cache:
        cache.mkdir(parents=True, exist_ok=True)
    processed = load_processed(processed_index)
    rows: list[dict] = []
    projects: list[dict] = []

    def one(accession: str) -> tuple[dict, str]:
        try:
            result = project_files(fetch, accession, cache), ""
        except Exception as err:  # one bad project must not stop the crawl
            result = ({"accession": accession, "expected": None, "listed": 0, "unique_names": 0,
                       "files": [], "verified": False}, type(err).__name__)
        if sleep:
            time.sleep(sleep)
        return result

    with ThreadPoolExecutor(max_workers=max(1, workers)) as pool:
        for i, (project, error) in enumerate(pool.map(one, accessions), 1):
            projects.append({k: project[k] for k in ("accession", "expected", "listed", "unique_names",
                                                     "verified")} | {"error": error})
            rows.extend(decide(project, processed))
            if progress:
                progress(i, len(accessions))
    stage = outdir.with_name(f".{outdir.name}.partial-{os.getpid()}")
    stage.mkdir(parents=True)
    write_tsv(stage / "inventory.tsv", FIELDS, rows)
    queue = [r for r in rows if r["decision"].startswith("QUEUE_")]
    write_tsv(stage / "queue.tsv", ["accession", "file_name", "file_format", "file_size_bytes", "checksum",
                                    "memory_gb", "decision"], queue)
    write_tsv(stage / "projects.tsv", ["accession", "expected", "listed", "unique_names", "verified", "error"],
              projects)

    def count(key: Callable[[dict], Any], items: list[dict]) -> dict[str, int]:
        out: dict[str, int] = {}
        for item in items:
            out[str(key(item))] = out.get(str(key(item)), 0) + 1
        return dict(sorted(out.items()))

    def tb(items: list[dict]) -> float:
        return round(sum(r["file_size_bytes"] or 0 for r in items) / 1e12, 3)

    spectrum = [r for r in rows if r["file_format"] in SPECTRUM_FORMATS]
    summary = {"schema_version": 1, "api": API, "projects": len(projects),
               "projects_verified": sum(1 for p in projects if p["verified"]),
               "projects_with_error": sum(1 for p in projects if p["error"]),
               "files": len(rows), "spectrum_files": len(spectrum), "spectrum_tb": tb(spectrum),
               "by_format": count(lambda r: r["file_format"], spectrum),
               "by_decision": count(lambda r: r["decision"], rows),
               "by_reason": count(lambda r: r["reason"], [r for r in rows if r["decision"] == "SKIP"
                                                          and r["reason"] != "SKIP_NOT_A_SPECTRUM_FILE"]),
               "queue_files": len(queue), "queue_tb": tb(queue),
               "queue_memory_gb": count(lambda r: r["memory_gb"], queue),
               "queue_projects": len({r["accession"] for r in queue}),
               "complete_study_asserted": False, "representation_evidence": "stem_match_unverified"}
    (stage / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    stage.rename(outdir)
    return summary


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", type=Path, required=True, help="New output directory (must not exist)")
    ap.add_argument("--cache", type=Path, help="Per-project JSON cache (resumable crawl)")
    ap.add_argument("--accession", action="append", default=[], help="Limit to these accessions")
    ap.add_argument("--max-projects", type=int, help="Take the first N projects in API listing order")
    ap.add_argument("--processed-index", type=Path, help="export_run_index.py index.jsonl for resume")
    ap.add_argument("--sleep", type=float, default=0.2, help="Seconds after each project, per worker")
    ap.add_argument("--workers", type=int, default=1, help="Parallel project fetches (keep small, <= 4)")
    args = ap.parse_args()
    try:
        accessions = sorted({a.upper() for a in args.accession}) or list_projects(http_fetch, args.max_projects)
        if any(not ACCESSION.fullmatch(a) for a in accessions):
            raise ValueError("Invalid accession")

        def progress(i: int, n: int) -> None:
            if i % 50 == 0 or i == n:
                print(f"JCODE_PROGRESS {json.dumps({'done': i, 'total': n})}", file=sys.stderr, flush=True)

        print(json.dumps(build(accessions, http_fetch, args.outdir, cache=args.cache,
                               processed_index=args.processed_index, sleep=args.sleep,
                               workers=min(max(args.workers, 1), 4), progress=progress),
                         indent=2, sort_keys=True))
        return 0
    except (ValueError, OSError) as err:
        print(json.dumps({"inventory": "BLOCKED", "reason": f"{type(err).__name__}: {err}"}), file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
