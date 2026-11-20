#!/usr/bin/env python3
"""List the PRIDE projects with complete prideQC metrics for all of their raw data. No DB writes.

Inputs:
  --inventory   pride_inventory.py inventory.tsv (one row per project file, with a decision)
  --results-root  Nextflow results with QC receipts (same rules as export_run_index.py: the
                  receipt SHA-256 must match, one canonical run per project file)

A project is COMPLETE when all of these are true:
  1. The inventory has at least one QC input for the project: a Thermo .raw file
     (QUEUE_NATIVE_RAW) or an mzML file that is the run representation (QUEUE_MZML_NO_VENDOR:
     no vendor file of the same run, or a vendor file without a supported reader).
  2. Each QC input has a canonical prideQC result (a success, also after an earlier failure).
  3. No spectrum file of the project is outside QC: a file without a supported reader (WIFF,
     Bruker, Shimadzu, mzXML, .raw.gz, archives) is accepted only when a QC input of the same
     run stem has QC (for example a WIFF file with an mzML of the same run).
  4. The project listing is verified and no spectrum file was skipped for a bad name or size.
  With --archives (inspect_archives.py output), an inspected ZIP is replaced by its members:
  each run member is a QC input, a member without a reader is a file without a reader, and
  a ZIP with only results or peak lists is not raw data.

Peak lists (MGF), result files and other files are not raw data and do not count.
Run stems match by file name; this is a heuristic (stem_match_unverified in the inventory).

Outputs:
  --output    TSV of the complete projects
  --excluded  optional TSV of the other projects with QC inputs, with the first reason
A JSON summary goes to stdout.
Python 3.9 compatible (Codon host python3).
"""
from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import re
import sys
from collections import Counter
from multiprocessing import Pool
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
_SPEC = importlib.util.spec_from_file_location("prideqc_complete_idx", HERE / "export_run_index.py")
IDX = importlib.util.module_from_spec(_SPEC)  # type: ignore[arg-type]
_SPEC.loader.exec_module(IDX)  # type: ignore[union-attr]

DEFAULT_TAGS = r"^prideqc-(prod|flight|mzml-screen)"
QC_INPUTS = {"QUEUE_NATIVE_RAW", "QUEUE_MZML_NO_VENDOR", "DONE", "FAILED_PREVIOUSLY", "ARCHIVE_MEMBER"}
# spectrum files that are not a QC input but are covered by a QC input of the same run
COVERED_OK = {"SKIP_MZML_VENDOR_RAW_PREFERRED"}
# files that are not raw data: they never block a project
NOT_RAW = {"SKIP_NOT_A_SPECTRUM_FILE", "SKIP_PEAK_LIST_NOT_A_RUN_REPRESENTATION"}
# the project cannot be complete when one of its spectrum files has this reason
ALWAYS_BLOCK = {"SKIP_PROJECT_LISTING_NOT_VERIFIED": "listing_not_verified",
                "SKIP_UNSAFE_FILE_NAME": "unsafe_file_name",
                "SKIP_SIZE_UNKNOWN_OR_ZERO": "size_unknown"}
NO_READER = {
    "SKIP_SCIEX_WIFF_NO_QUALIFIED_READER": "no_reader_sciex_wiff",
    "SKIP_SCIEX_WIFF_SCAN_COMPANION": "no_reader_sciex_wiff",
    "SKIP_BRUKER_D_NO_QUALIFIED_READER": "no_reader_bruker",
    "SKIP_BRUKER_NO_QUALIFIED_READER": "no_reader_bruker",
    "SKIP_SHIMADZU_NO_QUALIFIED_READER": "no_reader_shimadzu",
    "SKIP_MZXML_READER_NOT_QUALIFIED": "no_reader_mzxml",
    "SKIP_COMPRESSED_RAW_NOT_QUALIFIED": "no_reader_raw_gz",
    "SKIP_ARCHIVE_CONTENTS_NOT_INSPECTED": "archive_not_inspected",
    "SKIP_ARCHIVE_MEMBER_NO_READER": "no_reader_archive_member",
    "SKIP_UNCLASSIFIED": "unclassified_spectrum_file",
}
COLUMNS = ["pride_accession", "runs_with_qc", "vendor_raw_runs", "mzml_runs", "mzml_runs_with_unsupported_vendor",
           "covered_other_files", "runs_from_archives", "representation", "source_gb", "prideqc_versions", "run_tags"]
EXCLUDED_COLUMNS = ["pride_accession", "reason", "qc_inputs", "with_qc", "failed", "not_processed",
                    "files_without_reader"]


def _collect(item: tuple) -> dict:
    root, receipt = item
    return IDX.collect_one(Path(root), Path(receipt))


def qc_status(roots: list, workers: int, tag_regex: str) -> dict:
    """(accession, file casefold) -> canonical record, or the failed record when there is no success."""
    pattern = re.compile(tag_regex)
    items = [(str(root), str(receipt)) for root, receipt in IDX.receipts(roots)
             if pattern.search(receipt.parents[5].name)]
    with Pool(workers) as pool:
        records = pool.map(_collect, items, chunksize=200)
    IDX.resolve(records)
    status: dict = {}
    for rec in records:
        role = rec.get("role")
        key = (rec.get("accession"), str(rec.get("source_file", "")).casefold())
        if role == "canonical":
            status[key] = rec
        elif role == "failed" and key not in status:
            status[key] = rec
    return status


def read_inventory(path: Path) -> dict:
    """accession -> list of the inventory rows that are spectrum files (not results or peak lists)."""
    projects: dict = {}
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            reason = row["reason"].split(";")[0]
            if reason in NOT_RAW:
                continue
            projects.setdefault(row["accession"], []).append(dict(
                file=row["file_name"], stem=row["run_stem"], decision=row["decision"], reason=reason,
                size=int(row["file_size_bytes"]) if row["file_size_bytes"].isdigit() else 0))
    return projects


def read_archive_members(folder: Path) -> tuple:
    """inspect_archives.py output -> ({(accession, archive casefold): decision}, {same key: [members]})."""
    decisions: dict = {}
    members: dict = {}
    with (folder / "archives.tsv").open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            decisions[(row["accession"], row["file_name"].casefold())] = row["decision"]
    with (folder / "members.tsv").open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            members.setdefault((row["accession"], row["archive_file"].casefold()), []).append(row)
    return decisions, members


def expand_archives(projects: dict, decisions: dict, members: dict) -> dict:
    """Replace each inspected ZIP row with the rows of its members.

    An archive without run members and without unsupported members (results, peak lists)
    is not raw data and is dropped. A run member is a QC input (ARCHIVE_MEMBER). A member
    without a reader (Bruker .d, Waters .raw folder, WIFF, .raw.gz) blocks the project unless
    a QC input of the same run stem has QC. An archive that was not inspected stays blocking.
    """
    out: dict = {}
    for accession, rows in projects.items():
        new = []
        for r in rows:
            key = (accession, r["file"].casefold())
            decision = decisions.get(key) if r["reason"] == "SKIP_ARCHIVE_CONTENTS_NOT_INSPECTED" else None
            if decision is None or decision in ("SKIP_ARCHIVE_UNREADABLE", "SKIP_ARCHIVE_FORMAT_NOT_INSPECTED"):
                new.append(r)
                continue
            for m in members.get(key, []):
                name = m["member_file"]
                if m["kind"] == "run":
                    new.append(dict(file=name, stem=_stem(name), decision="ARCHIVE_MEMBER", archive=r["file"],
                                    reason="VENDOR_RAW_ARCHIVE_MEMBER" if m["member_format"] == "raw"
                                    else "MZML_ARCHIVE_MEMBER", size=int(m["member_size_bytes"] or 0)))
                else:
                    new.append(dict(file=name, stem=_stem(name), decision="SKIP", archive=r["file"],
                                    reason="SKIP_ARCHIVE_MEMBER_NO_READER", size=0))
        out[accession] = new
    return out


def _stem(name: str) -> str:
    lower = name.casefold().rstrip("/")
    for suffix in (".mzml.gz", ".raw.gz", ".wiff.scan", ".wiff2", ".wiff", ".mzml", ".raw", ".d"):
        if lower.endswith(suffix):
            return lower[: -len(suffix)]
    return lower


def classify(accession: str, rows: list, status: dict) -> tuple:
    """('complete', summary row) or ('excluded', excluded row)."""
    inputs = [r for r in rows if r["decision"] in QC_INPUTS]
    ok, failed, pending = [], [], []
    for r in inputs:
        rec = status.get((accession, r["file"].casefold()))
        if rec is None and r.get("archive"):
            # a failed archive download or listing has one receipt for the archive, not per member
            rec = status.get((accession, r["archive"].casefold()))
            rec = rec if rec is not None and rec.get("role") != "canonical" else None
        if rec is None:
            pending.append(r)
        elif rec.get("role") == "canonical":
            ok.append((r, rec))
        else:
            failed.append(r)
    covered = {r["stem"] for r, _ in ok}
    blocking: Counter = Counter()
    other_ok = 0
    for r in rows:
        if r["decision"] in QC_INPUTS:
            continue
        if r["reason"] in ALWAYS_BLOCK:
            blocking[ALWAYS_BLOCK[r["reason"]]] += 1
        elif r["reason"] in COVERED_OK or (r["reason"] in NO_READER and r["stem"] in covered):
            other_ok += 1
        elif r["reason"] in NO_READER:
            blocking[NO_READER[r["reason"]]] += 1
        else:
            blocking["unclassified_spectrum_file"] += 1
    if not inputs or failed or pending or blocking:
        if not inputs:
            reason = "no_qc_input"
        elif failed:
            reason = "failed_run"
        elif blocking:
            reason = sorted(blocking.items(), key=lambda kv: (-kv[1], kv[0]))[0][0]
        else:
            reason = "not_processed_yet"
        return "excluded", {"pride_accession": accession, "reason": reason, "qc_inputs": len(inputs),
                            "with_qc": len(ok), "failed": len(failed), "not_processed": len(pending),
                            "files_without_reader": sum(blocking.values())}
    raw = sum(1 for r, _ in ok if r["reason"].startswith("VENDOR_RAW"))
    from_archives = sum(1 for r, _ in ok if r.get("archive"))
    mzml = len(ok) - raw
    unsupported = sum(1 for r, _ in ok if r["reason"].startswith("VENDOR_FILE_HAS_NO_SUPPORTED_READER"))
    kinds = [k for k, n in (("vendor_raw", raw), ("mzml", mzml)) if n]
    return "complete", {
        "pride_accession": accession, "runs_with_qc": len(ok), "vendor_raw_runs": raw, "mzml_runs": mzml,
        "mzml_runs_with_unsupported_vendor": unsupported, "covered_other_files": other_ok,
        "runs_from_archives": from_archives,
        "representation": "+".join(kinds), "source_gb": round(sum(r["size"] for r, _ in ok) / 1e9, 3),
        "prideqc_versions": ";".join(sorted({str(rec.get("prideqc_version")) for _, rec in ok})),
        "run_tags": ";".join(sorted({str(rec.get("run_tag")) for _, rec in ok}))}


def write_tsv(path: Path, columns: list, rows: list) -> None:
    tmp = path.with_name(path.name + ".tmp")
    with tmp.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    tmp.replace(path)


def run(inventory: Path, roots: list, output: Path, excluded_path: Path | None = None, workers: int = 8,
        tag_regex: str = DEFAULT_TAGS, archives: Path | None = None) -> dict:
    projects = read_inventory(inventory)
    if archives is not None:
        projects = expand_archives(projects, *read_archive_members(archives))
    status = qc_status(roots, workers, tag_regex)
    complete, excluded = [], []
    for accession in sorted(projects):
        kind, row = classify(accession, projects[accession], status)
        (complete if kind == "complete" else excluded).append(row)
    write_tsv(output, COLUMNS, complete)
    if excluded_path:
        write_tsv(excluded_path, EXCLUDED_COLUMNS, excluded)
    with_inputs = [e for e in excluded if e["reason"] != "no_qc_input"]
    summary: dict[str, Any] = {
        "inventory_projects": len(projects),
        "projects_with_qc_inputs": len(complete) + len(with_inputs),
        "complete_projects": len(complete),
        "complete_runs": sum(r["runs_with_qc"] for r in complete),
        "complete_by_representation": dict(Counter(r["representation"] for r in complete)),
        "excluded_by_reason": dict(Counter(e["reason"] for e in excluded)),
        "output": str(output)}
    return summary


def main(argv: list | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--inventory", type=Path, required=True, help="pride_inventory.py inventory.tsv")
    ap.add_argument("--results-root", type=Path, action="append", required=True)
    ap.add_argument("--output", type=Path, required=True, help="TSV of the complete projects")
    ap.add_argument("--excluded", type=Path, help="TSV of the excluded projects with a reason")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--tag-regex", default=DEFAULT_TAGS)
    ap.add_argument("--archives", type=Path, help="inspect_archives.py output folder (archives.tsv, members.tsv)")
    args = ap.parse_args(argv)
    try:
        summary = run(args.inventory, args.results_root, args.output, args.excluded, args.workers, args.tag_regex,
                      args.archives)
    except (ValueError, OSError, KeyError) as err:
        print(json.dumps({"event": "BLOCKED", "reason": f"{type(err).__name__}: {err}"}), file=sys.stderr)
        return 2
    print(json.dumps(dict(summary, event="DONE"), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
