#!/usr/bin/env python3
"""Publish prideQC results to a public folder (one folder per PRIDE project). No DB writes.

Reads receipts under one or more results roots (the same rules as export_run_index.py:
the receipt SHA-256 must match, one canonical run per project file) and writes:

    <dest>/README.md
    <dest>/manifest.json                 counts, generation time, sources
    <dest>/stats/archive-stats.json      archive-wide statistics (qc_archive_stats.py)
    <dest>/stats/projects.ndjson         one summary per project (qc_project_summary.py)
    <dest>/PXDnnnnnn/project.json        project summary, one compact entry per run, failed runs
    <dest>/PXDnnnnnn/<file>.mzQC         the mzQC of one run (verbatim)
    <dest>/PXDnnnnnn/<file>.obo          the local vocabulary of that mzQC (verbatim)

Each file is written to a temporary name and then renamed, so a reader never sees a part
of a file. An mzQC that is already in the destination with the same SHA-256 (from the
previous project.json) is not written again. Nothing is deleted.

--dry-run reads and checks everything but writes nothing.
Python 3.9 compatible (Codon host python3).
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import importlib.util
import json
import os
import re
import sys
from multiprocessing import Pool
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent


def _load(name: str, file: str) -> Any:
    spec = importlib.util.spec_from_file_location(name, HERE / file)
    module = importlib.util.module_from_spec(spec)  # type: ignore[arg-type]
    spec.loader.exec_module(module)  # type: ignore[union-attr]
    return module


IDX = _load("prideqc_publish_idx", "export_run_index.py")
QS = _load("prideqc_publish_stats", "qc_archive_stats.py")
PS = _load("prideqc_publish_projects", "qc_project_summary.py")

SCHEMA_VERSION = 1
DEFAULT_TAGS = r"^prideqc-(prod|flight|mzml-screen)"
FAILED_KEYS = ("failure_stage", "error_type", "reader_failure_kind", "download_failure_kind",
               "download_http_status")
README = """# prideQC: run QC of PRIDE Archive datasets

prideQC (https://github.com/PRIDE-Archive/prideQC) is an identification-free QC pipeline.
It makes one mzQC file (HUPO-PSI mzQC 1.0.0) for each MS run of a PRIDE dataset.

The QC is a preview. It is not yet scientifically validated (scientificApproval is false).

## Folder structure

    manifest.json                 counts, generation time and sources of this release
    stats/archive-stats.json      statistics over all runs and projects (histograms, categories)
    stats/projects.ndjson         one line for each project: the project summary of project.json
    PXDnnnnnn/project.json        the project summary and one entry for each run
    PXDnnnnnn/<file>.mzQC         the mzQC of one run (<file> is the PRIDE file name)
    PXDnnnnnn/<file>.obo          the local vocabulary (QCPRIDE terms) that this mzQC uses

## project.json

- summary: number of runs with QC, acquisition methods, MS2 resolutions, and the median
  run duration (min), MS1 spectra and MS2 spectra per run.
- runs: one entry for each run with QC. mzqc and obo are file names in the same folder.
  mzqcSha256 is the SHA-256 of the mzQC file. metrics holds the numeric metrics and the
  short form of the metric objects. The full metrics are in the mzQC file.
- failedRuns: files that prideQC tried but for which it has no QC, with the failure stage.

A project without QC has no folder. A project folder can get more runs in a later release.
"""


def _now() -> str:
    return dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")  # noqa: UP017


def _sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _write(path: Path, data: bytes) -> None:
    """Write data to path through a temporary file in the same folder (atomic rename)."""
    tmp = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    with tmp.open("wb") as handle:
        handle.write(data)
    os.replace(tmp, path)


def compact_value(value: Any) -> Any:
    """The short form of a metric value (the same rules as the pride-mongo-ws run summary).

    Numbers and lists of numbers are kept. From an object with suggested_tolerance or cv_term,
    only suggested_tolerance, unit, confidence, classification, resolution_regime and cv_term
    are kept. Anything else returns None (the metric is left out).
    """
    if isinstance(value, bool):
        return None
    if isinstance(value, (int, float)):
        return value
    if isinstance(value, list) and value and all(
            isinstance(x, (int, float)) and not isinstance(x, bool) for x in value):
        return value
    if isinstance(value, dict) and ("suggested_tolerance" in value or "cv_term" in value):
        keys = ("suggested_tolerance", "unit", "confidence", "classification", "resolution_regime", "cv_term")
        return {k: value[k] for k in keys if k in value}
    return None


def compact_metrics(metrics: list) -> list:
    out = []
    for m in metrics:
        value = compact_value(m.get("value"))
        if value is None:
            continue
        item = {"accession": m.get("accession"), "name": m.get("name"), "value": value}
        if m.get("unit") is not None:
            item["unit"] = m.get("unit")
        out.append(item)
    return out


def _metrics(mzqc: dict) -> list:
    return [{"accession": m.get("accession"), "name": m.get("name"), "value": m.get("value"), "unit": m.get("unit")}
            for q in mzqc.get("runQualities", []) for m in q.get("qualityMetrics", [])]


def _collect(item: tuple) -> dict:
    root, receipt = item
    return IDX.collect_one(Path(root), Path(receipt))


def _previous(dest_dir: Path) -> dict:
    try:
        doc = json.loads((dest_dir / "project.json").read_text(encoding="utf-8"))
        return {r["sourceFile"]: r.get("mzqcSha256") for r in doc.get("runs", [])}
    except (OSError, ValueError, KeyError, TypeError):
        return {}


def publish_project(job: tuple) -> dict:
    """Write (or check, on a dry run) one project. Returns its summary and the stats input."""
    accession, canonical, failed, dest, dry_run, generated = job
    target = Path(dest) / accession
    previous = _previous(target)
    if not dry_run:
        target.mkdir(parents=True, exist_ok=True)
    runs, stats_rows, written, nbytes = [], [], 0, 0
    values: dict = {k: [] for k in PS.BUCKETS}
    cats: dict = {"acquisition_method": set(), "ms2_resolution": set()}
    for rec in sorted(canonical, key=lambda r: r["source_file"].casefold()):
        name = rec["source_file"]
        data = Path(rec["_mzqc"]).read_bytes()
        if _sha(data) != rec["mzqc_sha256"]:
            raise ValueError(f"mzQC changed after the check: {accession}/{name}")
        obo = Path(rec["_obo"]).read_bytes()
        if _sha(obo) != rec["obo_sha256"]:
            raise ValueError(f"OBO changed after the check: {accession}/{name}")
        metrics = _metrics(json.loads(data)["mzQC"])
        nums, run_cats, ptms = QS.run_values({"metrics": metrics})
        mass_shift = any(m.get("accession") == QS.MASS_SHIFTS for m in metrics)
        stats_rows.append((nums, run_cats, sorted(ptms, key=lambda p: (p[0], p[1])), mass_shift))
        for k in values:
            if k in nums:
                values[k].append(nums[k])
        for k in cats:
            cats[k].add(run_cats[k])
        mzqc_path, obo_path = target / f"{name}.mzQC", target / f"{name}.obo"
        same = previous.get(name) == rec["mzqc_sha256"] and mzqc_path.is_file() \
            and mzqc_path.stat().st_size == len(data) and obo_path.is_file()
        if not same and not dry_run:
            _write(obo_path, obo)
            _write(mzqc_path, data)
        if not same:
            written += 1
            nbytes += len(data) + len(obo)
        runs.append({
            "sourceFile": name, "mzqc": f"{name}.mzQC", "obo": f"{name}.obo",
            "mzqcSha256": rec["mzqc_sha256"], "oboSha256": rec["obo_sha256"],
            "sourceSha256": rec.get("source_sha256"), "sourceSizeBytes": rec.get("source_size_bytes"),
            "mzqcCreationDate": rec.get("mzqc_creation_date"), "mzqcVersion": rec.get("mzqc_version"),
            "prideqcVersion": rec.get("prideqc_version"), "pyopenmsVersion": rec.get("pyopenms_version"),
            "runTag": rec.get("run_tag"),
            "putativeModifications": [{"name": n, "unimod": u, "deltaMassDa": d} for n, u, d in
                                      sorted(ptms, key=lambda p: (p[0], p[1]))],
            "metrics": compact_metrics(metrics),
        })
    summary = PS._summary(accession, len(runs), values, cats)
    latest: dict = {}
    for r in sorted(failed, key=lambda r: str(r.get("run_tag") or "")):
        latest[r["source_file"].casefold()] = r  # one entry per file: the newest batch wins
    failed_runs = [dict({"sourceFile": r["source_file"]},
                        **{k: r.get(k) for k in FAILED_KEYS if r.get(k) is not None})
                   for r in sorted(latest.values(), key=lambda r: r["source_file"].casefold())]
    doc = {"schemaVersion": SCHEMA_VERSION, "pxdAccession": accession, "generatedAt": generated,
           "scientificApproval": False, "summary": summary, "runs": runs, "failedRuns": failed_runs}
    if not dry_run:
        _write(target / "project.json",
               (json.dumps(doc, ensure_ascii=False, allow_nan=False, separators=(",", ":")) + "\n").encode())
    return {"accession": accession, "summary": summary, "stats": stats_rows, "written": written,
            "bytes": nbytes, "failed": len(failed_runs)}


def plan(roots: list, workers: int, tag_regex: str, accessions: set) -> tuple:
    """Check all receipts and choose the canonical run of each project file."""
    pattern = re.compile(tag_regex)
    items = [(str(root), str(receipt)) for root, receipt in IDX.receipts(roots)
             if pattern.search(receipt.parents[5].name)
             and (not accessions or receipt.parents[4].name in accessions)]
    with Pool(workers) as pool:
        records = pool.map(_collect, items, chunksize=200)
    IDX.resolve(records)
    projects: dict = {}
    for rec in records:
        role = rec.get("role")
        if role in ("canonical", "failed"):
            projects.setdefault(rec["accession"], ([], []))[0 if role == "canonical" else 1].append(rec)
    roles: dict = {}
    for rec in records:
        roles[rec.get("role")] = roles.get(rec.get("role"), 0) + 1
    return projects, roles, len(items)


def publish(roots: list, dest: Path, workers: int = 8, tag_regex: str = DEFAULT_TAGS,
            accessions: set | None = None, dry_run: bool = False) -> dict:
    dest = Path(dest)
    if not dry_run and not dest.is_dir():
        raise ValueError(f"Destination folder does not exist: {dest}")
    generated = _now()
    projects, roles, receipts = plan(roots, workers, tag_regex, accessions or set())
    with_qc = sorted(a for a, (canonical, _) in projects.items() if canonical)
    jobs = [(a, projects[a][0], projects[a][1], str(dest), dry_run, generated) for a in with_qc]
    acc = QS.Accumulator()
    summaries, written, nbytes, failed_files = [], 0, 0, 0
    with Pool(workers) as pool:
        for result in pool.imap(publish_project, jobs, chunksize=1):
            for nums, run_cats, ptms, mass_shift in result["stats"]:
                acc.add_values(result["accession"], nums, run_cats, {tuple(p) for p in ptms}, mass_shift)
            summary = dict(result["summary"], schemaVersion=SCHEMA_VERSION, generatedAt=generated,
                           scientificApproval=False)
            summaries.append(summary)
            written += result["written"]
            nbytes += result["bytes"]
            failed_files += result["failed"]
    stats = acc.finish()
    stats_doc = {"id": "prideqc-archive-stats-" + generated.replace("-", "").replace(":", ""),
                 "generatedAt": generated, "source": "pride_qc", "scientificApproval": False, **stats}
    runs = sum(s["qcRuns"] for s in summaries)
    manifest = {"schemaVersion": SCHEMA_VERSION, "generatedAt": generated, "scientificApproval": False,
                "projects": len(summaries), "runs": runs,
                "failedFilesWithoutQc": failed_files,
                "projectsWithoutQc": len(projects) - len(with_qc),
                "receipts": receipts, "roles": dict(sorted((str(k), v) for k, v in roles.items())),
                "tagRegex": tag_regex,
                "mzqcFilesWritten": written, "bytesWritten": nbytes, "dryRun": dry_run}
    if not dry_run:
        (dest / "stats").mkdir(exist_ok=True)
        _write(dest / "stats" / "archive-stats.json",
               json.dumps(stats_doc, allow_nan=False, separators=(",", ":")).encode())
        _write(dest / "stats" / "projects.ndjson",
               "".join(json.dumps(s, allow_nan=False, separators=(",", ":")) + "\n" for s in summaries).encode())
        _write(dest / "README.md", README.encode())
        _write(dest / "manifest.json", (json.dumps(manifest, indent=2, sort_keys=True) + "\n").encode())
    return manifest


def main(argv: list | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results-root", type=Path, action="append", required=True,
                    help="Directory containing <TAG>/<PXD>/... result trees (repeatable)")
    ap.add_argument("--dest", type=Path, required=True, help="Public folder (must exist)")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--tag-regex", default=DEFAULT_TAGS, help=f"run tags to read (default {DEFAULT_TAGS})")
    ap.add_argument("--accession", action="append", default=[], help="only these projects (repeatable)")
    ap.add_argument("--dry-run", action="store_true", help="check everything, write nothing")
    args = ap.parse_args(argv)
    try:
        manifest = publish(args.results_root, args.dest, args.workers, args.tag_regex,
                           set(args.accession), args.dry_run)
    except (ValueError, OSError) as err:
        print(json.dumps({"event": "BLOCKED", "reason": f"{type(err).__name__}: {err}"}), file=sys.stderr)
        return 2
    print(json.dumps(dict(manifest, event="DONE"), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
