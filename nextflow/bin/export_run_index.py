#!/usr/bin/env python3
"""Build a run-level QC export and index from PRIDEQC Nextflow results. No DB writes.

Reads receipts under one or more results roots:
    <root>/<TAG>/<PXD>/_run_outputs/ready/units/<UNIT>/result.json
and writes a NEW export directory (refuses an existing one):
    <export>/<PXD>/<source_file>.mzQC          verbatim copy (receipt SHA256 kept)
    <export>/<PXD>/<source_file>.obo           verbatim copy
    <export>/<PXD>/<source_file>.summary.json  companion annotations, if present
    <export>/<PXD>/accession.json              runs of this accession in the export
    <export>/index.jsonl                       one row per receipt (canonical/superseded/failed/rejected)
    <export>/runs.ndjson                       one Mongo-ready document per canonical run
    <export>/export-manifest.json              counts and SHA256 of the index files

Rules:
- A success receipt is exported only if the mzQC and OBO SHA256 match the receipt.
- One canonical run per (accession, source file). Prefer receipts with computed
  source SHA256 evidence, then the newest mzQC creationDate, then the run tag name.
  Other successes for the same file become `superseded` (kept in the index only).
- A failure is `failed` unless the same file has a success elsewhere (`superseded_failure`).
- The same source SHA256 under a different (accession, file) is flagged, not merged.
- Accession completeness is NOT asserted: every accession is `partial_unverified`.
Python 3.9 compatible (Codon host python3).
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import sys
from pathlib import Path
from typing import Any

SCHEMA_VERSION = 1
ACCESSION = re.compile(r"PXD\d{6,}\Z")
SAFE_NAME = re.compile(r"[A-Za-z0-9_.+\-]{1,255}\Z")
MS1_COUNT, MS2_COUNT = "MS:4000059", "MS:4000060"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _json(path: Path) -> Any:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None


def _software(mzqc: dict) -> dict:
    found = {}
    for run in mzqc.get("runQualities", [])[:1]:
        for item in run.get("metadata", {}).get("analysisSoftware", []):
            if item.get("name") == "prideqc":
                found["prideqc_version"] = item.get("version")
            elif item.get("value") == "pyOpenMS":
                found["pyopenms_version"] = item.get("version")
    return found


def _metric_value(mzqc: dict, accession: str) -> Any:
    for run in mzqc.get("runQualities", [])[:1]:
        for metric in run.get("qualityMetrics", []):
            if metric.get("accession") == accession:
                return metric.get("value")
    return None


def collect(roots: list) -> list:
    """Return one candidate record per receipt. Never raises on a bad receipt."""
    records = []
    for root in roots:
        for receipt in sorted(Path(root).glob("*/*/_run_outputs/ready/units/*/result.json")):
            unit_dir = receipt.parent
            accession_dir = unit_dir.parents[3]
            tag = accession_dir.parent.name
            data = _json(receipt)
            base: dict[str, Any] = {"run_tag": tag, "results_root": str(root),
                                    "receipt": str(receipt.relative_to(root))}
            if not isinstance(data, dict):
                records.append(dict(base, role="rejected", reject_reason="unreadable_receipt"))
                continue
            accession = str(data.get("accession", ""))
            source_file = Path(str(data.get("source_member") or data.get("source_archive") or "")).name
            rec: dict[str, Any] = dict(base, accession=accession, source_file=source_file,
                       unit_id=data.get("unit_id"), task_id=data.get("task_id"),
                       status=data.get("status"), failure_stage=data.get("failure_stage"),
                       error_type=data.get("error_type"), exit_code=data.get("exit_code"),
                       reader_failure_kind=data.get("reader_failure_kind"))
            if not ACCESSION.fullmatch(accession) or accession_dir.name != accession:
                records.append(dict(rec, role="rejected", reject_reason="accession_mismatch"))
                continue
            if not SAFE_NAME.fullmatch(source_file) or source_file in {".", ".."}:
                records.append(dict(rec, role="rejected", reject_reason="unsafe_source_name"))
                continue
            status_file = accession_dir / "_archive_status" / f"{data.get('task_id')}.prepare-status.json"
            evidence = (_json(status_file) or {}).get("source_evidence") or {}
            rec["source_sha256"] = evidence.get("sha256")
            rec["source_size_bytes"] = evidence.get("size_bytes")
            rec["source_evidence"] = "computed_on_scratch" if evidence.get("sha256") else "absent"
            if data.get("status") != "success":
                records.append(dict(rec, role="failed"))
                continue
            qc = unit_dir / "qc"
            mzqc_path = qc / str(data.get("mzqc_file", ""))
            obo_path = qc / str(data.get("obo_file", ""))
            if not (mzqc_path.is_file() and obo_path.is_file()) or mzqc_path.is_symlink() or obo_path.is_symlink():
                records.append(dict(rec, role="rejected", reject_reason="missing_qc_files"))
                continue
            if sha256(mzqc_path) != data.get("mzqc_sha256") or sha256(obo_path) != data.get("obo_sha256"):
                records.append(dict(rec, role="rejected", reject_reason="qc_sha256_mismatch"))
                continue
            mzqc = (_json(mzqc_path) or {}).get("mzQC") or {}
            summary = sorted(qc.glob("*.summary.json"))
            rec.update(_software(mzqc))
            rec.update(mzqc_sha256=data["mzqc_sha256"], obo_sha256=data["obo_sha256"],
                       mzqc_creation_date=mzqc.get("creationDate"),
                       mzqc_version=mzqc.get("version"),
                       metric_count=sum(len(r.get("qualityMetrics", [])) for r in mzqc.get("runQualities", [])),
                       ms1_spectra=_metric_value(mzqc, MS1_COUNT),
                       ms2_spectra=_metric_value(mzqc, MS2_COUNT),
                       _mzqc=str(mzqc_path), _obo=str(obo_path),
                       _summary=str(summary[0]) if len(summary) == 1 else None)
            records.append(dict(rec, role="candidate"))
    return records


def _rank(rec: dict) -> tuple:
    return (rec.get("source_evidence") == "computed_on_scratch",
            rec.get("mzqc_creation_date") or "", rec.get("run_tag") or "")


def resolve(records: list) -> list:
    """Assign canonical/superseded roles and duplicate-bytes flags in place."""
    groups: dict = {}
    for rec in records:
        if rec.get("role") in ("candidate", "failed"):
            groups.setdefault((rec["accession"], rec["source_file"].casefold()), []).append(rec)
    for items in groups.values():
        successes = sorted((r for r in items if r["role"] == "candidate"), key=_rank, reverse=True)
        for i, rec in enumerate(successes):
            rec["role"] = "canonical" if i == 0 else "superseded"
            if i:
                rec["superseded_by"] = successes[0]["mzqc_sha256"]
        if successes:
            for rec in items:
                if rec["role"] == "failed":
                    rec["role"] = "superseded_failure"
                    rec["superseded_by"] = successes[0]["mzqc_sha256"]
    by_bytes: dict = {}
    for rec in records:
        if rec.get("role") == "canonical" and rec.get("source_sha256"):
            by_bytes.setdefault(rec["source_sha256"], []).append(rec)
    for items in by_bytes.values():
        if len(items) > 1:
            keys = sorted(f"{r['accession']}/{r['source_file']}" for r in items)
            for rec in items:
                rec["duplicate_source_bytes"] = [k for k in keys if k != f"{rec['accession']}/{rec['source_file']}"]
    return records


def _public(rec: dict) -> dict:
    return {k: v for k, v in sorted(rec.items()) if not k.startswith("_")}


def _write_lines(path: Path, rows: list) -> None:
    with path.open("x", encoding="utf-8") as handle:
        for row in rows:
            handle.write(json.dumps(row, sort_keys=True, ensure_ascii=False, allow_nan=False) + "\n")


def export(roots: list, outdir: Path) -> dict:
    outdir = Path(outdir)
    if outdir.exists() or outdir.is_symlink():
        raise ValueError(f"Refusing existing export directory: {outdir}")
    records = resolve(collect([Path(r) for r in roots]))
    stage = outdir.with_name(f".{outdir.name}.partial-{os.getpid()}")
    stage.mkdir(parents=True)
    try:
        runs = []
        per_accession: dict = {}
        for rec in records:
            if rec.get("role") != "canonical":
                continue
            target = stage / rec["accession"]
            target.mkdir(exist_ok=True)
            name = rec["source_file"]
            shutil.copyfile(rec["_mzqc"], target / f"{name}.mzQC")
            shutil.copyfile(rec["_obo"], target / f"{name}.obo")
            if sha256(target / f"{name}.mzQC") != rec["mzqc_sha256"]:
                raise ValueError(f"Copy integrity failure: {rec['accession']}/{name}")
            rec["mzqc_path"] = f"{rec['accession']}/{name}.mzQC"
            rec["obo_path"] = f"{rec['accession']}/{name}.obo"
            if rec.get("_summary"):
                shutil.copyfile(rec["_summary"], target / f"{name}.summary.json")
                rec["summary_path"] = f"{rec['accession']}/{name}.summary.json"
                rec["summary_sha256"] = sha256(target / f"{name}.summary.json")
            mzqc = json.loads((target / f"{name}.mzQC").read_text(encoding="utf-8"))["mzQC"]
            metrics = [{"accession": m.get("accession"), "name": m.get("name"),
                        "value": m.get("value"), "unit": m.get("unit")}
                       for q in mzqc.get("runQualities", []) for m in q.get("qualityMetrics", [])]
            run_id = hashlib.sha256(f"{rec['accession']}|{name}".encode()).hexdigest()
            runs.append({"_id": run_id, "schema_version": SCHEMA_VERSION,
                         "pxd_accession": rec["accession"], "source_file": name,
                         "source_sha256": rec.get("source_sha256"),
                         "source_size_bytes": rec.get("source_size_bytes"),
                         "mzqc_sha256": rec["mzqc_sha256"], "mzqc_path": rec["mzqc_path"],
                         "prideqc_version": rec.get("prideqc_version"),
                         "pyopenms_version": rec.get("pyopenms_version"),
                         "mzqc_creation_date": rec.get("mzqc_creation_date"),
                         "run_tag": rec["run_tag"], "metrics": metrics,
                         "scientific_approval": False})
            per_accession.setdefault(rec["accession"], []).append(rec)
        for rec in records:
            if rec.get("role") == "failed":
                per_accession.setdefault(rec["accession"], []).append(rec)
                (stage / rec["accession"]).mkdir(exist_ok=True)
        for accession, items in sorted(per_accession.items()):
            doc = {"schema_version": SCHEMA_VERSION, "pxd_accession": accession,
                   "inventory_scope": "partial_unverified", "complete_study": False,
                   "runs": [_public(r) for r in sorted(items, key=lambda r: r["source_file"].casefold())]}
            (stage / accession / "accession.json").write_text(
                json.dumps(doc, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        rows = [_public(r) for r in sorted(records, key=lambda r: (r.get("accession", ""), r.get("source_file", ""),
                                                                  r.get("role", ""), r.get("run_tag", "")))]
        _write_lines(stage / "index.jsonl", rows)
        _write_lines(stage / "runs.ndjson", sorted(runs, key=lambda r: (r["pxd_accession"], r["source_file"])))
        roles: dict = {}
        for rec in records:
            roles[rec.get("role")] = roles.get(rec.get("role"), 0) + 1
        manifest = {"schema_version": SCHEMA_VERSION, "results_roots": [str(r) for r in roots],
                    "receipts": len(records), "roles": dict(sorted(roles.items())),
                    "canonical_runs": len(runs), "accessions": len(per_accession),
                    "canonical_with_source_sha256": sum(1 for r in runs if r["source_sha256"]),
                    "duplicate_source_bytes_groups": len({r["source_sha256"] for r in records
                                                          if r.get("duplicate_source_bytes")}),
                    "index_sha256": sha256(stage / "index.jsonl"),
                    "runs_ndjson_sha256": sha256(stage / "runs.ndjson"),
                    "mongo_write_performed": False, "complete_study_asserted": False}
        (stage / "export-manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n",
                                                    encoding="utf-8")
        stage.rename(outdir)
    except BaseException:
        shutil.rmtree(stage, ignore_errors=True)
        raise
    return manifest


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results-root", type=Path, action="append", required=True,
                    help="Directory containing <TAG>/<PXD>/... result trees (repeatable)")
    ap.add_argument("--outdir", type=Path, required=True, help="New export directory (must not exist)")
    args = ap.parse_args()
    try:
        print(json.dumps(export(args.results_root, args.outdir), indent=2, sort_keys=True))
        return 0
    except (ValueError, OSError) as err:
        print(json.dumps({"export": "BLOCKED", "reason": f"{type(err).__name__}: {err}"}), file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
