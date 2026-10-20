#!/usr/bin/env python3
"""Construct a fail-closed, one-representation-per-MS-run Nextflow manifest.

Input is a *curated complete per-project candidate inventory*, NOT the output of
PRIDE's RAW-only inventory endpoint. In mixed vendor/open-format datasets the
operator must review source MSRun/SDRF/conversion evidence and supply a stable
logical_run_id. Filename-stem equality alone is NOT evidence of equivalence.

No complete-study attestation is issued here, and no files are downloaded.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from collections import defaultdict
from pathlib import Path

ACC = re.compile(r"PXD[0-9]{6,}$")
NAME = re.compile(r"[A-Za-z0-9_.+\-]+$")
RUN_ID = re.compile(r"[A-Za-z0-9_.+\-]{1,128}$")
HEX = re.compile(r"[0-9a-f]{64}$")
FORMATS = (".mzml.gz", ".raw.gz", ".mzml", ".raw")
REQUIRED = {"pxd_accession", "file_name", "file_category", "size_bytes", "logical_run_id",
            "equivalence_evidence", "native_reader_verified"}
FIELDS = ["task_id", "pxd_accession", "archive_file", "sdrf_path", "needs_file_map",
          "logical_run_id", "representation_policy", "representation_inventory_sha256"]


def format_kind(filename: str) -> str:
    lowered = filename.casefold()
    if lowered.endswith((".mzml", ".mzml.gz")):
        return "mzml"
    if lowered.endswith((".raw", ".raw.gz")):
        return "native_raw"
    raise ValueError(f"Unsupported non-single-file format: {filename}")


def source_stem(filename: str) -> str:
    lower = filename.casefold()
    for suffix in FORMATS:
        if lower.endswith(suffix):
            return lower[:-len(suffix)]
    raise ValueError("Unsupported source suffix")


def file_hash(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as src:
        for chunk in iter(lambda: src.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _read(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as src:
        reader = csv.DictReader(src, delimiter="\t")
        if not REQUIRED.issubset(reader.fieldnames or []):
            raise ValueError(f"Missing inventory columns: {sorted(REQUIRED-set(reader.fieldnames or []))}")
        rows = list(reader)
    if not rows:
        raise ValueError("Empty candidate inventory")
    return rows


def _tsv(path: Path, columns: list[str], rows: list[dict]) -> None:
    with path.open("x", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def build(inventory: Path, outdir: Path, *, prefer: str = "native", max_size_gib: int = 100) -> dict:
    if prefer not in ("native", "mzml") or not (1 <= max_size_gib <= 500):
        raise ValueError("Invalid preference/size bound")
    if outdir.exists():
        raise ValueError("Preserving existing manifest directory")
    rows = _read(inventory)
    size_cap = max_size_gib * 1024**3
    per_run: dict[tuple[str, str], list[dict]] = defaultdict(list)
    per_stem: dict[tuple[str, str], set[str]] = defaultdict(set)
    seen_names: set[tuple[str, str]] = set()
    kinds: dict[str, set[str]] = defaultdict(set)
    for row in rows:
        pxd, name, run_id = (row[key].strip() for key in
                             ("pxd_accession", "file_name", "logical_run_id"))
        if not ACC.fullmatch(pxd) or not NAME.fullmatch(name) or name in (".", ".."):
            raise ValueError(f"Unsafe accession/source: {pxd} {name}")
        if not RUN_ID.fullmatch(run_id):
            raise ValueError(f"Missing/unsafe reviewed logical_run_id: {pxd}/{name}")
        cat = row["file_category"].strip().upper()
        if cat not in ("RAW", "PEAK"):
            raise ValueError(f"Source must be PRIDE RAW or PEAK category: {pxd}/{name}")
        nbytes = row["size_bytes"].strip()
        if not nbytes.isascii() or not nbytes.isdecimal() or not 0 < int(nbytes) <= size_cap:
            raise ValueError(f"Missing/out-of-range file size: {pxd}/{name}")
        verified = row["native_reader_verified"].strip().lower()
        if verified not in ("0", "1"):
            raise ValueError(f"native_reader_verified must be 0/1: {pxd}/{name}")
        kind = format_kind(name)
        if kind == "mzml" and verified == "1":
            raise ValueError("native_reader_verified is only for native vendor files")
        key = pxd, name.casefold()
        if key in seen_names:
            raise ValueError(f"Duplicate source filename: {pxd}/{name}")
        seen_names.add(key)
        group = pxd, run_id
        per_stem[(pxd, source_stem(name))].add(run_id)
        kinds[pxd].add(kind)
        per_run[group].append({**row, "pxd_accession": pxd, "file_name": name,
                               "logical_run_id": run_id, "format_kind": kind,
                               "native_reader_verified": verified})
    # Identical base-name files assigned to different acquisitions are ambiguous,
    # e.g. run.raw and run.mzML with conflicting user-supplied run IDs.
    if any(len(ids) != 1 for ids in per_stem.values()):
        raise ValueError("Same-stem representations assigned to different run IDs")
    # In mixed-format studies, every row MUST have a reviewed link to the
    # acquisition; a basename match or empty evidence can never authorize it.
    for pxd, _run_id in per_run:
        if kinds[pxd] == {"mzml", "native_raw"}:
            for row in per_run[(pxd, _run_id)]:
                evidence = row["equivalence_evidence"].strip()
                if not evidence.startswith(("MSRUN:", "SDRF:", "CONVERSION:", "CURATED:")) or len(evidence) < 8:
                    raise ValueError(f"Mixed-format project requires reviewed evidence per file: {pxd}")
    choices: list[dict] = []
    audit: list[dict] = []
    for (pxd, run_id), members in sorted(per_run.items()):
        # Multiple members in one group MUST share the same explicit evidence;
        # groups with independently asserted evidence are not deduplicated.
        proofs = {m["equivalence_evidence"].strip() for m in members}
        if len(members) > 1 and (len(proofs) != 1 or not next(iter(proofs))):
            raise ValueError(f"Unverified equivalent representations: {pxd}/{run_id}")
        ranked = sorted(members, key=lambda m: (
            0 if (m["format_kind"] == "native_raw" and m["native_reader_verified"] == "1" and prefer == "native")
            else 1 if m["format_kind"] == "mzml" else 2,
            m["file_name"].casefold()))
        chosen = ranked[0]
        if chosen["format_kind"] == "native_raw" and chosen["native_reader_verified"] != "1":
            raise ValueError(f"No verified native reader or supported mzML alternative: {pxd}/{run_id}")
        choices.append(chosen)
        for item in members:
            audit.append({"pxd_accession": pxd, "logical_run_id": run_id,
                          "file_name": item["file_name"], "format_kind": item["format_kind"],
                          "selected": "1" if item is chosen else "0",
                          "selected_source": chosen["file_name"],
                          "equivalence_evidence": item["equivalence_evidence"].strip()})
    checksum = file_hash(inventory)
    manifest = [{"task_id": i, "pxd_accession": r["pxd_accession"],
                 "archive_file": r["file_name"], "sdrf_path": "", "needs_file_map": "0",
                 "logical_run_id": r["logical_run_id"],
                 "representation_policy": "curated_single_primary_v1",
                 "representation_inventory_sha256": checksum}
                for i, r in enumerate(choices, 1)]
    accessions = sorted(kinds)
    reports = [{"report_id": i, "pxd_accession": pxd,
                "expected_file_tasks": sum(x["pxd_accession"] == pxd for x in manifest)}
               for i, pxd in enumerate(accessions, 1)]
    outdir.mkdir(parents=True)
    _tsv(outdir / "file-manifest.tsv", FIELDS, manifest)
    _tsv(outdir / "accession-reports.tsv", ["report_id", "pxd_accession", "expected_file_tasks"], reports)
    _tsv(outdir / "representation-decisions.tsv", ["pxd_accession", "logical_run_id", "file_name",
                                                   "format_kind", "selected", "selected_source",
                                                   "equivalence_evidence"], audit)
    info = {"schema_version": 1, "inventory_sha256": checksum,
            "selection_policy": "curated_single_primary_v1", "complete_study_inventory_verified": False,
            "require_operator_review": True, "study_merge_allowed": False,
            "source_rows": len(rows), "logical_runs": len(manifest),
            "suppressed_alternates": len(rows) - len(manifest), "accessions": accessions,
            "preferred_format": prefer, "max_file_size_gib": max_size_gib}
    (outdir / "representation-provenance.json").write_text(json.dumps(info, indent=2, sort_keys=True) + "\n")
    return info


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--inventory", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--prefer", choices=("native", "mzml"), default="native")
    ap.add_argument("--max-size-gib", type=int, default=100)
    args = ap.parse_args()
    print(json.dumps(build(args.inventory, args.outdir, prefer=args.prefer,
                           max_size_gib=args.max_size_gib), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
