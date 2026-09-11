#!/usr/bin/env python3
"""Resolve annotated SDRFs into one PRIDE archive file per Slurm task.

The SDRF's comment[data file] may describe a converted analysis file while
comment[file uri] points at the actual archive object.  Preserve both names so
prideQC can download the archive file and map the result back to the SDRF row.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import sys
import urllib.error
import urllib.parse
import urllib.request
from collections import defaultdict
from collections.abc import Iterable
from pathlib import Path

MISSING = {"", "not available", "not provided", "not applicable"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gt", required=True, type=Path)
    parser.add_argument("--sdrf-root", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--meta", required=True, type=Path)
    parser.add_argument(
        "--manifest-path-base",
        type=Path,
        help=(
            "Optional base path to strip from persisted file paths. "
            "Use this when manifest preparation runs inside a container bind mount."
        ),
    )
    parser.add_argument("--repo", default="bigbio/sdrf-annotated-datasets")
    parser.add_argument("--ref", default="main")
    return parser.parse_args()


def request_bytes(url: str, timeout: int = 120) -> bytes:
    request = urllib.request.Request(
        url,
        headers={"Accept": "application/vnd.github+json", "User-Agent": "prideQC-HPC"},
    )
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return response.read()


def unique_in_order(values: Iterable[str]) -> list[str]:
    seen: set[str] = set()
    result: list[str] = []
    for value in values:
        if value not in seen:
            seen.add(value)
            result.append(value)
    return result


def safe_basename(value: str, *, label: str, source: Path) -> str:
    value = value.strip()
    if value in {".", ".."} or any(char in value for char in "/\\:\x00\r\n\t"):
        raise ValueError(f"Unsafe/non-basename {label} in {source}: {value!r}")
    if not value:
        raise ValueError(f"Empty {label} in {source}")
    return value


def archive_name_from_uri(uri: str, *, source: Path) -> str | None:
    uri = uri.strip()
    if uri.casefold() in MISSING:
        return None
    parsed = urllib.parse.urlparse(uri)
    if parsed.scheme and parsed.scheme.casefold() not in {"ftp", "http", "https"}:
        return None
    basename = Path(urllib.parse.unquote(parsed.path)).name
    if not basename:
        return None
    return safe_basename(basename, label="comment[file uri] basename", source=source)


def load_gt(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "pxd_accession" not in reader.fieldnames:
            raise ValueError("Ground-truth TSV must contain pxd_accession.")
        rows = [dict(row) for row in reader]
    accessions = unique_in_order(
        row["pxd_accession"].strip().upper()
        for row in rows
        if row.get("pxd_accession", "").strip()
    )
    bad = [value for value in accessions if not re.fullmatch(r"PXD\d{6,}", value)]
    if bad:
        raise ValueError(f"Unexpected PRIDE accessions in ground truth: {bad}")
    return rows, accessions


def discover_urls(pxd: str, repo: str, ref: str) -> list[str]:
    quoted_ref = urllib.parse.quote(ref, safe="")
    url = f"https://api.github.com/repos/{repo}/contents/datasets/{pxd}?ref={quoted_ref}"
    entries = json.loads(request_bytes(url, timeout=60))
    urls = sorted(
        str(item["download_url"])
        for item in entries
        if item.get("type") == "file"
        and str(item.get("name", "")).lower().endswith(".sdrf.tsv")
        and item.get("download_url")
    )
    if not urls:
        raise RuntimeError(f"No annotated SDRF files found for {pxd}.")
    return urls


def fallback_urls(rows: list[dict[str, str]], pxd: str) -> list[str]:
    return unique_in_order(
        row.get("source_sdrf_url", "").strip()
        for row in rows
        if row.get("pxd_accession", "").strip().upper() == pxd
        and row.get("source_sdrf_url", "").strip()
    )


def infer_template(path: Path, gt_rows: list[dict[str, str]], url: str) -> str:
    if "dia" in path.name.casefold():
        return "dia-acquisition"
    for row in gt_rows:
        if row.get("source_sdrf_url", "").strip() == url:
            if row.get("acquisition_method", "").strip().casefold() == "dia":
                return "dia-acquisition"
    try:
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            fields = reader.fieldnames or []
            acquisition_fields = [
                field for field in fields if "acquisition method" in field.casefold()
            ]
            for index, row in enumerate(reader):
                if index >= 1000:
                    break
                for field in acquisition_fields:
                    value = (row.get(field) or "").casefold()
                    if "data-independent" in value or re.search(
                        r"(^|[^a-z])dia([^a-z]|$)", value
                    ):
                        return "dia-acquisition"
    except (OSError, csv.Error):
        pass
    return "ms-proteomics"


def file_records(path: Path) -> list[dict[str, str]]:
    """Return unique SDRF-file/archive-file mappings in SDRF order."""
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"SDRF has no header: {path}")
        lookup = {field.casefold(): field for field in reader.fieldnames}
        data_field = lookup.get("comment[data file]")
        uri_field = lookup.get("comment[file uri]")
        if data_field is None:
            raise ValueError(f"SDRF has no comment[data file] column: {path}")

        by_data_file: dict[str, dict[str, str]] = {}
        order: list[str] = []
        for row in reader:
            sdrf_name = (row.get(data_field) or "").strip()
            if sdrf_name.casefold() in MISSING:
                continue
            sdrf_name = safe_basename(
                sdrf_name, label="comment[data file]", source=path
            )
            file_uri = (row.get(uri_field) or "").strip() if uri_field else ""
            archive_name = archive_name_from_uri(file_uri, source=path) or sdrf_name
            record = {
                "sdrf_data_file": sdrf_name,
                "archive_file": archive_name,
                "file_uri": file_uri,
                "needs_file_map": "1" if archive_name != sdrf_name else "0",
            }
            previous = by_data_file.get(sdrf_name)
            if previous is None:
                by_data_file[sdrf_name] = record
                order.append(sdrf_name)
            elif previous["archive_file"] != archive_name:
                raise ValueError(
                    f"Ambiguous archive filenames for SDRF data file {sdrf_name!r} in {path}: "
                    f"{previous['archive_file']!r} vs {archive_name!r}"
                )

    records = [by_data_file[name] for name in order]
    if not records:
        raise ValueError(f"SDRF contains no usable comment[data file] values: {path}")
    return records


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def manifest_path(path: Path, base: Path | None) -> str:
    """Return a host-portable path for storage in the manifest."""
    resolved = path.resolve()
    if base is None:
        return str(resolved)
    resolved_base = base.resolve()
    try:
        return resolved.relative_to(resolved_base).as_posix()
    except ValueError as exc:
        raise ValueError(
            f"Manifest path {resolved} is outside manifest path base {resolved_base}"
        ) from exc


def main() -> int:
    args = parse_args()
    rows, accessions = load_gt(args.gt)
    args.sdrf_root.mkdir(parents=True, exist_ok=True)
    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    args.meta.parent.mkdir(parents=True, exist_ok=True)

    manifest_rows: list[dict[str, str]] = []
    per_pxd: dict[str, int] = defaultdict(int)
    mapped_count = 0
    sdrf_count = 0
    discovery_notes: list[str] = []

    for pxd in accessions:
        try:
            urls = discover_urls(pxd, args.repo, args.ref)
            discovery_notes.append(f"{pxd}\tgithub-api\t{len(urls)}")
        except (
            OSError,
            RuntimeError,
            urllib.error.HTTPError,
            urllib.error.URLError,
            json.JSONDecodeError,
        ) as exc:
            urls = fallback_urls(rows, pxd)
            discovery_notes.append(
                f"{pxd}\tgt-fallback\t{len(urls)}\t{type(exc).__name__}"
            )
        if not urls:
            raise RuntimeError(f"Could not resolve any annotated SDRF URL for {pxd}.")

        for url in urls:
            basename = Path(urllib.parse.urlparse(url).path).name
            if not basename.lower().endswith(".sdrf.tsv"):
                raise ValueError(f"Unexpected SDRF URL for {pxd}: {url}")
            target_dir = args.sdrf_root / pxd
            target_dir.mkdir(parents=True, exist_ok=True)
            target = target_dir / basename
            payload = request_bytes(url)
            if not payload:
                raise RuntimeError(f"Downloaded empty SDRF for {pxd}: {url}")
            target.write_bytes(payload)
            template = infer_template(target, rows, url)
            digest = sha256(target)
            records = file_records(target)
            sdrf_count += 1
            for record in records:
                if record["needs_file_map"] == "1":
                    mapped_count += 1
                manifest_rows.append(
                    {
                        "task_id": str(len(manifest_rows) + 1),
                        "pxd_accession": pxd,
                        "sdrf_basename": basename,
                        "sdrf_path": manifest_path(target, args.manifest_path_base),
                        "sdrf_sha256": digest,
                        "source_sdrf_url": url,
                        "sdrf_template": template,
                        **record,
                    }
                )
                per_pxd[pxd] += 1

    fieldnames = [
        "task_id",
        "pxd_accession",
        "sdrf_basename",
        "sdrf_path",
        "sdrf_sha256",
        "source_sdrf_url",
        "sdrf_template",
        "sdrf_data_file",
        "archive_file",
        "file_uri",
        "needs_file_map",
    ]
    tmp_manifest = args.manifest.with_suffix(args.manifest.suffix + ".tmp")
    with tmp_manifest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(manifest_rows)
    os.replace(tmp_manifest, args.manifest)

    with args.meta.open("w", encoding="utf-8") as handle:
        handle.write(f"gt={args.gt.resolve()}\n")
        handle.write(f"gt_sha256={sha256(args.gt)}\n")
        handle.write(f"sdrf_repo={args.repo}\n")
        handle.write(f"sdrf_ref={args.ref}\n")
        handle.write(f"accession_count={len(accessions)}\n")
        handle.write(f"sdrf_count={sdrf_count}\n")
        handle.write(f"file_task_count={len(manifest_rows)}\n")
        handle.write(f"archive_alias_count={mapped_count}\n")
        for pxd in accessions:
            handle.write(f"files_{pxd}={per_pxd[pxd]}\n")
        handle.write("\n[discovery]\n")
        handle.write("\n".join(discovery_notes) + "\n")

    print(f"accessions={len(accessions)}")
    print(f"sdrfs={sdrf_count}")
    print(f"file_tasks={len(manifest_rows)}")
    print(f"archive_aliases={mapped_count}")
    for pxd in accessions:
        print(f"{pxd}\t{per_pxd[pxd]}")
    print(f"manifest={args.manifest}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, RuntimeError, ValueError, csv.Error) as exc:
        print(f"prepare_prideqc_file_manifest: {exc}", file=sys.stderr)
        raise SystemExit(2) from exc
