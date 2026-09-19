#!/usr/bin/env python3
"""Build a prideQC file-array manifest directly from public PRIDE metadata.

This helper is intentionally metadata-only.  It is for public accessions where
no trusted SDRF is available and therefore leaves every SDRF provenance column
blank.  The output uses the existing prideQC file-array manifest schema so the
same Slurm launcher can analyze explicit repository files without fabricating
sample annotations.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import sys
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any

PRIDE_API = "https://www.ebi.ac.uk/pride/ws/archive/v3"
MANIFEST_FIELDS = [
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
INVENTORY_FIELDS = [
    "file_name",
    "file_category",
    "file_size_bytes",
    "ftp_uri",
    "aspera_uri",
    "selected",
]
NATURAL_PART_RE = re.compile(r"(\d+)")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--accession", required=True)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--meta", required=True, type=Path)
    parser.add_argument(
        "--inventory",
        type=Path,
        help="Optional complete PRIDE file inventory TSV (defaults next to the manifest).",
    )
    parser.add_argument(
        "--file",
        action="append",
        default=[],
        dest="files",
        help=(
            "Exact public RAW filename to include; repeat for a smoke subset. "
            "If omitted, every current PRIDE RAW file is included."
        ),
    )
    parser.add_argument("--pride-api", default=PRIDE_API)
    return parser.parse_args()


def request_json_with_headers(
    url: str,
    timeout: int = 300,
) -> tuple[object, dict[str, str]]:
    request = urllib.request.Request(
        url,
        headers={"Accept": "application/json", "User-Agent": "prideQC-HPC"},
    )
    with urllib.request.urlopen(request, timeout=timeout) as response:
        payload = json.loads(response.read())
        headers = {key.casefold(): value for key, value in response.headers.items()}
    return payload, headers


def _page_records(payload: object, accession: str) -> tuple[list[dict[str, Any]], int | None, int | None]:
    """Accept current bare-list v3 responses and legacy HAL wrappers."""
    if isinstance(payload, list):
        records = payload
        number = None
        total_pages = None
    elif isinstance(payload, dict):
        embedded = payload.get("_embedded")
        if isinstance(embedded, dict) and isinstance(embedded.get("files"), list):
            records = embedded["files"]
        elif isinstance(payload.get("files"), list):
            records = payload["files"]
        else:
            raise RuntimeError(f"Unexpected PRIDE project-file payload for {accession}.")
        page = payload.get("page")
        number = page.get("number") if isinstance(page, dict) else None
        total_pages = page.get("totalPages") if isinstance(page, dict) else None
    else:
        raise RuntimeError(f"Unexpected PRIDE project-file payload for {accession}.")

    typed: list[dict[str, Any]] = []
    for record in records:
        if not isinstance(record, dict) or not str(record.get("fileName") or "").strip():
            raise RuntimeError(f"PRIDE returned a project file without fileName for {accession}.")
        typed.append(record)
    return typed, _optional_int(number), _optional_int(total_pages)


def _optional_int(value: object) -> int | None:
    if value is None:
        return None
    try:
        return int(value)
    except (TypeError, ValueError) as exc:
        raise RuntimeError(f"Invalid PRIDE pagination value: {value!r}") from exc


def _total_records(headers: dict[str, str]) -> int | None:
    value = headers.get("total_records") or headers.get("total-records")
    return _optional_int(value)


def pride_project_files(accession: str, api_base: str = PRIDE_API) -> list[dict[str, Any]]:
    records_by_name: dict[str, dict[str, Any]] = {}
    page = 0
    page_size = 100
    expected_total: int | None = None
    while True:
        url = (
            f"{api_base.rstrip('/')}/projects/{urllib.parse.quote(accession)}/files"
            f"?pageSize={page_size}&page={page}"
        )
        payload, headers = request_json_with_headers(url)
        records, page_number, total_pages = _page_records(payload, accession)
        for record in records:
            name = str(record["fileName"]).strip()
            previous = records_by_name.get(name)
            if previous is not None and previous != record:
                raise RuntimeError(f"PRIDE returned conflicting metadata for {accession}/{name}.")
            records_by_name[name] = record

        header_total = _total_records(headers)
        if header_total is not None:
            expected_total = header_total
            if len(records_by_name) >= expected_total:
                break

        if total_pages is not None:
            next_page = (page_number + 1) if page_number is not None else page + 1
            if next_page >= total_pages:
                break
            page = next_page
            continue

        if not records:
            break
        page += 1

    if not records_by_name:
        raise RuntimeError(f"PRIDE returned no project files for {accession}.")
    if expected_total is not None and len(records_by_name) != expected_total:
        raise RuntimeError(
            f"PRIDE file pagination for {accession} returned {len(records_by_name)} unique "
            f"files but total_records={expected_total}."
        )
    return list(records_by_name.values())


def _category(record: dict[str, Any]) -> str:
    category = record.get("fileCategory")
    if isinstance(category, dict):
        value = category.get("value") or category.get("name")
    else:
        value = category
    return str(value or "").strip().upper()


def _public_location(record: dict[str, Any], protocol_name: str) -> str:
    locations = record.get("publicFileLocations") or []
    if not isinstance(locations, list):
        return ""
    for location in locations:
        if not isinstance(location, dict):
            continue
        if str(location.get("name") or "").strip().casefold() == protocol_name.casefold():
            return str(location.get("value") or "").strip()
    return ""


def _file_size(record: dict[str, Any]) -> int:
    value = record.get("fileSizeBytes")
    if value in (None, ""):
        return 0
    try:
        result = int(value)
    except (TypeError, ValueError) as exc:
        raise RuntimeError(
            f"Invalid fileSizeBytes for {record.get('fileName')!r}: {value!r}"
        ) from exc
    if result < 0:
        raise RuntimeError(f"Negative fileSizeBytes for {record.get('fileName')!r}.")
    return result


def _is_raw(record: dict[str, Any]) -> bool:
    category = _category(record)
    if category:
        return category == "RAW"
    return str(record.get("fileName") or "").casefold().endswith(".raw")


def _natural_key(name: str) -> tuple[tuple[int, object], ...]:
    return tuple(
        (0, int(part)) if part.isdigit() else (1, part.casefold())
        for part in NATURAL_PART_RE.split(name)
        if part
    )


def select_raw_files(
    records: list[dict[str, Any]], requested: list[str]
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    raw_records = [record for record in records if _is_raw(record)]
    raw_by_name = {str(record["fileName"]).strip(): record for record in raw_records}
    if not raw_records:
        raise RuntimeError("PRIDE metadata contains no RAW files.")

    if requested:
        duplicates = sorted({name for name in requested if requested.count(name) > 1})
        if duplicates:
            raise ValueError(f"Duplicate --file selections: {duplicates}")
        missing = [name for name in requested if name not in raw_by_name]
        if missing:
            raise ValueError(
                "Requested file(s) are not current PRIDE RAW files: " + ", ".join(missing)
            )
        selected = [raw_by_name[name] for name in requested]
    else:
        selected = sorted(raw_records, key=lambda record: _natural_key(str(record["fileName"])))
    return raw_records, selected


def write_outputs(
    accession: str,
    records: list[dict[str, Any]],
    requested: list[str],
    manifest: Path,
    meta: Path,
    inventory: Path,
    api_base: str,
) -> None:
    raw_records, selected = select_raw_files(records, requested)
    selected_names = {str(record["fileName"]).strip() for record in selected}

    manifest.parent.mkdir(parents=True, exist_ok=True)
    meta.parent.mkdir(parents=True, exist_ok=True)
    inventory.parent.mkdir(parents=True, exist_ok=True)

    tmp_manifest = manifest.with_suffix(manifest.suffix + ".tmp")
    with tmp_manifest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=MANIFEST_FIELDS, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for task_id, record in enumerate(selected, start=1):
            name = str(record["fileName"]).strip()
            writer.writerow(
                {
                    "task_id": task_id,
                    "pxd_accession": accession,
                    "sdrf_basename": "",
                    "sdrf_path": "",
                    "sdrf_sha256": "",
                    "source_sdrf_url": "",
                    "sdrf_template": "",
                    "sdrf_data_file": "",
                    "archive_file": name,
                    "file_uri": _public_location(record, "FTP Protocol"),
                    "needs_file_map": "0",
                }
            )
    os.replace(tmp_manifest, manifest)

    with inventory.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=INVENTORY_FIELDS, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for record in sorted(records, key=lambda item: _natural_key(str(item["fileName"]))):
            name = str(record["fileName"]).strip()
            writer.writerow(
                {
                    "file_name": name,
                    "file_category": _category(record),
                    "file_size_bytes": _file_size(record),
                    "ftp_uri": _public_location(record, "FTP Protocol"),
                    "aspera_uri": _public_location(record, "Aspera Protocol"),
                    "selected": "1" if name in selected_names else "0",
                }
            )

    raw_bytes = sum(_file_size(record) for record in raw_records)
    selected_bytes = sum(_file_size(record) for record in selected)
    with meta.open("w", encoding="utf-8") as handle:
        handle.write("manifest_mode=pride-project-metadata\n")
        handle.write(f"pxd_accession={accession}\n")
        handle.write(f"pride_api={api_base.rstrip('/')}\n")
        handle.write(f"project_file_count={len(records)}\n")
        handle.write(f"raw_file_count={len(raw_records)}\n")
        handle.write(f"selected_file_count={len(selected)}\n")
        handle.write(f"raw_bytes_total={raw_bytes}\n")
        handle.write(f"selected_bytes_total={selected_bytes}\n")
        handle.write(f"inventory={inventory.resolve()}\n")
        handle.write("trusted_sdrf_used=0\n")
        handle.write("\n[selected_files]\n")
        for record in selected:
            handle.write(str(record["fileName"]).strip() + "\n")


def main() -> int:
    args = parse_args()
    accession = args.accession.strip().upper()
    if not re.fullmatch(r"PXD\d{6,}", accession):
        raise ValueError(f"Invalid PRIDE accession: {args.accession!r}")
    inventory = args.inventory or args.manifest.with_name(
        args.manifest.stem + ".inventory.tsv"
    )
    records = pride_project_files(accession, args.pride_api)
    write_outputs(
        accession,
        records,
        args.files,
        args.manifest,
        args.meta,
        inventory,
        args.pride_api,
    )
    raw_records, selected = select_raw_files(records, args.files)
    print(f"pxd_accession={accession}")
    print(f"project_files={len(records)}")
    print(f"raw_files={len(raw_records)}")
    print(f"selected_files={len(selected)}")
    print(f"manifest={args.manifest}")
    print(f"inventory={inventory}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (
        OSError,
        RuntimeError,
        ValueError,
        csv.Error,
        json.JSONDecodeError,
        urllib.error.HTTPError,
        urllib.error.URLError,
    ) as exc:
        print(f"prepare_prideqc_accession_manifest: {exc}", file=sys.stderr)
        raise SystemExit(2) from exc
