#!/usr/bin/env python3
"""Inspect one failed PRIDE RAW's RAW-category metadata without downloading bytes.

WARNING: PRIDE PEAK-category mzML is not necessarily returned by this API;
absence here DOES NOT establish that a matching open-format file is absent.
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

ACC = re.compile(r"PXD[0-9]{6,}$")


def summarize(records: list[dict], source_name: str) -> dict:
    name = source_name.casefold()
    matches = [r for r in records if isinstance(r, dict) and
               str(r.get("fileName", "")).casefold() == name]
    if len(matches) != 1:
        return {"status": "source_missing_or_ambiguous", "raw_category_matches": len(matches),
                "open_format_inventory_complete": False}
    raw = matches[0]
    size = next((raw[k] for k in ("fileSizeBytes", "fileSizeInBytes", "sizeBytes", "fileSize", "size")
                 if k in raw), None)
    checksum = {k: raw.get(k) for k in ("checksum", "checksumType", "fileChecksum", "md5", "sha256")
                if k in raw and raw[k] is not None}
    base = name[:-4] if name.endswith(".raw") else Path(name).stem
    siblings = [str(r.get("fileName", "")) for r in records if isinstance(r, dict) and
                str(r.get("fileName", "")).casefold().startswith(base) and
                str(r.get("fileName", "")).casefold().endswith((".mzml", ".mzml.gz"))]
    return {"status": "metadata_found", "source_filename": raw["fileName"],
            "published_size_bytes": size, "published_checksums": checksum,
            "raw_category_mzml_same_stem_candidates": sorted(siblings),
            "open_format_inventory_complete": False,
            "reader_error_explanation": "undetermined_by_metadata"}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--accession", required=True)
    ap.add_argument("--source", required=True)
    args = ap.parse_args()
    if not ACC.fullmatch(args.accession) or "/" in args.source or "\\" in args.source:
        ap.error("Only uppercase PXD and a simple filename supported")
    from pridepy.download.client import Client
    records = Client().get_all_raw_file_list(args.accession)
    if not isinstance(records, list):
        raise ValueError("Unexpected PRIDE RAW inventory response")
    print(json.dumps({"accession": args.accession, **summarize(records, args.source)},
                     sort_keys=True, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
