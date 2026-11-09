#!/usr/bin/env python3
"""Upload prideQC canonical runs (export runs.ndjson) to pride-mongo-ws.

Default is a dry run: it maps and validates every document and makes no network call.
A real upload needs --execute and the API key in the environment variable named by
--api-key-env (default PRIDE_MONGO_WS_API_KEY). The key is never printed.

The document id is sha256("PXD|source_file"), so a repeated upload replaces documents
and never makes duplicates.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import sys
import time
import urllib.error
import urllib.request
from collections.abc import Iterator
from pathlib import Path

SHA256 = re.compile(r"^[0-9a-f]{64}$")
FIELDS = {
    "pxd_accession": "projectAccession",
    "source_file": "sourceFile",
    "source_sha256": "sourceSha256",
    "source_size_bytes": "sourceSizeBytes",
    "mzqc_sha256": "mzqcSha256",
    "mzqc_path": "mzqcPath",
    "mzqc_creation_date": "mzqcCreationDate",
    "prideqc_version": "prideqcVersion",
    "pyopenms_version": "pyopenmsVersion",
    "run_tag": "runTag",
    "scientific_approval": "scientificApproval",
    "schema_version": "schemaVersion",
}
MAX_DOC_BYTES = 12 * 1024 * 1024
DEFAULT_URL = "https://pride-mongo-ws-hx-webadmin-11.wp-k8s.ebi.ac.uk/pride/mongo/ws"


def run_id(accession: str, source_file: str) -> str:
    return hashlib.sha256(f"{accession}|{source_file}".encode()).hexdigest()


def to_document(run: dict) -> dict:
    """Map one runs.ndjson record to a MongoPrideQcRun JSON document."""
    doc: dict = {"id": run["_id"]}
    for src, dst in FIELDS.items():
        if src in run:
            doc[dst] = run[src]
    doc["metrics"] = [{"accession": m.get("accession"), "name": m.get("name"),
                       "value": m.get("value"), "unit": m.get("unit")} for m in run["metrics"]]
    expected = run_id(doc["projectAccession"], doc["sourceFile"])
    if doc["id"] != expected:
        raise ValueError(f"id mismatch for {doc['projectAccession']}/{doc['sourceFile']}")
    if not SHA256.match(doc.get("mzqcSha256") or ""):
        raise ValueError(f"bad mzqcSha256 for {doc['projectAccession']}/{doc['sourceFile']}")
    if not doc["metrics"]:
        raise ValueError(f"no metrics for {doc['projectAccession']}/{doc['sourceFile']}")
    if doc.get("scientificApproval") is not False:
        raise ValueError("scientificApproval must be false until scientific validation")
    if len(json.dumps(doc, allow_nan=False).encode()) > MAX_DOC_BYTES:
        raise ValueError(f"document too large: {doc['projectAccession']}/{doc['sourceFile']}")
    return doc


def iter_documents(path: Path, accessions: set | None) -> Iterator[dict]:
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            run = json.loads(line)
            if accessions and run["pxd_accession"] not in accessions:
                continue
            yield to_document(run)


def batches(items: Iterator[dict], size: int) -> Iterator[list[dict]]:
    batch: list[dict] = []
    for item in items:
        batch.append(item)
        if len(batch) == size:
            yield batch
            batch = []
    if batch:
        yield batch


def post(url: str, header: str, key: str, payload: list, retries: int = 3) -> list:
    body = json.dumps(payload, allow_nan=False).encode()
    for attempt in range(1, retries + 1):
        req = urllib.request.Request(url, data=body, method="POST",
                                     headers={"Content-Type": "application/json", header: key})
        try:
            with urllib.request.urlopen(req, timeout=120) as resp:
                return json.loads(resp.read() or b"[]")
        except urllib.error.HTTPError as exc:
            if exc.code < 500 or attempt == retries:
                raise RuntimeError(f"HTTP {exc.code} from {url}: {exc.read()[:300]!r}") from None
        except urllib.error.URLError:
            if attempt == retries:
                raise
        time.sleep(5 * attempt)
    raise AssertionError("unreachable")


def main(argv: list | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs_ndjson", type=Path)
    ap.add_argument("--base-url", default=DEFAULT_URL)
    ap.add_argument("--api-key-header", default="mongo-ws-api-key")
    ap.add_argument("--api-key-env", default="PRIDE_MONGO_WS_API_KEY")
    ap.add_argument("--accession", action="append", default=[], help="only this PXD (repeatable)")
    ap.add_argument("--batch-size", type=int, default=10)
    ap.add_argument("--execute", action="store_true", help="really write to pride-mongo-ws")
    args = ap.parse_args(argv)
    if not 1 <= args.batch_size <= 20:
        ap.error("--batch-size must be 1..20 (server limit)")
    accessions = set(args.accession) or None
    key = os.environ.get(args.api_key_env, "")
    if args.execute and not key:
        ap.error(f"--execute needs the API key in ${args.api_key_env}")
    url = args.base_url.rstrip("/") + "/qc-runs"
    docs = runs = 0
    for batch in batches(iter_documents(args.runs_ndjson, accessions), args.batch_size):
        runs += len(batch)
        if args.execute:
            saved = post(url, args.api_key_header, key, batch)
            ids = {d["id"] for d in batch}
            if {s.get("id") for s in saved} != ids:
                raise RuntimeError(f"server saved {len(saved)} of {len(batch)} documents")
            docs += len(saved)
            print(json.dumps({"event": "BATCH", "saved": docs}), flush=True)
    print(json.dumps({"event": "DONE", "mode": "execute" if args.execute else "dry_run",
                      "validated": runs, "saved": docs, "url": url}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
