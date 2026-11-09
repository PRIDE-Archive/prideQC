#!/usr/bin/env python3
"""Per-project prideQC summaries for the PRIDE Archive search filters.

Reads an export runs.ndjson (one canonical run per line, sorted by project) and
writes one JSON object per project (NDJSON). Each summary holds:

- qcRuns: number of runs with QC;
- acquisitionMethods, ms2Resolutions: the determined values over the runs
  ("Not determined" only when no run has a value);
- medians of run duration (min), MS1 spectra and MS2 spectra per run;
- one bucket label per median. The search page filters on these labels.

The labels are search filter values: they contain no ",", "==", "*" or "?",
because the search API splits filters on "," and "==" and treats "*" and "?"
as wildcards.

Default is a dry run: compute and write --output, no network call. --execute
also uploads the summaries to pride-mongo-ws POST /qc-project/batch (API key
from the environment variable named by --api-key-env; the key is never printed).
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import importlib.util
import json
import os
import statistics
import sys
import time
import urllib.error
import urllib.request
from collections.abc import Iterator
from pathlib import Path
from typing import Any

UTC = dt.timezone.utc  # noqa: UP017 (Codon host runs Python 3.9)
SCHEMA_VERSION = 1
DEFAULT_URL = "https://pride-mongo-ws-hx-webadmin-11.wp-k8s.ebi.ac.uk/pride/mongo/ws"
BATCH_MAX = 200
NOT_DETERMINED = "Not determined"

# qc_archive_stats.py owns the metric parsing (accessions, units, categories): reuse it
_SPEC = importlib.util.spec_from_file_location(
    "prideqc_qc_archive_stats", Path(__file__).resolve().parent / "qc_archive_stats.py")
assert _SPEC is not None and _SPEC.loader is not None
QS = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(QS)

# (upper bound, label); the last bucket has no upper bound. Bounds follow the
# v8 project medians (duration p05-p95 34-236 min, MS1 0.9k-20k, MS2 2.9k-119k).
BUCKETS: dict[str, list[tuple[float | None, str]]] = {
    "run_duration_min": [
        (30, "under 30 min"), (60, "30 to 60 min"), (120, "60 to 120 min"),
        (180, "120 to 180 min"), (None, "180 min or more"),
    ],
    "ms1_spectra": [
        (1000, "under 1k"), (5000, "1k to 5k"), (10000, "5k to 10k"),
        (20000, "10k to 20k"), (None, "20k or more"),
    ],
    "ms2_spectra": [
        (10000, "under 10k"), (30000, "10k to 30k"), (60000, "30k to 60k"),
        (120000, "60k to 120k"), (None, "120k or more"),
    ],
}


def bucket(kind: str, value: float | None) -> str | None:
    if value is None:
        return None
    for upper, label in BUCKETS[kind]:
        if upper is None or value < upper:
            return label
    return None  # pragma: no cover (last bucket has no bound)


def _determined(values: set[str]) -> list[str]:
    known = sorted(v for v in values if v != NOT_DETERMINED)
    return known or [NOT_DETERMINED]


def _summary(accession: str, runs: int, values: dict[str, list[float]], cats: dict[str, set[str]]) -> dict:
    def median(kind: str) -> float | None:
        vals = values.get(kind) or []
        return round(statistics.median(vals), 3) if vals else None

    duration, ms1, ms2 = median("run_duration_min"), median("ms1_spectra"), median("ms2_spectra")
    return {
        "id": accession,
        "projectAccession": accession,
        "qcRuns": runs,
        "acquisitionMethods": _determined(cats["acquisition_method"]),
        "ms2Resolutions": _determined(cats["ms2_resolution"]),
        "medianRunDurationMin": duration,
        "medianMs1Spectra": ms1,
        "medianMs2Spectra": ms2,
        "runDurationBucket": bucket("run_duration_min", duration),
        "ms1SpectraBucket": bucket("ms1_spectra", ms1),
        "ms2SpectraBucket": bucket("ms2_spectra", ms2),
    }


def summaries(path: Path) -> Iterator[dict]:
    """Stream one summary per project; runs.ndjson must be sorted by project."""
    current: str | None = None
    done: set[str] = set()
    runs = 0
    values: dict[str, list[float]] = {}
    cats: dict[str, set[str]] = {}

    def reset() -> None:
        nonlocal runs, values, cats
        runs = 0
        values = {k: [] for k in BUCKETS}
        cats = {"acquisition_method": set(), "ms2_resolution": set()}

    reset()
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            run = json.loads(line)
            accession = str(run.get("pxd_accession") or "")
            if not accession:
                raise ValueError("run without pxd_accession: " + str(run.get("source_file")))
            if accession != current:
                if accession in done:
                    raise ValueError("runs.ndjson must be sorted by project; " + accession + " appears twice")
                if current is not None:
                    yield _summary(current, runs, values, cats)
                    done.add(current)
                reset()
                current = accession
            nums, run_cats, _ = QS.run_values(run)
            runs += 1
            for k in BUCKETS:
                if k in nums:
                    values[k].append(nums[k])
            for k in cats:
                cats[k].add(run_cats[k])
    if current is not None:
        yield _summary(current, runs, values, cats)


def stamp(doc: dict, generated_at: str, source: str, source_sha256: str) -> dict:
    return {**doc, "schemaVersion": SCHEMA_VERSION, "generatedAt": generated_at, "source": source,
            "sourceSha256": source_sha256, "scientificApproval": False}


def post(url: str, header: str, key: str, payload: list, attempts: int = 4) -> Any:
    body = json.dumps(payload, allow_nan=False).encode()
    for attempt in range(1, attempts + 1):
        req = urllib.request.Request(url, data=body, method="POST",
                                     headers={"Content-Type": "application/json", header: key})
        try:
            with urllib.request.urlopen(req, timeout=120) as resp:
                return json.loads(resp.read() or b"[]")
        except urllib.error.HTTPError as exc:
            if exc.code < 500 or attempt == attempts:
                raise RuntimeError(f"HTTP {exc.code} from {url}: {exc.read()[:300]!r}") from None
        except (urllib.error.URLError, TimeoutError):
            if attempt == attempts:
                raise
        time.sleep(2 ** attempt)
    return None  # pragma: no cover


def main(argv: list | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs_ndjson", type=Path)
    ap.add_argument("--output", type=Path, required=True, help="where to write the summaries (NDJSON)")
    ap.add_argument("--base-url", default=DEFAULT_URL)
    ap.add_argument("--api-key-header", default="mongo-ws-api-key")
    ap.add_argument("--api-key-env", default="PRIDE_MONGO_WS_API_KEY")
    ap.add_argument("--batch-size", type=int, default=100)
    ap.add_argument("--execute", action="store_true", help="also upload to pride-mongo-ws /qc-project/batch")
    args = ap.parse_args(argv)
    if not 1 <= args.batch_size <= BATCH_MAX:
        ap.error(f"--batch-size must be 1 to {BATCH_MAX}")
    key = os.environ.get(args.api_key_env, "")
    if args.execute and not key:
        ap.error(f"--execute needs the API key in ${args.api_key_env}")

    digest = hashlib.sha256()
    with args.runs_ndjson.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    generated = dt.datetime.now(UTC).replace(microsecond=0).isoformat().replace("+00:00", "Z")
    source = args.runs_ndjson.parent.name + "/" + args.runs_ndjson.name

    docs = [stamp(d, generated, source, digest.hexdigest()) for d in summaries(args.runs_ndjson)]
    tmp = args.output.with_suffix(args.output.suffix + ".tmp")
    with tmp.open("w", encoding="utf-8") as out:
        for d in docs:
            out.write(json.dumps(d, allow_nan=False, separators=(",", ":")) + "\n")
    tmp.replace(args.output)

    url = args.base_url.rstrip("/") + "/qc-project/batch"
    saved = 0
    if args.execute:
        for i in range(0, len(docs), args.batch_size):
            batch = docs[i:i + args.batch_size]
            result = post(url, args.api_key_header, key, batch)
            if not isinstance(result, list) or len(result) != len(batch):
                raise RuntimeError("server did not confirm batch at offset " + str(i))
            saved += len(batch)
    print(json.dumps({"event": "DONE", "mode": "execute" if args.execute else "dry_run",
                      "projects": len(docs), "runs": sum(d["qcRuns"] for d in docs), "saved": saved,
                      "output": str(args.output), "url": url}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
