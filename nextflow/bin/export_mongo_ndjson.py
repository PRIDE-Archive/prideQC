#!/usr/bin/env python3
"""Prepare idempotent MongoDB documents WITHOUT connecting to a database.

One metric document per run/metric, plus an accession/revision manifest.
Never silently truncates metrics. Validate output with user before DB integration.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

LIMIT = 12 * 1024 * 1024  # Conservative BSON safety margin; no metric truncation


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def export(mzqc_file: Path, outdir: Path) -> dict:
    sidecar = mzqc_file.with_suffix(".manifest.json")
    if not sidecar.is_file():
        raise ValueError("Missing accession merger manifest")
    manifest = json.loads(sidecar.read_text())
    if manifest.get("status") != "complete":
        raise ValueError("Refusing incomplete mzQC artifact")
    digest = sha256(mzqc_file)
    if digest != manifest.get("mzqc_sha256"):
        raise ValueError("mzQC SHA256 mismatch")
    obo = mzqc_file.with_suffix(".obo")
    if sha256(obo) != manifest.get("obo_sha256"):
        raise ValueError("OBO SHA256 mismatch")
    root = json.loads(mzqc_file.read_text())["mzQC"]
    pxd = manifest["accession"]
    if (
        len(root["runQualities"]) != manifest["merged_runs"]
        or manifest["expected_runs"] != manifest["merged_runs"]
    ):
        raise ValueError("Run-count mismatch")
    outdir.mkdir(parents=True, exist_ok=True)
    target = outdir / f"{pxd}.metrics.ndjson"
    temp = target.with_suffix(".ndjson.tmp")
    n = 0
    labels: set[str] = set()
    try:
        with temp.open("w", encoding="utf-8") as handle:
            for run in root["runQualities"]:
                label = run["metadata"]["label"]
                if label in labels:
                    raise ValueError(f"Duplicate run label: {label}")
                labels.add(label)
                input_names = [f["name"] for f in run["metadata"]["inputFiles"]]
                for idx, metric in enumerate(run["qualityMetrics"]):
                    record = {
                        "_id": hashlib.sha256(f"{pxd}|{digest}|{label}|{idx}".encode()).hexdigest(),
                        "schema_version": 1,
                        "pxd_accession": pxd,
                        "mzqc_sha256": digest,
                        "run_label": label,
                        "input_files": input_names,
                        "metric_index": idx,
                        "metric_accession": metric["accession"],
                        "metric_name": metric["name"],
                        "metric_value": metric["value"],
                        "metric_description": metric.get("description"),
                    }
                    data = json.dumps(record, ensure_ascii=False, allow_nan=False)
                    if len(data.encode("utf-8")) > LIMIT:
                        raise ValueError(
                            f"Oversized metric; store in external object storage: {label}/{idx}"
                        )
                    handle.write(data + "\n")
                    n += 1
        temp.replace(target)
    except Exception:
        temp.unlink(missing_ok=True)
        raise
    receipt = {
        "schema_version": 1,
        "pxd_accession": pxd,
        "mzqc_sha256": digest,
        "source_manifest_sha256": sha256(sidecar),
        "observed_runs": len(labels),
        "expected_runs": manifest["expected_runs"],
        "metric_documents": n,
        "metrics_ndjson_sha256": sha256(target),
        "mongo_write_performed": False,
    }
    (outdir / f"{pxd}.ingestion-manifest.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n"
    )
    return receipt


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--mzqc", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    args = ap.parse_args()
    print(json.dumps(export(args.mzqc, args.outdir), sort_keys=True))


if __name__ == "__main__":
    main()
