#!/usr/bin/env python3
"""Read-only mzQC evidence inventory for scientific qualification planning.

This is NOT a scientific correctness test. Presence of a metric does not
validate its estimate, ontology semantics or experimental calibration.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

GROUP_TOKENS = {
    "mass_error": ("mass error", "search tolerance", "mass_error"),
    "mass_shift_ptm": ("mass shift", "modification mass shift", "mass_shift"),
    "tmt": ("TMT reporter", "TMT experimental", "tmt_reporter", "tmt_profile"),
}


def audit(run_dir: Path) -> dict:
    items = []
    for receipt_path in sorted(run_dir.rglob("result.json")):
        receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
        pxd = receipt.get("accession")
        entry = {"accession": pxd, "unit_id": receipt.get("unit_id"),
                 "receipt_status": receipt.get("status"),
                 "requested_estimators": receipt.get("requested_estimators", "unknown_old_receipt"),
                 "metric_count": 0, "metric_groups_present": [], "mzqc_sha256_verified": False,
                 "scientific_validity_verified": False}
        if receipt.get("status") == "success":
            qc = receipt_path.parent / "qc" / str(receipt.get("mzqc_file", ""))
            if not qc.is_file():
                raise ValueError(f"Missing mzQC for {pxd}")
            if hashlib.sha256(qc.read_bytes()).hexdigest() != receipt.get("mzqc_sha256"):
                raise ValueError(f"QC digest mismatch for {pxd}")
            data = json.loads(qc.read_text(encoding="utf-8"))
            body = data.get("mzQC", {})
            runs = body.get("runQualities", [])
            if body.get("version") != "1.0.0" or not isinstance(runs, list) or len(runs) != 1:
                raise ValueError(f"Unexpected per-run mzQC shape for {pxd}")
            metrics = runs[0].get("qualityMetrics", [])
            if not isinstance(metrics, list) or not all(
                isinstance(m, dict) and isinstance(m.get("name"), str) and
                isinstance(m.get("accession"), str) and "value" in m for m in metrics
            ):
                raise ValueError(f"Malformed qualityMetrics for {pxd}")
            entry["metric_count"] = len(metrics)
            names = "\n".join(str(m.get("name", "")) + " " + str(m.get("description", ""))
                              for m in metrics).casefold()
            entry["metric_groups_present"] = sorted(group for group, phrases in GROUP_TOKENS.items()
                                                    if any(phrase.casefold() in names for phrase in phrases))
            entry["mzqc_sha256_verified"] = True
        items.append(entry)
    if not items:
        raise ValueError("No QC receipts in this run directory")
    return {"schema_version": 1, "purpose": "evidence_inventory_not_scientific_qualification",
            "units": len(items), "successful_receipts": sum(x["receipt_status"] == "success" for x in items),
            "scientific_validity_verified": False, "records": items}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(audit(args.run), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
