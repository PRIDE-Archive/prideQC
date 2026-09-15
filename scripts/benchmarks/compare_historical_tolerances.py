#!/usr/bin/env python3
"""Compare RAW-derived v19 fragment tolerance candidates with SDRF annotations.

The script does not tune the estimator. It reports whether a candidate and a
historical SDRF fragment-tolerance annotation are directly comparable. Ratios
and differences are calculated only for matching units.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from collections import Counter, defaultdict
from pathlib import Path
from statistics import median

TOLERANCE_RE = re.compile(
    r"(?P<value>[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s*(?P<unit>ppm|Da)\b",
    re.IGNORECASE,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root", required=True, type=Path)
    parser.add_argument("--files-output", required=True, type=Path)
    parser.add_argument("--accessions-output", required=True, type=Path)
    parser.add_argument(
        "--gt",
        type=Path,
        help="Optional curated GT TSV. Used as a fallback historical source when a run SDRF cannot be parsed.",
    )
    return parser.parse_args()


def parse_json_value(row: dict[str, str] | None):
    if not row:
        return None
    value = row.get("value", "")
    if not value or value == "null":
        return None
    try:
        return json.loads(value)
    except json.JSONDecodeError:
        return None


def read_annotations(path: Path) -> dict[str, dict[str, str]]:
    result: dict[str, dict[str, str]] = {}
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            result[row["field"]] = row
    return result


def candidate_from_annotations(
    annotations: dict[str, dict[str, str]],
) -> tuple[str, float | None, dict | None]:
    for unit, field in (
        ("ppm", "suggested_fragment_search_tolerance_ppm"),
        ("Da", "suggested_fragment_search_tolerance_da"),
    ):
        row = annotations.get(field)
        payload = parse_json_value(row)
        if row and row.get("evidence") == "inferred" and isinstance(payload, dict):
            value = payload.get("suggested_tolerance")
            try:
                value = float(value)
            except (TypeError, ValueError):
                value = None
            if value is not None:
                return unit, value, payload
    return "", None, None


def detect_data_file_column(fieldnames: list[str]) -> str | None:
    preferred = (
        "comment[data file]",
        "comment[data file name]",
        "assay name",
    )
    for name in preferred:
        if name in fieldnames:
            return name
    for name in fieldnames:
        lowered = name.lower()
        if "data file" in lowered:
            return name
    return None


def detect_fragment_tolerance_column(fieldnames: list[str]) -> str | None:
    preferred = (
        "comment[fragment mass tolerance]",
        "comment[fragment mass tolerance unit]",
    )
    for name in preferred:
        if name in fieldnames:
            return name
    for name in fieldnames:
        lowered = name.lower()
        if "fragment" in lowered and "tolerance" in lowered:
            return name
    return None


def parse_historical_tolerance(value: str) -> tuple[float | None, str]:
    value = value.strip()
    match = TOLERANCE_RE.search(value)
    if not match:
        return None, ""
    try:
        parsed = float(match.group("value"))
    except ValueError:
        return None, ""
    unit = "ppm" if match.group("unit").lower() == "ppm" else "Da"
    return parsed, unit


def historical_from_sdrf(sdrf_path: Path, raw_file: str) -> tuple[float | None, str, str, str]:
    if not sdrf_path.exists():
        return None, "", "", "missing-sdrf"

    with sdrf_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = reader.fieldnames or []
        data_col = detect_data_file_column(fieldnames)
        tolerance_col = detect_fragment_tolerance_column(fieldnames)
        if not data_col or not tolerance_col:
            return None, "", "", "sdrf-column-not-found"
        for row in reader:
            if row.get(data_col, "").strip() != raw_file:
                continue
            raw_value = row.get(tolerance_col, "")
            value, unit = parse_historical_tolerance(raw_value)
            return value, unit, raw_value, "run-sdrf"
    return None, "", "", "sdrf-row-not-found"


def load_gt(path: Path | None) -> dict[tuple[str, str], dict[str, str]]:
    if path is None or not path.exists():
        return {}
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return {
        (row["pxd_accession"], row["raw_file"]): row
        for row in rows
    }


def fmt(value: float | None) -> str:
    return "" if value is None else f"{value:.10g}"


def median_or_blank(values: list[float]) -> float | None:
    return median(values) if values else None


def build_row(
    result_dir: Path,
    gt: dict[tuple[str, str], dict[str, str]],
) -> dict[str, str]:
    info = {}
    info_path = result_dir / "hpc-run-info.txt"
    for line in info_path.read_text(encoding="utf-8").splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            info[key] = value

    raw_file = info.get("archive_file", "")
    pxd = info.get("pxd", "")
    annotations = read_annotations(result_dir / "annotations.tsv")
    candidate_unit, candidate_value, candidate_payload = candidate_from_annotations(annotations)

    historical_value = None
    historical_unit = ""
    historical_raw = ""
    historical_source = ""

    sdrf_path = result_dir / "input.sdrf.tsv"
    historical_value, historical_unit, historical_raw, historical_source = historical_from_sdrf(
        sdrf_path,
        raw_file,
    )

    if historical_value is None:
        gt_row = gt.get((pxd, raw_file))
        if gt_row:
            historical_raw = gt_row.get("fragment_tolerance_value", "")
            gt_unit = gt_row.get("fragment_tolerance_unit", "")
            historical_value = (
                float(historical_raw) if historical_raw else None
            )
            historical_unit = gt_unit if gt_unit in {"ppm", "Da"} else ""
            historical_source = "ground-truth"

    if not candidate_unit:
        comparison = "candidate-unavailable"
    elif historical_value is None or not historical_unit:
        comparison = "historical-unavailable"
    elif candidate_unit != historical_unit:
        comparison = "unit-incompatible"
    else:
        comparison = "same-unit"

    ratio = None
    difference = None
    percent = None
    if comparison == "same-unit" and historical_value:
        ratio = candidate_value / historical_value
        difference = candidate_value - historical_value
        percent = difference / historical_value * 100

    precision = parse_json_value(
        annotations.get(
            "estimated_fragment_mass_error_ppm"
            if candidate_unit == "ppm"
            else "estimated_fragment_mass_error_da"
        )
    )

    return {
        "task_id": info.get("task_id", ""),
        "pxd": pxd,
        "file": raw_file,
        "task_outcome": info.get("task_outcome", ""),
        "historical_source": historical_source,
        "historical_raw": historical_raw,
        "historical_unit": historical_unit,
        "historical_value": fmt(historical_value),
        "candidate_unit": candidate_unit,
        "candidate_value": fmt(candidate_value),
        "comparison": comparison,
        "ratio_candidate_over_historical": fmt(ratio),
        "absolute_difference": fmt(difference),
        "percent_difference": fmt(percent),
        "candidate_confidence": str(candidate_payload.get("confidence", "")) if candidate_payload else "",
        "resolution_regime": str(candidate_payload.get("resolution_regime", "")) if candidate_payload else "",
        "fragment_match_window_da": fmt(
            float(candidate_payload["fragment_match_window_da"])
            if candidate_payload and candidate_payload.get("fragment_match_window_da") is not None
            else None
        ),
        "window_censored": str(candidate_payload.get("window_censored", "")) if candidate_payload else "",
        "robust_inlier_fraction": fmt(
            float(precision["robust_inlier_fraction"])
            if isinstance(precision, dict) and precision.get("robust_inlier_fraction") is not None
            else None
        ),
    }


def main() -> None:
    args = parse_args()
    gt = load_gt(args.gt)

    result_dirs = sorted(
        path.parent for path in args.run_root.rglob("hpc-run-info.txt")
        if (path.parent / "annotations.tsv").exists()
    )
    if not result_dirs:
        raise SystemExit(f"No completed annotation result directories found under {args.run_root}")

    rows = [build_row(path, gt) for path in result_dirs]
    args.files_output.parent.mkdir(parents=True, exist_ok=True)
    args.accessions_output.parent.mkdir(parents=True, exist_ok=True)

    fields = list(rows[0])
    with args.files_output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(sorted(rows, key=lambda r: (r["pxd"], r["file"])))

    by_pxd: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        by_pxd[row["pxd"]].append(row)

    summaries = []
    for pxd, group in sorted(by_pxd.items()):
        same = [r for r in group if r["comparison"] == "same-unit"]
        inferred = [r for r in group if r["candidate_value"]]
        candidates = [float(r["candidate_value"]) for r in inferred]
        historical = [
            float(r["historical_value"])
            for r in same
            if r["historical_value"]
        ]
        ratios = [
            float(r["ratio_candidate_over_historical"])
            for r in same
            if r["ratio_candidate_over_historical"]
        ]
        pcts = [
            float(r["percent_difference"])
            for r in same
            if r["percent_difference"]
        ]
        summaries.append({
            "pxd": pxd,
            "files": str(len(group)),
            "candidate_inferred": str(len(inferred)),
            "candidate_unavailable": str(len(group) - len(inferred)),
            "historical_available": str(sum(bool(r["historical_value"]) for r in group)),
            "same_unit_files": str(len(same)),
            "unit_incompatible_files": str(sum(r["comparison"] == "unit-incompatible" for r in group)),
            "candidate_unavailable_files": str(sum(r["comparison"] == "candidate-unavailable" for r in group)),
            "median_historical": fmt(median_or_blank(historical)),
            "median_candidate": fmt(median_or_blank(candidates)),
            "median_ratio_candidate_over_historical": fmt(median_or_blank(ratios)),
            "min_ratio_candidate_over_historical": fmt(min(ratios) if ratios else None),
            "max_ratio_candidate_over_historical": fmt(max(ratios) if ratios else None),
            "median_percent_difference": fmt(median_or_blank(pcts)),
            "candidate_units": ",".join(
                f"{unit}:{sum(r['candidate_unit'] == unit for r in group)}"
                for unit in ("ppm", "Da")
                if any(r["candidate_unit"] == unit for r in group)
            ),
            "historical_units": ",".join(
                f"{unit}:{sum(r['historical_unit'] == unit for r in group)}"
                for unit in ("ppm", "Da")
                if any(r["historical_unit"] == unit for r in group)
            ),
            "windows_da": ",".join(
                f"{window}:{sum(r['fragment_match_window_da'] == window for r in group)}"
                for window in ("0.2", "0.5", "1")
                if any(r["fragment_match_window_da"] == window for r in group)
            ),
            "confidence": ",".join(
                f"{confidence}:{sum(r['candidate_confidence'] == confidence for r in group)}"
                for confidence in ("high", "moderate", "low")
                if any(r["candidate_confidence"] == confidence for r in group)
            ),
        })

    summary_fields = list(summaries[0])
    with args.accessions_output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=summary_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(summaries)

    counts = Counter(row["comparison"] for row in rows)
    print(f"files={len(rows)}")
    print(f"accessions={len(summaries)}")
    print("comparison_counts=" + ",".join(f"{k}:{counts[k]}" for k in sorted(counts)))
    print(f"files_output={args.files_output}")
    print(f"accessions_output={args.accessions_output}")


if __name__ == "__main__":
    main()
