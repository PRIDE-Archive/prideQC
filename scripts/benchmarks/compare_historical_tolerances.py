#!/usr/bin/env python3
"""Compare RAW-derived fragment-tolerance candidates with historical SDRF values.

v21 is a validation/stabilization layer. It does not tune the v19 estimator.
Primary comparisons are same-unit only. For ppm/Da mismatches, an optional
mass-context conversion is reported using observed MS2 m/z statistics rather
than being used to redefine the primary comparison.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
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
        help="Optional curated GT TSV. Used as a fallback historical source.",
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
            if value is not None and math.isfinite(value) and value > 0:
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
        if "data file" in name.lower():
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
    match = TOLERANCE_RE.search(value.strip())
    if not match:
        return None, ""
    parsed = float(match.group("value"))
    unit = "ppm" if match.group("unit").lower() == "ppm" else "Da"
    if not math.isfinite(parsed) or parsed <= 0:
        return None, ""
    return parsed, unit


def historical_from_sdrf(
    sdrf_path: Path, raw_file: str
) -> tuple[float | None, str, str, str]:
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
    return {(row["pxd_accession"], row["raw_file"]): row for row in rows}


def read_metrics(path: Path) -> dict[str, str]:
    result: dict[str, str] = {}
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if len(row) >= 3:
                result[row[1]] = row[2]
    return result


def observed_ms2_mass_context(
    metrics: dict[str, str],
) -> tuple[float | None, float | None, float | None]:
    raw = metrics.get("ObservedMzRange_MS2")
    if not raw:
        return None, None, None
    try:
        values = json.loads(raw)
        low, high = float(values[0]), float(values[1])
    except (TypeError, ValueError, IndexError, json.JSONDecodeError):
        return None, None, None
    if not (math.isfinite(low) and math.isfinite(high)) or low <= 0 or high <= 0 or high < low:
        return None, None, None
    midpoint = (low + high) / 2.0
    return low, midpoint, high


def da_to_ppm(tolerance_da: float, mz: float) -> float:
    return tolerance_da / mz * 1_000_000.0


def ppm_to_da(tolerance_ppm: float, mz: float) -> float:
    return tolerance_ppm * mz / 1_000_000.0


def mass_context_conversion(
    historical_value: float | None,
    historical_unit: str,
    candidate_value: float | None,
    candidate_unit: str,
    mz_range: tuple[float | None, float | None, float | None],
) -> dict[str, str]:
    low, mid, high = mz_range
    if historical_value is None or candidate_value is None or not low or not mid or not high:
        return {
            "mass_context_unit": "",
            "mass_context_median_mz": "",
            "historical_equivalent_median_unit": "",
            "historical_equivalent_median": "",
            "candidate_equivalent_median_unit": "",
            "candidate_equivalent_median": "",
            "mass_context_ratio_candidate_over_historical": "",
        }

    if historical_unit == candidate_unit:
        return {
            "mass_context_unit": historical_unit,
            "mass_context_median_mz": f"{mid:.10g}",
            "historical_equivalent_median_unit": historical_unit,
            "historical_equivalent_median": f"{historical_value:.10g}",
            "candidate_equivalent_median_unit": candidate_unit,
            "candidate_equivalent_median": f"{candidate_value:.10g}",
            "mass_context_ratio_candidate_over_historical": (
                f"{candidate_value / historical_value:.10g}"
                if historical_value
                else ""
            ),
        }

    if historical_unit == "Da" and candidate_unit == "ppm":
        historical_eq = da_to_ppm(historical_value, mid)
        candidate_eq = candidate_value
        target_unit = "ppm"
    elif historical_unit == "ppm" and candidate_unit == "Da":
        historical_eq = historical_value
        candidate_eq = ppm_to_da(candidate_value, mid)
        target_unit = "Da"
    else:
        return {k: "" for k in (
            "mass_context_unit", "mass_context_median_mz",
            "historical_equivalent_median_unit", "historical_equivalent_median",
            "candidate_equivalent_median_unit", "candidate_equivalent_median",
            "mass_context_ratio_candidate_over_historical")}

    return {
        "mass_context_unit": target_unit,
        "mass_context_median_mz": f"{mid:.10g}",
        "historical_equivalent_median_unit": target_unit,
        "historical_equivalent_median": f"{historical_eq:.10g}",
        "candidate_equivalent_median_unit": target_unit,
        "candidate_equivalent_median": f"{candidate_eq:.10g}",
        "mass_context_ratio_candidate_over_historical": (
            f"{candidate_eq / historical_eq:.10g}"
            if historical_eq
            else ""
        ),
    }


def build_row(result_dir: Path, gt: dict[tuple[str, str], dict[str, str]]) -> dict[str, str]:
    info = {}
    for line in (result_dir / "hpc-run-info.txt").read_text(encoding="utf-8").splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            info[key] = value

    raw_file = info.get("archive_file", "")
    pxd = info.get("pxd", "")
    annotations = read_annotations(result_dir / "annotations.tsv")
    metrics = read_metrics(result_dir / "metrics.tsv")
    candidate_unit, candidate_value, candidate_payload = candidate_from_annotations(annotations)

    historical_value, historical_unit, historical_raw, historical_source = historical_from_sdrf(
        result_dir / "input.sdrf.tsv", raw_file
    )
    if historical_value is None:
        gt_row = gt.get((pxd, raw_file))
        if gt_row:
            raw_value = gt_row.get("fragment_tolerance_value", "")
            unit = gt_row.get("fragment_tolerance_unit", "")
            try:
                historical_value = float(raw_value) if raw_value else None
            except ValueError:
                historical_value = None
            historical_unit = unit if unit in {"ppm", "Da"} else ""
            historical_raw = gt_row.get("fragment_tolerance_raw", raw_value)
            historical_source = "ground-truth"

    if not candidate_unit:
        comparison = "candidate-unavailable"
    elif historical_value is None or not historical_unit:
        comparison = "historical-unavailable"
    elif candidate_unit != historical_unit:
        comparison = "unit-incompatible"
    else:
        comparison = "same-unit"

    ratio = difference = percent = None
    if comparison == "same-unit" and historical_value:
        ratio = candidate_value / historical_value
        difference = candidate_value - historical_value
        percent = difference / historical_value * 100.0

    context = mass_context_conversion(
        historical_value,
        historical_unit,
        candidate_value,
        candidate_unit,
        observed_ms2_mass_context(metrics),
    )

    precision_field = (
        "estimated_fragment_mass_error_ppm"
        if candidate_unit == "ppm"
        else "estimated_fragment_mass_error_da"
    )
    precision = parse_json_value(annotations.get(precision_field))
    inlier_fraction = (
        float(precision["robust_inlier_fraction"])
        if isinstance(precision, dict)
        and precision.get("robust_inlier_fraction") is not None
        else None
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
        "candidate_confidence": (
            str(candidate_payload.get("confidence", ""))
            if candidate_payload
            else ""
        ),
        "resolution_regime": (
            str(candidate_payload.get("resolution_regime", ""))
            if candidate_payload
            else ""
        ),
        "fragment_match_window_da": fmt(
            float(candidate_payload["fragment_match_window_da"])
            if candidate_payload
            and candidate_payload.get("fragment_match_window_da") is not None
            else None
        ),
        "window_censored": (
            str(candidate_payload.get("window_censored", ""))
            if candidate_payload
            else ""
        ),
        "robust_inlier_fraction": fmt(inlier_fraction),
        **context,
    }


def fmt(value: float | None) -> str:
    return "" if value is None else f"{value:.10g}"


def main() -> None:
    args = parse_args()
    gt = load_gt(args.gt)
    result_dirs = sorted(
        path.parent
        for path in args.run_root.rglob("hpc-run-info.txt")
        if (path.parent / "annotations.tsv").exists() and (path.parent / "metrics.tsv").exists()
    )
    if not result_dirs:
        raise SystemExit(
            "No result directories with annotations and metrics found under "
            f"{args.run_root}"
        )

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
        ratios = [
            float(r["ratio_candidate_over_historical"])
            for r in same
            if r["ratio_candidate_over_historical"]
        ]
        pcts = [float(r["percent_difference"]) for r in same if r["percent_difference"]]
        context_ratios = [
            float(r["mass_context_ratio_candidate_over_historical"])
            for r in group
            if r["mass_context_ratio_candidate_over_historical"]
        ]
        candidate_values = [float(r["candidate_value"]) for r in group if r["candidate_value"]]
        historical_values = [float(r["historical_value"]) for r in group if r["historical_value"]]
        inferred = [r for r in group if r["candidate_value"]]
        summaries.append({
            "pxd": pxd,
            "files": str(len(group)),
            "candidate_inferred": str(len(inferred)),
            "candidate_unavailable": str(len(group) - len(inferred)),
            "historical_available": str(sum(bool(r["historical_value"]) for r in group)),
            "same_unit_files": str(len(same)),
            "unit_incompatible_files": str(
                sum(r["comparison"] == "unit-incompatible" for r in group)
            ),
            "candidate_unavailable_files": str(
                sum(r["comparison"] == "candidate-unavailable" for r in group)
            ),
            "historical_unavailable_files": str(
                sum(r["comparison"] == "historical-unavailable" for r in group)
            ),
            "median_historical": fmt(median(historical_values) if historical_values else None),
            "median_candidate": fmt(median(candidate_values) if candidate_values else None),
            "median_ratio_candidate_over_historical": fmt(median(ratios) if ratios else None),
            "min_ratio_candidate_over_historical": fmt(min(ratios) if ratios else None),
            "max_ratio_candidate_over_historical": fmt(max(ratios) if ratios else None),
            "median_percent_difference": fmt(median(pcts) if pcts else None),
            "median_mass_context_ratio": fmt(median(context_ratios) if context_ratios else None),
            "candidate_units": ",".join(
                f"{u}:{sum(r['candidate_unit'] == u for r in group)}"
                for u in ("ppm", "Da")
                if any(r["candidate_unit"] == u for r in group)
            ),
            "historical_units": ",".join(
                f"{u}:{sum(r['historical_unit'] == u for r in group)}"
                for u in ("ppm", "Da")
                if any(r["historical_unit"] == u for r in group)
            ),
            "windows_da": ",".join(
                f"{w}:{sum(r['fragment_match_window_da'] == w for r in group)}"
                for w in ("0.2", "0.5", "1")
                if any(r["fragment_match_window_da"] == w for r in group)
            ),
            "confidence": ",".join(
                f"{c}:{sum(r['candidate_confidence'] == c for r in group)}"
                for c in ("high", "moderate", "low")
                if any(r["candidate_confidence"] == c for r in group)
            ),
        })

    with args.accessions_output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summaries[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(summaries)

    counts = Counter(row["comparison"] for row in rows)
    print(f"files={len(rows)}")
    print(f"accessions={len(summaries)}")
    print("comparison_counts=" + ",".join(f"{key}:{counts[key]}" for key in sorted(counts)))
    print(f"files_output={args.files_output}")
    print(f"accessions_output={args.accessions_output}")


if __name__ == "__main__":
    main()
