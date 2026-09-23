"""Cohort-level QC synthesis for conservative SDRF refinement.

Per-run estimators remain the source of measurements.  This module combines
those measurements across compatible runs so that one reanalysis tolerance can
be written for each inferred experiment group and only exceptionally strong,
chemically unambiguous recurrent PTM mass families are eligible for SDRF
annotation.
"""

from __future__ import annotations

import csv
import math
from collections import defaultdict
from dataclasses import dataclass
from math import ceil, comb, log10
from pathlib import Path
from statistics import median
from typing import Any

import numpy as np

from prideqc.models import AnalysisResult, Annotation, CVTerm, EvidenceKind

EXPERIMENT_GROUP_FIELD = "prideqc_experiment_group"
COHORT_PRECURSOR_FIELD = "cohort_precursor_search_tolerance"
COHORT_FRAGMENT_FIELD = "cohort_fragment_search_tolerance"
COHORT_MODIFICATION_FIELD = "cohort_high_support_modification"
COHORT_SDRF_FIELDS = frozenset(
    {COHORT_PRECURSOR_FIELD, COHORT_FRAGMENT_FIELD, COHORT_MODIFICATION_FIELD}
)
COHORT_OVERWRITE_FIELDS = frozenset({COHORT_PRECURSOR_FIELD, COHORT_FRAGMENT_FIELD})

_NUMERIC_FEATURES = ("precursor", "fragment", "rt", "isolation", "ms1", "ms2")
_SEPARATION_THRESHOLDS = {
    "precursor": 1.50,
    "fragment": 1.50,
    "rt": 1.15,
    "isolation": 1.20,
    "ms1": 1.25,
    "ms2": 1.25,
}
_FEATURE_NAMES = {
    "precursor": "precursor tolerance",
    "fragment": "fragment tolerance",
    "rt": "chromatography duration",
    "isolation": "MS2 isolation width",
    "ms1": "MS1 spectra",
    "ms2": "MS2 spectra",
}


@dataclass(frozen=True, slots=True)
class SemanticEvidence:
    status: str
    sources: tuple[str, ...] = ()
    notes: tuple[str, ...] = ()


@dataclass(frozen=True, slots=True)
class CohortSynthesis:
    assignments: dict[str, str]
    groups: dict[str, dict[str, Any]]
    evidence: list[str]
    ptm_families: list[dict[str, Any]]
    ptm_review_families: list[dict[str, Any]]

    def to_dict(self) -> dict[str, Any]:
        return {
            "assignments": self.assignments,
            "groups": self.groups,
            "evidence": self.evidence,
            "ptm_families": self.ptm_families,
            "ptm_review_families": self.ptm_review_families,
        }


def read_semantic_evidence(path: Path) -> dict[tuple[str, str], SemanticEvidence]:
    """Read independent accession/UniMod evidence using the v4 benchmark schema."""
    allowed = {"supported", "conflicting", "not-found"}
    grouped: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"pxd_accession", "unimod_accession", "evidence_status"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(
                "PTM study evidence TSV is missing required columns: "
                + ", ".join(sorted(missing))
            )
        for row in reader:
            accession = (row.get("pxd_accession") or "").strip().upper()
            unimod = (row.get("unimod_accession") or "").strip().upper()
            status = (row.get("evidence_status") or "").strip().casefold()
            if not accession or not unimod:
                raise ValueError(
                    "PTM study evidence rows require pxd_accession and unimod_accession"
                )
            if status not in allowed:
                raise ValueError(
                    "PTM study evidence status must be supported, conflicting, or not-found"
                )
            grouped[(accession, unimod)].append(row)

    output: dict[tuple[str, str], SemanticEvidence] = {}
    for key, rows in grouped.items():
        statuses = {(row.get("evidence_status") or "").strip().casefold() for row in rows}
        if "supported" in statuses and "conflicting" in statuses:
            status = "mixed"
        elif "conflicting" in statuses:
            status = "conflicting"
        elif "supported" in statuses:
            status = "supported"
        else:
            status = "not-found"
        sources = tuple(
            sorted(
                {
                    (row.get("evidence_source") or "").strip()
                    for row in rows
                    if (row.get("evidence_source") or "").strip()
                }
            )
        )
        notes = tuple(
            sorted(
                {
                    (row.get("evidence_note") or "").strip()
                    for row in rows
                    if (row.get("evidence_note") or "").strip()
                }
            )
        )
        output[key] = SemanticEvidence(status=status, sources=sources, notes=notes)
    return output


def _annotation(result: AnalysisResult, field: str) -> Annotation | None:
    return next((item for item in result.annotations if item.field == field), None)


def _payload(result: AnalysisResult, field: str) -> dict[str, Any] | None:
    item = _annotation(result, field)
    return item.value if item is not None and isinstance(item.value, dict) else None


def _metric(result: AnalysisResult, key: str) -> float | None:
    for item in result.metrics:
        if item.key != key or not isinstance(item.value, (int, float)):
            continue
        value = float(item.value)
        if math.isfinite(value) and value > 0:
            return value
    return None


def _positive_number(value: Any) -> float | None:
    if not isinstance(value, (int, float)):
        return None
    result = float(value)
    return result if math.isfinite(result) and result > 0 else None


def _result_features(result: AnalysisResult) -> dict[str, Any]:
    precursor = _payload(result, "suggested_precursor_search_tolerance_ppm")
    fragment_ppm = _payload(result, "suggested_fragment_search_tolerance_ppm")
    fragment_da = _payload(result, "suggested_fragment_search_tolerance_da")
    fragment = fragment_ppm or fragment_da
    acquisition = _annotation(result, "acquisition_method")
    instruments = sorted(
        {
            term.accession or term.name
            for term in result.metadata.instruments
            if term.accession or term.name
        }
    )
    return {
        "precursor": _positive_number(
            precursor.get("suggested_tolerance") if precursor else None
        ),
        "fragment": _positive_number(fragment.get("suggested_tolerance") if fragment else None),
        "fragment_unit": "ppm" if fragment_ppm else "Da" if fragment_da else "",
        "fragment_regime": str(fragment.get("resolution_regime") or "") if fragment else "",
        "rt": _metric(result, "ChromatographyDuration"),
        "isolation": _metric(result, "IsolationWidth_MS2_Median"),
        "ms1": _metric(result, "NumberOfSpectra_MS1"),
        "ms2": _metric(result, "NumberOfSpectra_MS2"),
        "instrument": ";".join(instruments),
        "acquisition": str(acquisition.value or "") if acquisition else "",
    }


def _robust_scaled_matrix(
    names: list[str], features: dict[str, dict[str, Any]]
) -> tuple[np.ndarray, list[str]] | None:
    minimum_observed = max(3, ceil(len(names) * 0.80))
    columns: list[list[float]] = []
    selected: list[str] = []
    for key in _NUMERIC_FEATURES:
        observed = [
            float(features[name][key])
            for name in names
            if _positive_number(features[name].get(key)) is not None
        ]
        if len(observed) < minimum_observed:
            continue
        if len({round(value, 12) for value in observed}) < 2:
            continue
        fill = median(observed)
        transformed = [
            log10(float(features[name][key]))
            if _positive_number(features[name].get(key)) is not None
            else log10(fill)
            for name in names
        ]
        med = median(transformed)
        q25, q75 = np.quantile(np.asarray(transformed, dtype=float), [0.25, 0.75])
        scale = float(q75 - q25)
        if scale <= 1e-12:
            mad = median(abs(value - med) for value in transformed)
            scale = 1.4826 * mad
        if scale <= 1e-12:
            continue
        columns.append([(value - med) / scale for value in transformed])
        selected.append(key)
    if not columns:
        return None
    return np.asarray(columns, dtype=float).T, selected


def _two_means(matrix: np.ndarray) -> np.ndarray | None:
    if matrix.shape[0] < 2:
        return None
    distances = np.sum((matrix[:, None, :] - matrix[None, :, :]) ** 2, axis=2)
    first, second = np.unravel_index(int(np.argmax(distances)), distances.shape)
    if first == second or distances[first, second] <= 1e-12:
        return None
    centers = np.vstack([matrix[first], matrix[second]]).astype(float)
    labels = np.zeros(matrix.shape[0], dtype=int)
    previous: np.ndarray | None = None
    for _ in range(100):
        squared = np.sum((matrix[:, None, :] - centers[None, :, :]) ** 2, axis=2)
        labels = np.argmin(squared, axis=1)
        if previous is not None and np.array_equal(labels, previous):
            break
        if len(set(labels.tolist())) < 2:
            return None
        previous = labels.copy()
        centers = np.vstack([matrix[labels == index].mean(axis=0) for index in range(2)])
    return labels


def _silhouette(matrix: np.ndarray, labels: np.ndarray) -> float:
    if matrix.shape[0] < 3 or len(set(labels.tolist())) < 2:
        return -1.0
    distances = np.sqrt(np.sum((matrix[:, None, :] - matrix[None, :, :]) ** 2, axis=2))
    values: list[float] = []
    for index, label in enumerate(labels):
        same = np.where(labels == label)[0]
        same = same[same != index]
        if same.size == 0:
            values.append(0.0)
            continue
        a = float(distances[index, same].mean())
        other_means = [
            float(distances[index, np.where(labels == other)[0]].mean())
            for other in sorted(set(labels.tolist()))
            if other != label
        ]
        b = min(other_means)
        denominator = max(a, b)
        values.append((b - a) / denominator if denominator > 0 else 0.0)
    return float(sum(values) / len(values))


def _numeric_values(
    group: list[str], features: dict[str, dict[str, Any]], key: str
) -> list[float]:
    return [
        float(features[name][key])
        for name in group
        if _positive_number(features[name].get(key)) is not None
    ]


def _supported_numeric_split(
    group: list[str], features: dict[str, dict[str, Any]], min_group: int
) -> tuple[list[list[str]], str] | None:
    matrix_and_keys = _robust_scaled_matrix(group, features)
    if matrix_and_keys is None:
        return None
    matrix, feature_keys = matrix_and_keys
    labels = _two_means(matrix)
    if labels is None:
        return None
    buckets = [
        [name for name, label in zip(group, labels, strict=True) if int(label) == cluster]
        for cluster in range(2)
    ]
    if any(len(bucket) < min_group for bucket in buckets):
        return None
    score = _silhouette(matrix, labels)
    if score < 0.45:
        return None
    differences: list[tuple[float, str]] = []
    for key in feature_keys:
        medians = [median(_numeric_values(bucket, features, key)) for bucket in buckets]
        if min(medians) <= 0:
            continue
        ratio = max(medians) / min(medians)
        if ratio >= _SEPARATION_THRESHOLDS[key]:
            differences.append((ratio, key))
    if not differences:
        return None
    buckets.sort(key=lambda bucket: min(group.index(name) for name in bucket))
    differences.sort(reverse=True)
    evidence = ", ".join(
        f"{_FEATURE_NAMES[key]} {ratio:.2f}x" for ratio, key in differences[:3]
    )
    return buckets, f"multivariate k=2, silhouette={score:.2f}; {evidence}"


def _supported_categorical_split(
    group: list[str], features: dict[str, dict[str, Any]], key: str, min_group: int
) -> list[list[str]] | None:
    values = {name: str(features[name].get(key) or "") for name in group}
    if any(not value for value in values.values()):
        return None
    buckets: dict[str, list[str]] = defaultdict(list)
    for name, value in values.items():
        buckets[value].append(name)
    if len(buckets) <= 1 or any(len(bucket) < min_group for bucket in buckets.values()):
        return None
    return [buckets[value] for value in sorted(buckets)]


def _supported_duration_split(
    group: list[str], features: dict[str, dict[str, Any]], min_group: int
) -> tuple[list[list[str]], str] | None:
    """Split clearly separated chromatography-duration regimes conservatively."""
    values: list[tuple[float, str]] = []
    for name in group:
        value = _positive_number(features[name].get("rt"))
        if value is None:
            return None
        values.append((value, name))
    if len(values) < 2 * min_group:
        return None
    values.sort()
    best: tuple[float, list[list[str]], float, float] | None = None
    for cut in range(min_group, len(values) - min_group + 1):
        left = values[:cut]
        right = values[cut:]
        left_values = [item[0] for item in left]
        right_values = [item[0] for item in right]
        median_ratio = median(right_values) / median(left_values)
        boundary_ratio = right_values[0] / left_values[-1]
        if median_ratio < _SEPARATION_THRESHOLDS["rt"] or boundary_ratio < 1.10:
            continue
        matrix = np.asarray([[log10(value)] for value, _ in values], dtype=float)
        labels = np.asarray([0] * len(left) + [1] * len(right), dtype=int)
        score = _silhouette(matrix, labels)
        if score < 0.60:
            continue
        buckets = [[name for _, name in left], [name for _, name in right]]
        if best is None or score > best[0]:
            best = (score, buckets, median_ratio, boundary_ratio)
    if best is None:
        return None
    score, buckets, median_ratio, boundary_ratio = best
    return (
        buckets,
        "univariate chromatography duration split, "
        f"silhouette={score:.2f}; median ratio={median_ratio:.2f}x; "
        f"boundary gap={boundary_ratio:.2f}x",
    )


def _experiment_groups(
    results: list[AnalysisResult], max_groups: int = 4
) -> tuple[dict[str, str], list[str], dict[str, dict[str, Any]]]:
    names = [result.input_path.name for result in results]
    features = {result.input_path.name: _result_features(result) for result in results}
    if not names:
        return {}, [], features
    min_group = max(3, ceil(len(names) * 0.04))
    groups = [names]
    evidence: list[str] = []
    for key in ("fragment_unit", "fragment_regime", "instrument", "acquisition"):
        updated: list[list[str]] = []
        for group in groups:
            split = _supported_categorical_split(group, features, key, min_group)
            if split is not None and len(updated) + len(split) <= max_groups:
                updated.extend(split)
                evidence.append(f"categorical {key} split")
            else:
                updated.append(group)
        groups = updated
    changed = True
    while changed and len(groups) < max_groups:
        changed = False
        updated = []
        for index, group in enumerate(groups):
            remaining = len(groups) - index - 1
            if len(updated) + remaining + 2 > max_groups:
                updated.append(group)
                continue
            numeric_split = _supported_duration_split(group, features, min_group)
            if numeric_split is None:
                numeric_split = _supported_numeric_split(group, features, min_group)
            if numeric_split is None:
                updated.append(group)
                continue
            buckets, reason = numeric_split
            updated.extend(buckets)
            evidence.append(reason)
            changed = True
        groups = updated
    groups.sort(key=lambda group: min(names.index(name) for name in group))
    assignments: dict[str, str] = {}
    for index, group in enumerate(groups, start=1):
        label = f"Experiment group {index}"
        assignments.update({name: label for name in group})
    return assignments, evidence, features


def _format_tolerance(value: float, unit: str) -> str:
    if unit == "ppm":
        return f"{math.ceil(value):g} ppm"
    # Preserve useful low-resolution precision while rounding upward so the
    # shared setting still covers the maximum supported per-run estimate.
    rounded = math.ceil(value * 1000.0 - 1e-12) / 1000.0
    return f"{rounded:g} Da"


def _group_tolerance_annotation(
    results: list[AnalysisResult], label: str, field: str
) -> Annotation | None:
    if field == COHORT_PRECURSOR_FIELD:
        source_fields: tuple[str, ...] = ("suggested_precursor_search_tolerance_ppm",)
        sdrf_column = "comment[precursor mass tolerance]"
    else:
        source_fields = (
            "suggested_fragment_search_tolerance_ppm",
            "suggested_fragment_search_tolerance_da",
        )
        sdrf_column = "comment[fragment mass tolerance]"
    estimates: list[tuple[float, str]] = []
    source_measurement_support = 0
    for result in results:
        for source_field in source_fields:
            payload = _payload(result, source_field)
            if not payload:
                continue
            value = _positive_number(payload.get("suggested_tolerance"))
            unit = str(payload.get("unit") or "")
            if value is not None and unit in {"ppm", "Da"}:
                estimates.append((value, unit))
                annotation = _annotation(result, source_field)
                source_measurement_support += int(annotation.support or 0) if annotation else 0
                break
    required = max(1, ceil(len(results) * 0.80))
    if len(estimates) < required:
        return None
    units = {unit for _, unit in estimates}
    if len(units) != 1:
        return None
    unit = next(iter(units))
    values = [value for value, _ in estimates]
    common = max(values)
    formatted_value = _format_tolerance(common, unit)
    return Annotation(
        field,
        {
            "experiment_group": label,
            "unit": unit,
            "median": median(values),
            "common_max": common,
            "files_with_estimate": len(estimates),
            "group_files": len(results),
            "source_measurement_support": source_measurement_support,
        },
        EvidenceKind.INFERRED,
        "prideQC cohort reanalysis tolerance synthesis v1",
        (
            "Per-run repeat-spectrum precision estimates are combined within a conservatively "
            "inferred experiment group. The SDRF proposal uses the maximum supported per-run "
            "tolerance so one shared reanalysis setting covers all supported runs in the group; "
            "at least 80% of group runs must provide a compatible estimate."
        ),
        support=len(estimates),
        total=len(results),
        sdrf_column=sdrf_column,
        sdrf_value=formatted_value,
    )


def _mass_shift_records(result: AnalysisResult) -> list[dict[str, Any]]:
    item = _annotation(result, "putative_modification_mass_shifts")
    if item is None or not isinstance(item.value, list):
        return []
    return [record for record in item.value if isinstance(record, dict)]


def _beta_prevalence_probability(hits: int, runs: int, threshold: float = 0.10) -> float:
    if runs <= 0 or hits < 0 or hits > runs:
        return 0.0
    total = runs + 1
    probability = sum(
        comb(total, index) * threshold**index * (1.0 - threshold) ** (total - index)
        for index in range(hits + 1)
    )
    return min(1.0, max(0.0, probability))


def _cluster_mass_shift_records(
    observations: list[tuple[float, str, dict[str, Any]]], tolerance_da: float = 0.02
) -> list[list[tuple[float, str, dict[str, Any]]]]:
    observations.sort(key=lambda item: item[0])
    families: list[list[tuple[float, str, dict[str, Any]]]] = []
    for observation in observations:
        if not families:
            families.append([observation])
            continue
        center = median(item[0] for item in families[-1])
        if abs(observation[0] - center) <= tolerance_da:
            families[-1].append(observation)
        else:
            families.append([observation])
    return families


def _sdrf_candidates(record: dict[str, Any]) -> dict[str, tuple[str, str]]:
    """Return non-decoy/non-substitution UniMod candidates relevant to SDRF identity."""
    raw = record.get("diagnostic_unimod_candidates") or record.get("unimod_candidates") or []
    result: dict[str, tuple[str, str]] = {}
    if not isinstance(raw, list):
        return result
    for candidate in raw:
        if not isinstance(candidate, dict):
            continue
        category = str(candidate.get("candidate_category") or "")
        if category in {"decoy", "amino-acid-substitution"}:
            continue
        accession = str(candidate.get("unimod_accession") or "").strip()
        name = str(candidate.get("name") or "").strip()
        if accession and name:
            result[accession] = (name, category)
    return result


def _group_ptm_families(
    results: list[AnalysisResult],
    label: str,
    *,
    semantic_evidence: dict[tuple[str, str], SemanticEvidence] | None = None,
    project_accession: str | None = None,
) -> list[dict[str, Any]]:
    """Return strict RAW-supported PTM families with an independent semantic gate."""
    observations: list[tuple[float, str, dict[str, Any]]] = []
    for result in results:
        for record in _mass_shift_records(result):
            value = record.get("delta_mass_da")
            if isinstance(value, (int, float)) and math.isfinite(float(value)):
                observations.append((float(value), result.input_path.name, record))
    families = _cluster_mass_shift_records(observations)
    rows: list[dict[str, Any]] = []
    if len(results) < 3:
        return rows
    accession_key = (project_accession or "").strip().upper()
    evidence_lookup = semantic_evidence or {}
    for family in families:
        run_hits = {item[1] for item in family}
        if len(run_hits) < 3:
            continue
        prevalence = len(run_hits) / len(results) if results else 0.0
        probability = _beta_prevalence_probability(len(run_hits), len(results))
        high_support_runs = {
            run_name
            for _, run_name, record in family
            if str(record.get("confidence") or "") == "high-support"
        }
        high_support_fraction = len(high_support_runs) / len(run_hits) if run_hits else 0.0
        if prevalence < 0.90 or high_support_fraction < 0.80 or probability < 0.99:
            continue
        candidate_union: dict[str, tuple[str, str]] = {}
        candidate_runs: set[str] = set()
        for _, run_name, record in family:
            candidates = _sdrf_candidates(record)
            candidate_union.update(candidates)
            if candidates:
                candidate_runs.add(run_name)
        if len(candidate_union) != 1:
            continue
        accession, (name, category) = next(iter(candidate_union.items()))
        candidate_run_fraction = len(candidate_runs) / len(run_hits) if run_hits else 0.0
        if category != "biological-ptm" or candidate_run_fraction < 0.80:
            continue
        semantic = evidence_lookup.get((accession_key, accession.upper()))
        semantic_status = semantic.status if semantic is not None else "not-evaluated"
        sdrf_status = (
            "eligible" if semantic_status == "supported" else "hold-semantic-not-supported"
        )
        rows.append(
            {
                "experiment_group": label,
                "unimod_accession": accession,
                "unimod_name": name,
                "median_mass_da": median(item[0] for item in family),
                "family_runs": len(run_hits),
                "group_runs": len(results),
                "run_prevalence": prevalence,
                "high_support_run_fraction": high_support_fraction,
                "candidate_run_fraction": candidate_run_fraction,
                "raw_prevalence_probability": probability,
                "median_pair_support": median(
                    int(item[2].get("pair_support", 0) or 0) for item in family
                ),
                "mass_identity_ambiguous": False,
                "semantic_evidence_status": semantic_status,
                "semantic_evidence_sources": list(semantic.sources) if semantic else [],
                "semantic_evidence_notes": list(semantic.notes) if semantic else [],
                "sdrf_status": sdrf_status,
            }
        )
    rows.sort(
        key=lambda row: (
            float(row["run_prevalence"]),
            float(row["raw_prevalence_probability"]),
            int(row["family_runs"]),
        ),
        reverse=True,
    )
    return rows[:5]


def synthesize_cohort(
    results: list[AnalysisResult],
    *,
    semantic_evidence: dict[tuple[str, str], SemanticEvidence] | None = None,
    project_accession: str | None = None,
) -> CohortSynthesis:
    """Attach conservative cohort-level SDRF proposals to successful results."""
    if not results:
        return CohortSynthesis({}, {}, [], [], [])
    assignments, evidence, features = _experiment_groups(results)
    by_name = {result.input_path.name: result for result in results}
    grouped: dict[str, list[AnalysisResult]] = defaultdict(list)
    for name, label in assignments.items():
        grouped[label].append(by_name[name])

    group_summaries: dict[str, dict[str, Any]] = {}
    eligible_ptm_families: list[dict[str, Any]] = []
    review_ptm_families: list[dict[str, Any]] = []
    for label, members in grouped.items():
        durations = _numeric_values([item.input_path.name for item in members], features, "rt")
        precursor = _group_tolerance_annotation(members, label, COHORT_PRECURSOR_FIELD)
        fragment = _group_tolerance_annotation(members, label, COHORT_FRAGMENT_FIELD)
        ptm_review = _group_ptm_families(
            members,
            label,
            semantic_evidence=semantic_evidence,
            project_accession=project_accession,
        )
        ptms = [item for item in ptm_review if item["sdrf_status"] == "eligible"]
        review_ptm_families.extend(ptm_review)
        eligible_ptm_families.extend(ptms)
        summary: dict[str, Any] = {
            "files": len(members),
            "members": [item.input_path.name for item in members],
        }
        if durations:
            summary.update(
                {
                    "run_duration_median_min": median(durations) / 60.0,
                    "run_duration_min_min": min(durations) / 60.0,
                    "run_duration_max_min": max(durations) / 60.0,
                }
            )
        if precursor is not None:
            summary["precursor_tolerance"] = precursor.value
        if fragment is not None:
            summary["fragment_tolerance"] = fragment.value
        summary["ptm_review_families"] = ptm_review
        summary["sdrf_eligible_ptms"] = ptms
        group_summaries[label] = summary

        for result in members:
            result.annotations.append(
                Annotation(
                    EXPERIMENT_GROUP_FIELD,
                    label,
                    EvidenceKind.INFERRED,
                    "prideQC conservative multi-feature experiment grouping v1",
                    (
                        "Grouping uses compatible fragment units/regimes, instrument/acquisition "
                        "metadata when complete, and robustly scaled precursor tolerance, fragment "
                        "tolerance, chromatography duration, MS2 isolation width and MS1/MS2 "
                        "counts. Strong, well-separated chromatography-duration regimes are "
                        "split before multivariate clustering. It indicates acquisition "
                        "heterogeneity and is not proof of separate biological experiments."
                    ),
                    support=len(members),
                    total=len(results),
                )
            )
            if precursor is not None:
                result.annotations.append(precursor)
            if fragment is not None:
                result.annotations.append(fragment)
            for ptm in ptms:
                result.annotations.append(
                    Annotation(
                        COHORT_MODIFICATION_FIELD,
                        ptm,
                        EvidenceKind.INFERRED,
                        "prideQC recurrent PTM mass-family synthesis v2",
                        (
                            "Automatic SDRF PTM write-back requires both strict recurrent RAW "
                            "mass-family support and independent semantic/study evidence that "
                            "supports the exact UniMod accession. RAW recurrence alone remains "
                            "review-only and does not establish peptide or site localization."
                        ),
                        support=int(ptm["family_runs"]),
                        total=int(ptm["group_runs"]),
                        sdrf_column="comment[modification parameters]",
                        sdrf_value=CVTerm(
                            str(ptm["unimod_accession"]), str(ptm["unimod_name"])
                        ).sdrf_value(),
                    )
                )

    return CohortSynthesis(
        assignments,
        group_summaries,
        evidence,
        eligible_ptm_families,
        review_ptm_families,
    )
