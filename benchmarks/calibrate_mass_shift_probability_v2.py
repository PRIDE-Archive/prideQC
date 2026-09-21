#!/usr/bin/env python3
"""Cluster-level PTM mass-annotation calibration for frozen prideQC mass shifts.

v2 corrects two failure modes found in the first target/decoy prototype:

* the statistical unit is one best biological-PTM mass annotation per recurrent
  cluster, not every mass-compatible UniMod row on that cluster;
* the null is sampled by repeatedly shifting each observed cluster mass against
  the *fixed* biological-PTM catalogue, preserving cluster evidence while
  estimating chance catalogue matches with a dense deterministic null.

``mass_annotation_probability`` is the empirical probability that the best
biological-PTM mass annotation for a recurrent cluster is enriched beyond this
shifted-mass null, conditional on the cluster evidence. It is not a peptide,
site-localization, or PTM-identity posterior probability.
"""

from __future__ import annotations

import argparse
import bisect
import csv
import importlib.util
import json
import math
import statistics
import sys
from collections import Counter, defaultdict
from collections.abc import Sequence
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

_V1_PATH = Path(__file__).with_name("calibrate_mass_shift_probability_v1.py")
_V1_SPEC = importlib.util.spec_from_file_location("prideqc_probability_v1", _V1_PATH)
if _V1_SPEC is None or _V1_SPEC.loader is None:  # pragma: no cover - import guard
    raise RuntimeError(f"Cannot load v1 probability helpers from {_V1_PATH}")
v1 = importlib.util.module_from_spec(_V1_SPEC)
sys.modules[_V1_SPEC.name] = v1
_V1_SPEC.loader.exec_module(v1)

MODEL_VERSION = "mass-shift-probability-v2.0"
DEFAULT_DECOY_REPLICATES = 128
DEFAULT_DECOY_MIN_OFFSET_DA = 0.05
DEFAULT_DECOY_MAX_OFFSET_DA = 5.0
DEFAULT_MINIMUM_MASS_DA = 0.5
DEFAULT_MAXIMUM_MASS_DA = 500.0
_GOLDEN_RATIO_CONJUGATE = 0.6180339887498949


@dataclass(frozen=True, slots=True)
class ClusterAnnotation:
    cluster: Any
    unimod_accession: str
    unimod_name: str
    theoretical_mass_da: float
    residual_da: float
    candidate_count: int
    candidate_dominance: float
    runner_up_residual_da: float | None
    residual_margin_da: float | None
    score: float


@dataclass(frozen=True, slots=True)
class CalibrationBin:
    minimum_score: float
    maximum_score: float
    target_count: float
    weighted_decoy_count: float
    local_fdr: float
    probability: float


@dataclass(frozen=True, slots=True)
class ProbabilityModel:
    model_version: str
    probability_semantics: str
    score_formula: str
    decoy_replicates: int
    decoy_min_offset_da: float
    decoy_max_offset_da: float
    decoy_weight: float
    catalog_sha256: str
    biological_catalog_records: int
    training_accessions: tuple[str, ...]
    training_target_clusters: int
    training_decoy_hits: int
    calibration_bins: tuple[CalibrationBin, ...]
    q_thresholds: tuple[tuple[float, float], ...]


def _candidate_weight(residual_da: float, tolerance_da: float) -> float:
    """Mass-only likelihood weight used solely to quantify candidate dominance."""

    sigma = max(abs(tolerance_da) / 4.0, 1e-6)
    z = abs(residual_da) / sigma
    return math.exp(-0.5 * z * z)


def annotation_score(
    cluster: Any,
    *,
    residual_da: float,
    candidate_dominance: float,
) -> float:
    """Evidence axis shared by target and shifted-null annotations."""

    base = v1.evidence_score(
        pair_support=cluster.pair_support,
        unique_spectrum_support=cluster.unique_spectrum_support,
        median_spectral_similarity=cluster.median_spectral_similarity,
        support_fraction=cluster.support_fraction,
        cluster_sigma_da=cluster.cluster_sigma_da,
        match_tolerance_da=cluster.match_tolerance_da,
        residual_da=residual_da,
    )
    # Dominance is a mass-ambiguity penalty, not an identity prior. A unique
    # candidate adds no penalty; an exact two-way isobaric tie costs ~1.04.
    return base + 1.5 * math.log(max(candidate_dominance, 1e-12))


def _best_candidate(
    *,
    cluster: Any,
    candidates: Sequence[tuple[str, str, float, float]],
) -> ClusterAnnotation:
    """Return the closest candidate plus a normalized mass-dominance measure.

    Candidate tuples are ``(accession, name, theoretical_mass, residual)``.
    """

    if not candidates:
        raise ValueError("At least one candidate is required.")
    ordered = sorted(
        candidates,
        key=lambda item: (abs(item[3]), item[0].casefold(), item[1].casefold()),
    )
    weights = [
        _candidate_weight(item[3], cluster.match_tolerance_da) for item in ordered
    ]
    total_weight = sum(weights)
    dominance = weights[0] / total_weight if total_weight > 0.0 else 1.0 / len(weights)
    best = ordered[0]
    runner = abs(ordered[1][3]) if len(ordered) > 1 else None
    best_residual = abs(best[3])
    margin = runner - best_residual if runner is not None else None
    return ClusterAnnotation(
        cluster=cluster,
        unimod_accession=best[0],
        unimod_name=best[1],
        theoretical_mass_da=best[2],
        residual_da=best[3],
        candidate_count=len(ordered),
        candidate_dominance=dominance,
        runner_up_residual_da=runner,
        residual_margin_da=margin,
        score=annotation_score(
            cluster,
            residual_da=best[3],
            candidate_dominance=dominance,
        ),
    )


def biological_cluster_annotations(
    candidates: Sequence[Any],
) -> list[ClusterAnnotation]:
    """Collapse displayed candidate rows to one best biological annotation/cluster."""

    grouped: dict[str, list[Any]] = defaultdict(list)
    for item in candidates:
        if (
            item.candidate_category == "biological-ptm"
            and item.cluster.classification == "putative-ptm"
        ):
            grouped[item.cluster.cluster_key].append(item)

    result: list[ClusterAnnotation] = []
    for rows in grouped.values():
        # Deduplicate repeated specificity rows if an upstream format ever emits them.
        unique: dict[tuple[str, str, float], Any] = {}
        for item in rows:
            key = (
                item.unimod_accession,
                item.unimod_name,
                round(float(item.theoretical_mass_da), 9),
            )
            current = unique.get(key)
            if current is None or abs(item.residual_da) < abs(current.residual_da):
                unique[key] = item
        values = list(unique.values())
        cluster = values[0].cluster
        result.append(
            _best_candidate(
                cluster=cluster,
                candidates=[
                    (
                        item.unimod_accession,
                        item.unimod_name,
                        float(item.theoretical_mass_da),
                        float(item.residual_da),
                    )
                    for item in values
                ],
            )
        )
    return sorted(
        result,
        key=lambda item: (
            item.cluster.accession,
            item.cluster.data_file,
            item.cluster.cluster_rank,
        ),
    )


def _catalog_index(catalog: Sequence[Any]) -> tuple[list[float], list[Any]]:
    biological = sorted(
        (item for item in catalog if item.category == "biological-ptm"),
        key=lambda item: (item.mass_da, item.accession, item.name),
    )
    return [float(item.mass_da) for item in biological], biological


def deterministic_decoy_offsets(
    replicates: int,
    *,
    minimum_offset_da: float,
    maximum_offset_da: float,
) -> tuple[float, ...]:
    if replicates < 1:
        raise ValueError("decoy replicates must be >= 1")
    if not 0.0 < minimum_offset_da < maximum_offset_da:
        raise ValueError("decoy offset bounds must satisfy 0 < minimum < maximum")
    offsets = []
    width = maximum_offset_da - minimum_offset_da
    for index in range(1, replicates + 1):
        fraction = (index * _GOLDEN_RATIO_CONJUGATE) % 1.0
        magnitude = minimum_offset_da + fraction * width
        sign = -1.0 if index % 2 else 1.0
        offsets.append(sign * magnitude)
    return tuple(offsets)


def _wrap_mass(mass: float, minimum: float, maximum: float) -> float:
    width = maximum - minimum
    if width <= 0.0:
        raise ValueError("mass range must have positive width")
    return minimum + ((mass - minimum) % width)


def shifted_decoy_scores(
    targets: Sequence[ClusterAnnotation],
    catalog: Sequence[Any],
    *,
    offsets: Sequence[float],
    minimum_mass_da: float,
    maximum_mass_da: float,
) -> list[tuple[str, float]]:
    """Sample a dense cluster-preserving null against the fixed PTM catalogue."""

    masses, records = _catalog_index(catalog)
    result: list[tuple[str, float]] = []
    for target in targets:
        cluster = target.cluster
        tolerance = abs(cluster.match_tolerance_da)
        for offset in offsets:
            shifted = _wrap_mass(
                cluster.delta_mass_da + offset,
                minimum_mass_da,
                maximum_mass_da,
            )
            lo = bisect.bisect_left(masses, shifted - tolerance)
            hi = bisect.bisect_right(masses, shifted + tolerance)
            if lo == hi:
                continue
            candidates = []
            for record in records[lo:hi]:
                residual = shifted - float(record.mass_da)
                candidates.append(
                    (
                        str(record.accession),
                        str(record.name),
                        float(record.mass_da),
                        residual,
                    )
                )
            decoy = _best_candidate(cluster=cluster, candidates=candidates)
            result.append((cluster.accession, decoy.score))
    return result


def _pava_non_decreasing(values: Sequence[float], weights: Sequence[float]) -> list[float]:
    if len(values) != len(weights):
        raise ValueError("PAVA values and weights must have equal length.")
    blocks: list[dict[str, Any]] = []
    for index, (value, weight) in enumerate(zip(values, weights, strict=True)):
        blocks.append({"start": index, "end": index, "weight": weight, "mean": value})
        while len(blocks) >= 2 and blocks[-2]["mean"] > blocks[-1]["mean"]:
            right = blocks.pop()
            left = blocks.pop()
            total_weight = left["weight"] + right["weight"]
            mean = (
                left["mean"] * left["weight"] + right["mean"] * right["weight"]
            ) / total_weight
            blocks.append(
                {
                    "start": left["start"],
                    "end": right["end"],
                    "weight": total_weight,
                    "mean": mean,
                }
            )
    fitted = [0.0] * len(values)
    for block in blocks:
        for index in range(block["start"], block["end"] + 1):
            fitted[index] = float(block["mean"])
    return fitted


def fit_probability_model(
    target_scores: Sequence[float],
    decoy_scores: Sequence[float],
    *,
    decoy_replicates: int,
    decoy_min_offset_da: float,
    decoy_max_offset_da: float,
    catalog_hash: str,
    biological_catalog_records: int,
    training_accessions: Sequence[str],
    minimum_targets_per_bin: int = 25,
) -> ProbabilityModel:
    if decoy_replicates < 1:
        raise ValueError("decoy_replicates must be >= 1")
    if not target_scores:
        raise ValueError("Cannot fit probability model without target clusters.")
    targets = sorted(float(x) for x in target_scores)
    decoys = sorted(float(x) for x in decoy_scores)
    decoy_weight = 1.0 / decoy_replicates

    target_score_counts = Counter(targets)
    score_counts = [
        (score, target_score_counts[score])
        for score in sorted(target_score_counts, reverse=True)
    ]
    bin_ranges: list[tuple[float, float]] = []
    bin_high: float | None = None
    bin_count = 0
    for score, count in score_counts:
        if bin_high is None:
            bin_high = score
        bin_count += count
        if bin_count >= max(1, minimum_targets_per_bin):
            bin_ranges.append((score, bin_high))
            bin_high = None
            bin_count = 0
    if bin_high is not None:
        if bin_ranges:
            _old_minimum, previous_maximum = bin_ranges[-1]
            bin_ranges[-1] = (score_counts[-1][0], previous_maximum)
        else:
            bin_ranges.append((score_counts[-1][0], bin_high))

    raw_fdr: list[float] = []
    weights: list[float] = []
    counts: list[tuple[float, float]] = []
    for minimum, maximum in bin_ranges:
        target_count = sum(minimum <= score <= maximum for score in targets)
        raw_decoys = sum(minimum <= score <= maximum for score in decoys)
        weighted_decoys = raw_decoys * decoy_weight
        # Jeffreys-like half-event floor is on the target-equivalent scale, not
        # divided by the number of null replicates. This prevents zero-decoy
        # bins from becoming effectively certain simply because the null was
        # oversampled many times.
        local_fdr = min(1.0, (weighted_decoys + 0.5) / max(target_count, 1))
        raw_fdr.append(local_fdr)
        weights.append(max(float(target_count), 1.0))
        counts.append((float(target_count), weighted_decoys))

    fitted_fdr = _pava_non_decreasing(raw_fdr, weights)
    bins = tuple(
        CalibrationBin(
            minimum_score=minimum,
            maximum_score=maximum,
            target_count=counts[index][0],
            weighted_decoy_count=counts[index][1],
            local_fdr=min(max(fitted_fdr[index], 0.0), 1.0),
            probability=1.0 - min(max(fitted_fdr[index], 0.0), 1.0),
        )
        for index, (minimum, maximum) in enumerate(bin_ranges)
    )

    # Higher score is better. A target observed at score s is included by any
    # threshold <= s, so its q-value is the minimum tail FDR over that threshold
    # and all less-stringent (lower-score) thresholds. Thresholds are stored
    # high -> low, therefore the cumulative minimum is taken from the end.
    thresholds = sorted(set(targets + decoys), reverse=True)
    raw_q: list[tuple[float, float]] = []
    for threshold in thresholds:
        target_tail = len(targets) - bisect.bisect_left(targets, threshold)
        decoy_tail = len(decoys) - bisect.bisect_left(decoys, threshold)
        weighted_decoys = decoy_tail * decoy_weight
        fdr = min(1.0, (weighted_decoys + 0.5) / max(target_tail, 1))
        raw_q.append((threshold, fdr))
    running = 1.0
    q_values = [1.0] * len(raw_q)
    for index in range(len(raw_q) - 1, -1, -1):
        running = min(running, raw_q[index][1])
        q_values[index] = running
    q_thresholds = [
        (raw_q[index][0], q_values[index]) for index in range(len(raw_q))
    ]

    return ProbabilityModel(
        model_version=MODEL_VERSION,
        probability_semantics=(
            "Empirical probability that the best biological-PTM mass annotation for a "
            "recurrent cluster is enriched beyond a dense cluster-preserving shifted-mass "
            "null, conditional on prideQC cluster evidence; not peptide/site/PTM-identity "
            "probability."
        ),
        score_formula=(
            "v1 cluster evidence score + 1.5*log(candidate_mass_dominance), with exactly "
            "one closest biological-PTM candidate per cluster"
        ),
        decoy_replicates=decoy_replicates,
        decoy_min_offset_da=decoy_min_offset_da,
        decoy_max_offset_da=decoy_max_offset_da,
        decoy_weight=decoy_weight,
        catalog_sha256=catalog_hash,
        biological_catalog_records=biological_catalog_records,
        training_accessions=tuple(sorted(training_accessions)),
        training_target_clusters=len(targets),
        training_decoy_hits=len(decoys),
        calibration_bins=bins,
        q_thresholds=tuple(q_thresholds),
    )


def probability_for_score(model: ProbabilityModel, score: float) -> float:
    if not model.calibration_bins:
        return 0.0
    for item in model.calibration_bins:
        if item.minimum_score <= score <= item.maximum_score:
            return item.probability
    if score > model.calibration_bins[0].maximum_score:
        return model.calibration_bins[0].probability
    return model.calibration_bins[-1].probability


def q_value_for_score(model: ProbabilityModel, score: float) -> float:
    if not model.q_thresholds:
        return 1.0
    for threshold, q_value in model.q_thresholds:
        if score >= threshold:
            return q_value
    return model.q_thresholds[-1][1]


def model_to_json(model: ProbabilityModel) -> dict[str, Any]:
    payload = asdict(model)
    payload["calibration_bins"] = [asdict(item) for item in model.calibration_bins]
    payload["q_thresholds"] = [list(item) for item in model.q_thresholds]
    return payload


def model_from_json(payload: dict[str, Any]) -> ProbabilityModel:
    return ProbabilityModel(
        model_version=str(payload["model_version"]),
        probability_semantics=str(payload["probability_semantics"]),
        score_formula=str(payload["score_formula"]),
        decoy_replicates=int(payload["decoy_replicates"]),
        decoy_min_offset_da=float(payload["decoy_min_offset_da"]),
        decoy_max_offset_da=float(payload["decoy_max_offset_da"]),
        decoy_weight=float(payload["decoy_weight"]),
        catalog_sha256=str(payload["catalog_sha256"]),
        biological_catalog_records=int(payload["biological_catalog_records"]),
        training_accessions=tuple(str(x) for x in payload["training_accessions"]),
        training_target_clusters=int(payload["training_target_clusters"]),
        training_decoy_hits=int(payload["training_decoy_hits"]),
        calibration_bins=tuple(CalibrationBin(**item) for item in payload["calibration_bins"]),
        q_thresholds=tuple((float(x[0]), float(x[1])) for x in payload["q_thresholds"]),
    )


def scored_annotation_row(
    item: ClusterAnnotation,
    model: ProbabilityModel,
    *,
    calibration_scope: str,
) -> dict[str, Any]:
    probability = probability_for_score(model, item.score)
    q_value = q_value_for_score(model, item.score)
    cluster = item.cluster
    return {
        "pxd_accession": cluster.accession,
        "data_file": cluster.data_file,
        "cluster_rank": cluster.cluster_rank,
        "delta_mass_da": cluster.delta_mass_da,
        "pair_support": cluster.pair_support,
        "unique_spectrum_support": cluster.unique_spectrum_support,
        "median_spectral_similarity": cluster.median_spectral_similarity,
        "support_fraction_of_accepted_pairs": cluster.support_fraction,
        "cluster_sigma_da": cluster.cluster_sigma_da,
        "match_tolerance_da": cluster.match_tolerance_da,
        "classification": cluster.classification,
        "support_tier": cluster.confidence,
        "diagnostic_unimod_candidate_count": cluster.diagnostic_candidate_count,
        "biological_candidate_count": item.candidate_count,
        "unimod_accession": item.unimod_accession,
        "unimod_name": item.unimod_name,
        "unimod_theoretical_delta_mass_da": item.theoretical_mass_da,
        "unimod_residual_da": item.residual_da,
        "candidate_mass_dominance": item.candidate_dominance,
        "runner_up_residual_da": (
            "" if item.runner_up_residual_da is None else item.runner_up_residual_da
        ),
        "residual_margin_da": "" if item.residual_margin_da is None else item.residual_margin_da,
        "evidence_score": item.score,
        "mass_annotation_probability": probability,
        "annotation_local_fdr": 1.0 - probability,
        "annotation_q_value": q_value,
        "probability_model_version": model.model_version,
        "probability_calibration_scope": calibration_scope,
        "source_tsv": cluster.source_tsv,
    }


def _score_cross_validated(
    targets: Sequence[ClusterAnnotation],
    decoy_scores: Sequence[tuple[str, float]],
    *,
    catalog_hash: str,
    biological_catalog_records: int,
    decoy_replicates: int,
    decoy_min_offset_da: float,
    decoy_max_offset_da: float,
    minimum_targets_per_bin: int,
) -> list[dict[str, Any]]:
    accessions = sorted({item.cluster.accession for item in targets})
    output: list[dict[str, Any]] = []
    for accession in accessions:
        train_targets = [item.score for item in targets if item.cluster.accession != accession]
        train_decoys = [score for pxd, score in decoy_scores if pxd != accession]
        train_accessions = sorted(
            {item.cluster.accession for item in targets if item.cluster.accession != accession}
        )
        if not train_targets:
            train_targets = [item.score for item in targets]
            train_decoys = [score for _pxd, score in decoy_scores]
            train_accessions = accessions
        model = fit_probability_model(
            train_targets,
            train_decoys,
            decoy_replicates=decoy_replicates,
            decoy_min_offset_da=decoy_min_offset_da,
            decoy_max_offset_da=decoy_max_offset_da,
            catalog_hash=catalog_hash,
            biological_catalog_records=biological_catalog_records,
            training_accessions=train_accessions,
            minimum_targets_per_bin=minimum_targets_per_bin,
        )
        for item in targets:
            if item.cluster.accession == accession:
                output.append(
                    scored_annotation_row(
                        item,
                        model,
                        calibration_scope="leave-one-accession-out",
                    )
                )
    return output


def build_sdrf_suggestions(
    scored_rows: Sequence[dict[str, Any]],
    *,
    successful_runs_by_accession: dict[str, set[str]],
    probability_threshold: float,
    q_value_threshold: float,
    minimum_candidate_dominance: float,
    minimum_confident_runs: int,
    minimum_confident_run_fraction: float,
    maximum_candidates_per_accession: int,
) -> list[dict[str, Any]]:
    grouped: dict[tuple[str, str, str], list[dict[str, Any]]] = defaultdict(list)
    for row in scored_rows:
        if float(row["mass_annotation_probability"]) < probability_threshold:
            continue
        if float(row["annotation_q_value"]) > q_value_threshold:
            continue
        if float(row["candidate_mass_dominance"]) < minimum_candidate_dominance:
            continue
        key = (
            str(row["pxd_accession"]),
            str(row["unimod_accession"]),
            str(row["unimod_name"]),
        )
        grouped[key].append(row)

    by_accession: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for (accession, unimod_accession, name), rows in grouped.items():
        runs = {str(row["data_file"]) for row in rows}
        total = len(successful_runs_by_accession.get(accession, set()))
        fraction = len(runs) / total if total else 0.0
        if len(runs) < minimum_confident_runs or fraction < minimum_confident_run_fraction:
            continue
        probabilities = [float(row["mass_annotation_probability"]) for row in rows]
        q_values = [float(row["annotation_q_value"]) for row in rows]
        residuals = [abs(float(row["unimod_residual_da"])) for row in rows]
        dominance = [float(row["candidate_mass_dominance"]) for row in rows]
        pair_support = [int(row["pair_support"]) for row in rows]
        similarities = [float(row["median_spectral_similarity"]) for row in rows]
        by_accession[accession].append(
            {
                "pxd_accession": accession,
                "unimod_accession": unimod_accession,
                "unimod_name": name,
                "confident_runs": len(runs),
                "successful_runs": total,
                "confident_run_fraction": fraction,
                "median_probability": statistics.median(probabilities),
                "minimum_probability": min(probabilities),
                "maximum_q_value": max(q_values),
                "median_candidate_mass_dominance": statistics.median(dominance),
                "minimum_candidate_mass_dominance": min(dominance),
                "median_abs_residual_da": statistics.median(residuals),
                "median_pair_support": statistics.median(pair_support),
                "median_spectral_similarity": statistics.median(similarities),
                "sdrf_column": "comment[modification parameters]",
                "sdrf_value": f"NT={name};AC={unimod_accession}",
                "sdrf_status": "review-required",
                "sdrf_warning": (
                    "RAW-derived best mass annotation only; does not establish searched "
                    "modification, fixed/variable status, target residue, localization, or "
                    "PTM identity."
                ),
            }
        )

    output: list[dict[str, Any]] = []
    for _accession, rows in sorted(by_accession.items()):
        rows.sort(
            key=lambda row: (
                -int(row["confident_runs"]),
                -float(row["median_probability"]),
                -float(row["median_candidate_mass_dominance"]),
                float(row["maximum_q_value"]),
                float(row["median_abs_residual_da"]),
                str(row["unimod_accession"]),
            )
        )
        output.extend(rows[:maximum_candidates_per_accession])
    return output


def _write_tsv(path: Path, rows: Sequence[dict[str, Any]], fields: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _scored_fields() -> list[str]:
    return [
        "pxd_accession", "data_file", "cluster_rank", "delta_mass_da", "pair_support",
        "unique_spectrum_support", "median_spectral_similarity",
        "support_fraction_of_accepted_pairs", "cluster_sigma_da", "match_tolerance_da",
        "classification", "support_tier", "diagnostic_unimod_candidate_count",
        "biological_candidate_count", "unimod_accession", "unimod_name",
        "unimod_theoretical_delta_mass_da", "unimod_residual_da",
        "candidate_mass_dominance", "runner_up_residual_da", "residual_margin_da",
        "evidence_score", "mass_annotation_probability", "annotation_local_fdr",
        "annotation_q_value", "probability_model_version", "probability_calibration_scope",
        "source_tsv",
    ]


def _suggestion_fields() -> list[str]:
    return [
        "pxd_accession", "unimod_accession", "unimod_name", "confident_runs",
        "successful_runs", "confident_run_fraction", "median_probability",
        "minimum_probability", "maximum_q_value", "median_candidate_mass_dominance",
        "minimum_candidate_mass_dominance", "median_abs_residual_da",
        "median_pair_support", "median_spectral_similarity", "sdrf_column", "sdrf_value",
        "sdrf_status", "sdrf_warning",
    ]


def calibrate(args: argparse.Namespace) -> None:
    _clusters, candidates = v1.read_mass_shift_tables(args.result_root.resolve())
    targets = biological_cluster_annotations(candidates)
    if args.catalog_source == "openms":
        catalog = v1.catalog_from_openms()
    else:
        catalog = v1.catalog_from_observed_candidates(candidates)
    biological_catalog = [item for item in catalog if item.category == "biological-ptm"]
    if not biological_catalog:
        raise SystemExit("No biological-PTM records are available for calibration.")

    offsets = deterministic_decoy_offsets(
        args.decoy_replicates,
        minimum_offset_da=args.decoy_min_offset_da,
        maximum_offset_da=args.decoy_max_offset_da,
    )
    decoy_scores = shifted_decoy_scores(
        targets,
        biological_catalog,
        offsets=offsets,
        minimum_mass_da=args.minimum_mass_da,
        maximum_mass_da=args.maximum_mass_da,
    )
    digest = v1.catalog_sha256(biological_catalog)
    accessions = sorted({item.cluster.accession for item in targets})

    cross_validated = _score_cross_validated(
        targets,
        decoy_scores,
        catalog_hash=digest,
        biological_catalog_records=len(biological_catalog),
        decoy_replicates=args.decoy_replicates,
        decoy_min_offset_da=args.decoy_min_offset_da,
        decoy_max_offset_da=args.decoy_max_offset_da,
        minimum_targets_per_bin=args.minimum_targets_per_bin,
    )
    full_model = fit_probability_model(
        [item.score for item in targets],
        [score for _accession, score in decoy_scores],
        decoy_replicates=args.decoy_replicates,
        decoy_min_offset_da=args.decoy_min_offset_da,
        decoy_max_offset_da=args.decoy_max_offset_da,
        catalog_hash=digest,
        biological_catalog_records=len(biological_catalog),
        training_accessions=accessions,
        minimum_targets_per_bin=args.minimum_targets_per_bin,
    )
    successful_runs = v1.successful_runs_by_accession(args.result_root.resolve())
    suggestions = build_sdrf_suggestions(
        cross_validated,
        successful_runs_by_accession=successful_runs,
        probability_threshold=args.probability_threshold,
        q_value_threshold=args.q_value_threshold,
        minimum_candidate_dominance=args.minimum_candidate_dominance,
        minimum_confident_runs=args.minimum_confident_runs,
        minimum_confident_run_fraction=args.minimum_confident_run_fraction,
        maximum_candidates_per_accession=args.maximum_candidates_per_accession,
    )

    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    model_path = output / "mass-shift-probability-v2.model.json"
    scored_path = output / "mass-shift-probability-v2.cross-validated-annotations.tsv"
    suggestions_path = output / "mass-shift-probability-v2.sdrf-review.tsv"
    metadata_path = output / "mass-shift-probability-v2.metadata.json"
    model_path.write_text(json.dumps(model_to_json(full_model), indent=2, sort_keys=True) + "\n")
    _write_tsv(scored_path, cross_validated, _scored_fields())
    _write_tsv(suggestions_path, suggestions, _suggestion_fields())
    metadata = {
        "model_version": MODEL_VERSION,
        "result_root": str(args.result_root.resolve()),
        "catalog_source": args.catalog_source,
        "catalog_sha256": digest,
        "biological_catalog_records": len(biological_catalog),
        "target_clusters": len(targets),
        "raw_decoy_hits": len(decoy_scores),
        "decoy_replicates": args.decoy_replicates,
        "decoy_weight": 1.0 / args.decoy_replicates,
        "decoy_min_offset_da": args.decoy_min_offset_da,
        "decoy_max_offset_da": args.decoy_max_offset_da,
        "calibration_scope": "leave-one-accession-out",
        "probability_threshold": args.probability_threshold,
        "q_value_threshold": args.q_value_threshold,
        "minimum_candidate_dominance": args.minimum_candidate_dominance,
        "minimum_confident_runs": args.minimum_confident_runs,
        "minimum_confident_run_fraction": args.minimum_confident_run_fraction,
        "maximum_candidates_per_accession": args.maximum_candidates_per_accession,
        "v1_findings_addressed": [
            "candidate-row multiplicity replaced by one best biological annotation per cluster",
            (
                "four sparse shifted-catalog decoys replaced by dense "
                "cluster-preserving shifted-mass null"
            ),
            "exact/isobaric ambiguity penalized by candidate mass dominance",
        ],
        "sdrf_policy": "review-only; no automatic modification-parameter assertion",
    }
    metadata_path.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")

    print(f"model={model_path}")
    print(f"cross_validated_annotations={scored_path}")
    print(f"sdrf_review={suggestions_path}")
    print(f"target_clusters={len(targets)}")
    print(f"raw_decoy_hits={len(decoy_scores)}")
    print(f"effective_decoy_hits={len(decoy_scores) / args.decoy_replicates:.6f}")
    print(f"sdrf_review_candidates={len(suggestions)}")


def apply_model(args: argparse.Namespace) -> None:
    model = model_from_json(json.loads(args.model.read_text(encoding="utf-8")))
    _clusters, candidates = v1.read_mass_shift_tables(args.result_root.resolve())
    targets = biological_cluster_annotations(candidates)
    scored = [
        scored_annotation_row(item, model, calibration_scope="fixed-model")
        for item in targets
    ]
    successful_runs = v1.successful_runs_by_accession(args.result_root.resolve())
    suggestions = build_sdrf_suggestions(
        scored,
        successful_runs_by_accession=successful_runs,
        probability_threshold=args.probability_threshold,
        q_value_threshold=args.q_value_threshold,
        minimum_candidate_dominance=args.minimum_candidate_dominance,
        minimum_confident_runs=args.minimum_confident_runs,
        minimum_confident_run_fraction=args.minimum_confident_run_fraction,
        maximum_candidates_per_accession=args.maximum_candidates_per_accession,
    )
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    scored_path = output / "mass-shift-probability-v2.scored-annotations.tsv"
    suggestions_path = output / "mass-shift-probability-v2.sdrf-review.tsv"
    _write_tsv(scored_path, scored, _scored_fields())
    _write_tsv(suggestions_path, suggestions, _suggestion_fields())
    print(f"scored_annotations={scored_path}")
    print(f"sdrf_review={suggestions_path}")
    print(f"target_clusters={len(targets)}")
    print(f"sdrf_review_candidates={len(suggestions)}")


def _common_filter_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--probability-threshold", type=float, default=0.95)
    parser.add_argument("--q-value-threshold", type=float, default=0.05)
    parser.add_argument("--minimum-candidate-dominance", type=float, default=0.80)
    parser.add_argument("--minimum-confident-runs", type=int, default=2)
    parser.add_argument("--minimum-confident-run-fraction", type=float, default=0.05)
    parser.add_argument("--maximum-candidates-per-accession", type=int, default=3)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Dense shifted-mass target/decoy calibration for one best biological PTM "
            "annotation per prideQC recurrent cluster."
        )
    )
    subparsers = parser.add_subparsers(dest="command", required=True)
    fit = subparsers.add_parser("calibrate")
    fit.add_argument("--result-root", required=True, type=Path)
    fit.add_argument("--output-dir", required=True, type=Path)
    fit.add_argument("--catalog-source", choices=("openms", "observed"), default="openms")
    fit.add_argument("--decoy-replicates", type=int, default=DEFAULT_DECOY_REPLICATES)
    fit.add_argument("--decoy-min-offset-da", type=float, default=DEFAULT_DECOY_MIN_OFFSET_DA)
    fit.add_argument("--decoy-max-offset-da", type=float, default=DEFAULT_DECOY_MAX_OFFSET_DA)
    fit.add_argument("--minimum-mass-da", type=float, default=DEFAULT_MINIMUM_MASS_DA)
    fit.add_argument("--maximum-mass-da", type=float, default=DEFAULT_MAXIMUM_MASS_DA)
    fit.add_argument("--minimum-targets-per-bin", type=int, default=25)
    _common_filter_arguments(fit)

    apply_parser = subparsers.add_parser("apply")
    apply_parser.add_argument("--result-root", required=True, type=Path)
    apply_parser.add_argument("--model", required=True, type=Path)
    apply_parser.add_argument("--output-dir", required=True, type=Path)
    _common_filter_arguments(apply_parser)
    return parser


def _validate_args(args: argparse.Namespace) -> None:
    for name in ("probability_threshold", "q_value_threshold", "minimum_candidate_dominance"):
        value = float(getattr(args, name))
        if not 0.0 <= value <= 1.0:
            raise SystemExit(f"--{name.replace('_', '-')} must be in [0, 1].")
    if args.minimum_confident_runs < 1:
        raise SystemExit("--minimum-confident-runs must be >= 1.")
    if not 0.0 <= args.minimum_confident_run_fraction <= 1.0:
        raise SystemExit("--minimum-confident-run-fraction must be in [0, 1].")
    if args.maximum_candidates_per_accession < 1:
        raise SystemExit("--maximum-candidates-per-accession must be >= 1.")
    if args.command == "calibrate":
        if args.decoy_replicates < 16:
            raise SystemExit("--decoy-replicates must be >= 16 for a dense null.")
        if not 0.02 < args.decoy_min_offset_da < args.decoy_max_offset_da:
            raise SystemExit("decoy offsets must satisfy 0.02 < minimum < maximum.")
        if args.minimum_mass_da >= args.maximum_mass_da:
            raise SystemExit("--minimum-mass-da must be smaller than --maximum-mass-da.")
        if args.minimum_targets_per_bin < 1:
            raise SystemExit("--minimum-targets-per-bin must be >= 1.")


def main() -> None:
    args = build_parser().parse_args()
    _validate_args(args)
    if args.command == "calibrate":
        calibrate(args)
    else:
        apply_model(args)


if __name__ == "__main__":
    main()
