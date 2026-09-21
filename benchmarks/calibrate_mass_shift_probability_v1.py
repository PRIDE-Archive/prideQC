#!/usr/bin/env python3
"""Calibrate identification-free PTM-candidate confidence from frozen mass-shift TSVs.

This is an offline calibration layer for the v22.1 recurrent mass-shift scout. It
never changes spectrum pairing, clustering, UniMod matching, or the frozen scout
output. Instead it scores existing biological-PTM candidate rows, estimates an
empirical target/decoy local FDR, and produces a conservative accession-level
shortlist for SDRF review.

Probability semantics
---------------------
``ptm_candidate_probability`` is an empirical calibration quantity: the estimated
probability that a biological-PTM candidate match is not explained by chance
matching to a mass-density-preserving shifted (decoy) modification catalogue,
given the recurrent-cluster evidence used by prideQC. It is *not* a peptide,
site-localization, or PTM-identity posterior probability.

The SDRF output is deliberately a review file. SDRF ``comment[modification
parameters]`` describes modifications searched; prideQC cannot infer search
configuration, fixed/variable status, target residues, or localization from RAW
mass-shift evidence alone.
"""

from __future__ import annotations

import argparse
import bisect
import csv
import hashlib
import json
import math
import re
import statistics
from collections import Counter, defaultdict
from collections.abc import Iterable, Sequence
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

MODEL_VERSION = "mass-shift-probability-v1.0"
DEFAULT_DECOY_OFFSETS_DA = (0.271828, 0.618034, 0.913271, 1.414214)
DEFAULT_MINIMUM_MASS_DA = 0.5
DEFAULT_MAXIMUM_MASS_DA = 500.0


@dataclass(frozen=True, slots=True)
class CatalogRecord:
    accession: str
    name: str
    mass_da: float
    category: str


@dataclass(frozen=True, slots=True)
class ClusterEvidence:
    accession: str
    data_file: str
    cluster_key: str
    cluster_rank: int
    delta_mass_da: float
    cluster_sigma_da: float
    cluster_min_da: float
    cluster_max_da: float
    pair_support: int
    unique_spectrum_support: int
    median_spectral_similarity: float
    classification: str
    confidence: str
    match_tolerance_da: float
    support_fraction: float
    diagnostic_candidate_count: int
    source_tsv: str


@dataclass(frozen=True, slots=True)
class CandidateObservation:
    cluster: ClusterEvidence
    unimod_accession: str
    unimod_name: str
    theoretical_mass_da: float
    residual_da: float
    candidate_category: str
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
    decoy_offsets_da: tuple[float, ...]
    decoy_multiplicity: int
    catalog_sha256: str
    biological_catalog_records: int
    training_accessions: tuple[str, ...]
    training_target_observations: int
    training_decoy_observations: int
    calibration_bins: tuple[CalibrationBin, ...]
    q_thresholds: tuple[tuple[float, float], ...]


def _float(value: str | None, default: float = 0.0) -> float:
    if value in (None, ""):
        return default
    try:
        result = float(value)
    except (TypeError, ValueError):
        return default
    return result if math.isfinite(result) else default


def _int(value: str | None, default: int = 0) -> int:
    if value in (None, ""):
        return default
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return default


def _candidate_category(name: str, source_classification: str) -> str:
    folded_name = name.casefold()
    folded_source = source_classification.casefold()
    if "decoy" in folded_name or "decoy" in folded_source:
        return "decoy"
    if "substitution" in folded_name or "substitution" in folded_source or "->" in name:
        return "amino-acid-substitution"
    if any(
        token in folded_source
        for token in ("post", "co-translational", "pre-translational", "glycosyl")
    ):
        return "biological-ptm"
    if any(token in folded_source for token in ("chemical", "artefact", "artifact")):
        return "sample-prep-or-artifact"
    return "other-modification"


def evidence_score(
    *,
    pair_support: int,
    unique_spectrum_support: int,
    median_spectral_similarity: float,
    support_fraction: float,
    cluster_sigma_da: float,
    match_tolerance_da: float,
    residual_da: float,
) -> float:
    """Transparent monotonic evidence score used only as a calibration axis.

    The absolute scale has no probability interpretation. Target/decoy calibration
    maps this score to local FDR / probability. Candidate ambiguity is kept out of
    the score so target and decoy observations are treated symmetrically; ambiguity
    is applied later as a conservative SDRF-filter guard.
    """

    tolerance = max(abs(match_tolerance_da), 1e-12)
    residual_ratio = min(abs(residual_da) / tolerance, 5.0)
    sigma_ratio = min(abs(cluster_sigma_da) / tolerance, 5.0)
    similarity = min(max(median_spectral_similarity, 0.0), 1.0)
    fraction = min(max(support_fraction, 0.0), 1.0)
    return (
        1.25 * math.log1p(max(pair_support, 0))
        + 1.00 * math.log1p(max(unique_spectrum_support, 0))
        + 3.00 * similarity
        + 0.75 * math.log1p(1000.0 * fraction)
        - 1.50 * residual_ratio
        - 0.25 * sigma_ratio
    )


def _accession_from_path(path: Path) -> str:
    for part in reversed(path.parts):
        if re.fullmatch(r"PXD\d+", part, flags=re.IGNORECASE):
            return part.upper()
    return "UNKNOWN"


def read_mass_shift_tables(
    result_root: Path,
) -> tuple[list[ClusterEvidence], list[CandidateObservation]]:
    """Read reported clusters once while preserving candidate-level rows."""

    clusters: list[ClusterEvidence] = []
    candidates: list[CandidateObservation] = []
    for path in sorted(result_root.glob("**/*.mass-shifts.tsv")):
        accession = _accession_from_path(path)
        with path.open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        grouped: dict[int, list[dict[str, str]]] = defaultdict(list)
        for row in rows:
            rank = _int(row.get("cluster_rank"), -1)
            if rank > 0:
                grouped[rank].append(row)
        for rank, group in sorted(grouped.items()):
            first = group[0]
            tolerance = _float(first.get("match_tolerance_da"), 0.02)
            cluster = ClusterEvidence(
                accession=accession,
                data_file=first.get("data_file", path.name.removesuffix(".mass-shifts.tsv")),
                cluster_key=f"{path}:{rank}",
                cluster_rank=rank,
                delta_mass_da=_float(first.get("delta_mass_da")),
                cluster_sigma_da=_float(first.get("cluster_sigma_da")),
                cluster_min_da=_float(first.get("cluster_min_da")),
                cluster_max_da=_float(first.get("cluster_max_da")),
                pair_support=_int(first.get("pair_support")),
                unique_spectrum_support=_int(first.get("unique_spectrum_support")),
                median_spectral_similarity=_float(first.get("median_spectral_similarity")),
                classification=first.get("classification", ""),
                confidence=first.get("confidence", ""),
                match_tolerance_da=tolerance,
                support_fraction=_float(first.get("support_fraction_of_accepted_pairs")),
                diagnostic_candidate_count=_int(first.get("diagnostic_unimod_candidate_count")),
                source_tsv=str(path),
            )
            clusters.append(cluster)
            for row in group:
                accession_value = (row.get("unimod_accession") or "").strip()
                name = (row.get("unimod_name") or "").strip()
                category = (row.get("unimod_candidate_category") or "").strip()
                theoretical = _float(row.get("unimod_theoretical_delta_mass_da"), math.nan)
                residual = _float(row.get("unimod_residual_da"), math.nan)
                if (
                    not accession_value
                    or not name
                    or not math.isfinite(theoretical)
                    or not math.isfinite(residual)
                ):
                    continue
                if not category:
                    category = _candidate_category(
                        name,
                        (row.get("unimod_source_classification") or "").strip(),
                    )
                candidates.append(
                    CandidateObservation(
                        cluster=cluster,
                        unimod_accession=accession_value,
                        unimod_name=name,
                        theoretical_mass_da=theoretical,
                        residual_da=residual,
                        candidate_category=category,
                        score=evidence_score(
                            pair_support=cluster.pair_support,
                            unique_spectrum_support=cluster.unique_spectrum_support,
                            median_spectral_similarity=cluster.median_spectral_similarity,
                            support_fraction=cluster.support_fraction,
                            cluster_sigma_da=cluster.cluster_sigma_da,
                            match_tolerance_da=cluster.match_tolerance_da,
                            residual_da=residual,
                        ),
                    )
                )
    return clusters, candidates




def successful_runs_by_accession(result_root: Path) -> dict[str, set[str]]:
    """Count every successful-run mass-shift table, including header-only TSVs."""

    runs: dict[str, set[str]] = defaultdict(set)
    for path in sorted(result_root.glob("**/*.mass-shifts.tsv")):
        accession = _accession_from_path(path)
        runs[accession].add(path.name.removesuffix(".mass-shifts.tsv"))
    return dict(runs)


def catalog_from_observed_candidates(
    candidates: Iterable[CandidateObservation],
) -> list[CatalogRecord]:
    """Fallback catalogue for tests/offline debugging; OpenMS is preferred."""

    unique: dict[tuple[str, str, float], CatalogRecord] = {}
    for item in candidates:
        key = (item.unimod_accession, item.unimod_name, round(abs(item.theoretical_mass_da), 9))
        unique[key] = CatalogRecord(
            item.unimod_accession,
            item.unimod_name,
            abs(item.theoretical_mass_da),
            item.candidate_category,
        )
    return sorted(unique.values(), key=lambda item: (item.mass_da, item.accession, item.name))


def catalog_from_openms() -> list[CatalogRecord]:
    """Load the full pinned OpenMS/UniMod catalogue used by prideQC."""

    from prideqc.mass_shift import load_openms_modifications

    result = []
    for item in load_openms_modifications():
        result.append(
            CatalogRecord(
                accession=item.accession,
                name=item.name,
                mass_da=abs(item.delta_mass_da),
                category=_candidate_category(item.name, item.source_classification),
            )
        )
    return result


def catalog_sha256(catalog: Sequence[CatalogRecord]) -> str:
    lines = [
        f"{item.accession}\t{item.name}\t{item.mass_da:.9f}\t{item.category}\n"
        for item in sorted(catalog, key=lambda x: (x.mass_da, x.accession, x.name))
    ]
    return hashlib.sha256("".join(lines).encode()).hexdigest()


def _wrapped_shift(mass: float, offset: float, minimum: float, maximum: float) -> float:
    width = maximum - minimum
    if width <= 0:
        raise ValueError("Decoy mass range must have positive width.")
    return minimum + ((mass - minimum + offset) % width)


def decoy_catalog_masses(
    catalog: Sequence[CatalogRecord],
    *,
    offsets: Sequence[float],
    minimum_mass_da: float,
    maximum_mass_da: float,
) -> list[float]:
    biological = [
        item.mass_da
        for item in catalog
        if item.category == "biological-ptm" and minimum_mass_da <= item.mass_da <= maximum_mass_da
    ]
    masses = [
        _wrapped_shift(mass, offset, minimum_mass_da, maximum_mass_da)
        for offset in offsets
        for mass in biological
    ]
    return sorted(masses)


def _candidate_decoy_scores(
    clusters: Sequence[ClusterEvidence],
    decoy_masses: Sequence[float],
) -> list[tuple[str, float]]:
    """Return accession/score for every decoy catalogue match."""

    result: list[tuple[str, float]] = []
    for cluster in clusters:
        tolerance = cluster.match_tolerance_da
        lo = bisect.bisect_left(decoy_masses, cluster.delta_mass_da - tolerance)
        hi = bisect.bisect_right(decoy_masses, cluster.delta_mass_da + tolerance)
        for mass in decoy_masses[lo:hi]:
            residual = cluster.delta_mass_da - mass
            result.append(
                (
                    cluster.accession,
                    evidence_score(
                        pair_support=cluster.pair_support,
                        unique_spectrum_support=cluster.unique_spectrum_support,
                        median_spectral_similarity=cluster.median_spectral_similarity,
                        support_fraction=cluster.support_fraction,
                        cluster_sigma_da=cluster.cluster_sigma_da,
                        match_tolerance_da=tolerance,
                        residual_da=residual,
                    ),
                )
            )
    return result


def _pava_non_decreasing(values: Sequence[float], weights: Sequence[float]) -> list[float]:
    """Weighted PAVA: return non-decreasing fitted values."""

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
    decoy_offsets: Sequence[float],
    catalog_hash: str,
    biological_catalog_records: int,
    training_accessions: Sequence[str],
    minimum_targets_per_bin: int = 25,
) -> ProbabilityModel:
    decoy_multiplicity = len(decoy_offsets)
    if decoy_multiplicity <= 0:
        raise ValueError("At least one decoy offset is required.")
    if not target_scores:
        raise ValueError("Cannot fit probability model without biological PTM target candidates.")

    targets = sorted(float(x) for x in target_scores)
    decoys = sorted(float(x) for x in decoy_scores)

    # Adaptive non-overlapping bins are defined from unique target-score levels.
    # Score bins are stored high-to-low so local FDR should be non-decreasing
    # down the list. Repeated scores are never split across adjacent bins.
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
            previous_minimum, previous_maximum = bin_ranges[-1]
            del previous_minimum
            bin_ranges[-1] = (score_counts[-1][0], previous_maximum)
        else:
            bin_ranges.append((score_counts[-1][0], bin_high))

    raw_fdr = []
    weights = []
    counts: list[tuple[float, float]] = []
    decoy_scale = 1.0 / decoy_multiplicity
    for minimum, maximum in bin_ranges:
        target_count = sum(minimum <= score <= maximum for score in targets)
        decoy_count = sum(minimum <= score <= maximum for score in decoys)
        weighted_decoys = decoy_count * decoy_scale
        # Half-decoy pseudo-count avoids spuriously perfect probability in sparse
        # high-score bins while becoming negligible in well-populated bins.
        local = min(1.0, (weighted_decoys + 0.5 * decoy_scale) / max(target_count, 1))
        raw_fdr.append(local)
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

    # Tail-FDR thresholds for q-values. Higher scores are better. Reverse
    # cumulative minima make q non-decreasing as score decreases.
    thresholds = sorted(set(targets + decoys), reverse=True)
    raw_q: list[tuple[float, float]] = []
    for threshold in thresholds:
        target_tail = len(targets) - bisect.bisect_left(targets, threshold)
        decoy_tail = len(decoys) - bisect.bisect_left(decoys, threshold)
        weighted_decoy_tail = decoy_tail * decoy_scale
        fdr = min(
            1.0,
            (weighted_decoy_tail + decoy_scale) / max(target_tail, 1),
        )
        raw_q.append((threshold, fdr))
    running = 1.0
    q_values = [0.0] * len(raw_q)
    for index in range(len(raw_q) - 1, -1, -1):
        running = min(running, raw_q[index][1])
        q_values[index] = running
    q_thresholds = tuple(
        (raw_q[index][0], q_values[index]) for index in range(len(raw_q))
    )

    return ProbabilityModel(
        model_version=MODEL_VERSION,
        probability_semantics=(
            "Empirical probability that a biological-PTM mass candidate is not explained by "
            "chance matching to a shifted UniMod decoy catalogue, conditional on prideQC "
            "recurrent-cluster evidence; not peptide/site/PTM-identity probability."
        ),
        score_formula=(
            "1.25*log1p(pair_support) + log1p(unique_spectrum_support) + "
            "3*median_spectral_similarity + 0.75*log1p(1000*support_fraction) - "
            "1.5*abs(candidate_residual)/match_tolerance - "
            "0.25*cluster_sigma/match_tolerance"
        ),
        decoy_offsets_da=tuple(float(x) for x in decoy_offsets),
        decoy_multiplicity=decoy_multiplicity,
        catalog_sha256=catalog_hash,
        biological_catalog_records=biological_catalog_records,
        training_accessions=tuple(sorted(training_accessions)),
        training_target_observations=len(targets),
        training_decoy_observations=len(decoys),
        calibration_bins=bins,
        q_thresholds=q_thresholds,
    )


def probability_for_score(model: ProbabilityModel, score: float) -> float:
    bins = model.calibration_bins
    if not bins:
        return 0.0
    for item in bins:
        if item.minimum_score <= score <= item.maximum_score:
            return item.probability
    if score > bins[0].maximum_score:
        return bins[0].probability
    return bins[-1].probability


def q_value_for_score(model: ProbabilityModel, score: float) -> float:
    thresholds = model.q_thresholds
    if not thresholds:
        return 1.0
    # Thresholds are high-to-low. Use the least stringent threshold not above
    # the observation so q cannot improve merely by interpolation.
    for threshold, q_value in thresholds:
        if score >= threshold:
            return q_value
    return thresholds[-1][1]


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
        decoy_offsets_da=tuple(float(x) for x in payload["decoy_offsets_da"]),
        decoy_multiplicity=int(payload["decoy_multiplicity"]),
        catalog_sha256=str(payload["catalog_sha256"]),
        biological_catalog_records=int(payload["biological_catalog_records"]),
        training_accessions=tuple(str(x) for x in payload["training_accessions"]),
        training_target_observations=int(payload["training_target_observations"]),
        training_decoy_observations=int(payload["training_decoy_observations"]),
        calibration_bins=tuple(CalibrationBin(**item) for item in payload["calibration_bins"]),
        q_thresholds=tuple((float(x[0]), float(x[1])) for x in payload["q_thresholds"]),
    )


def _biological_targets(candidates: Iterable[CandidateObservation]) -> list[CandidateObservation]:
    return [
        item
        for item in candidates
        if item.candidate_category == "biological-ptm"
        and item.cluster.classification == "putative-ptm"
    ]


def _score_candidates_cross_validated(
    targets: Sequence[CandidateObservation],
    decoy_scores: Sequence[tuple[str, float]],
    *,
    catalog_hash: str,
    biological_catalog_records: int,
    decoy_offsets: Sequence[float],
    minimum_targets_per_bin: int,
) -> list[dict[str, Any]]:
    accessions = sorted({item.cluster.accession for item in targets})
    output: list[dict[str, Any]] = []
    for accession in accessions:
        train_targets = [item.score for item in targets if item.cluster.accession != accession]
        train_decoys = [score for pxd, score in decoy_scores if pxd != accession]
        train_accessions = sorted(
            {
                item.cluster.accession
                for item in targets
                if item.cluster.accession != accession
            }
        )
        if not train_targets:
            train_targets = [item.score for item in targets]
            train_decoys = [score for _pxd, score in decoy_scores]
            train_accessions = accessions
        model = fit_probability_model(
            train_targets,
            train_decoys,
            decoy_offsets=decoy_offsets,
            catalog_hash=catalog_hash,
            biological_catalog_records=biological_catalog_records,
            training_accessions=train_accessions,
            minimum_targets_per_bin=minimum_targets_per_bin,
        )
        for item in targets:
            if item.cluster.accession != accession:
                continue
            output.append(
                scored_candidate_row(
                    item,
                    model,
                    calibration_scope="leave-one-accession-out",
                )
            )
    return output


def scored_candidate_row(
    item: CandidateObservation,
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
        "unimod_accession": item.unimod_accession,
        "unimod_name": item.unimod_name,
        "unimod_theoretical_delta_mass_da": item.theoretical_mass_da,
        "unimod_residual_da": item.residual_da,
        "unimod_candidate_category": item.candidate_category,
        "evidence_score": item.score,
        "ptm_candidate_probability": probability,
        "candidate_local_fdr": 1.0 - probability,
        "candidate_q_value": q_value,
        "probability_model_version": model.model_version,
        "probability_calibration_scope": calibration_scope,
        "source_tsv": cluster.source_tsv,
    }


def build_sdrf_suggestions(
    scored_rows: Sequence[dict[str, Any]],
    *,
    successful_runs_by_accession: dict[str, set[str]] | None = None,
    probability_threshold: float,
    q_value_threshold: float,
    maximum_candidate_ambiguity: int,
    minimum_confident_runs: int,
    minimum_confident_run_fraction: float,
    maximum_candidates_per_accession: int,
) -> list[dict[str, Any]]:
    """Collapse confident run-level evidence into a small review shortlist."""

    successful_runs: dict[str, set[str]] = defaultdict(set)
    if successful_runs_by_accession is not None:
        for accession, runs in successful_runs_by_accession.items():
            successful_runs[accession].update(runs)
    else:
        for row in scored_rows:
            successful_runs[str(row["pxd_accession"])].add(str(row["data_file"]))

    grouped: dict[tuple[str, str, str], list[dict[str, Any]]] = defaultdict(list)
    for row in scored_rows:
        if float(row["ptm_candidate_probability"]) < probability_threshold:
            continue
        if float(row["candidate_q_value"]) > q_value_threshold:
            continue
        if int(row["diagnostic_unimod_candidate_count"]) > maximum_candidate_ambiguity:
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
        total = len(successful_runs.get(accession, set()))
        fraction = len(runs) / total if total else 0.0
        if len(runs) < minimum_confident_runs or fraction < minimum_confident_run_fraction:
            continue
        probabilities = [float(row["ptm_candidate_probability"]) for row in rows]
        q_values = [float(row["candidate_q_value"]) for row in rows]
        residuals = [abs(float(row["unimod_residual_da"])) for row in rows]
        pair_support = [int(row["pair_support"]) for row in rows]
        similarities = [float(row["median_spectral_similarity"]) for row in rows]
        by_accession[accession].append(
            {
                "pxd_accession": accession,
                "unimod_accession": unimod_accession,
                "unimod_name": name,
                "confident_runs": len(runs),
                "successful_runs_with_scored_biological_candidates": total,
                "confident_run_fraction": fraction,
                "median_probability": statistics.median(probabilities),
                "minimum_probability": min(probabilities),
                "maximum_q_value": max(q_values),
                "median_abs_residual_da": statistics.median(residuals),
                "median_pair_support": statistics.median(pair_support),
                "median_spectral_similarity": statistics.median(similarities),
                "sdrf_column": "comment[modification parameters]",
                "sdrf_value": f"NT={name};AC={unimod_accession}",
                "sdrf_status": "review-required",
                "sdrf_warning": (
                    "RAW-derived candidate only; does not establish searched modification, "
                    "fixed/variable status, target residue, or localization."
                ),
            }
        )

    output: list[dict[str, Any]] = []
    for _accession, rows in sorted(by_accession.items()):
        rows.sort(
            key=lambda row: (
                -int(row["confident_runs"]),
                -float(row["median_probability"]),
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
        "pxd_accession",
        "data_file",
        "cluster_rank",
        "delta_mass_da",
        "pair_support",
        "unique_spectrum_support",
        "median_spectral_similarity",
        "support_fraction_of_accepted_pairs",
        "cluster_sigma_da",
        "match_tolerance_da",
        "classification",
        "support_tier",
        "diagnostic_unimod_candidate_count",
        "unimod_accession",
        "unimod_name",
        "unimod_theoretical_delta_mass_da",
        "unimod_residual_da",
        "unimod_candidate_category",
        "evidence_score",
        "ptm_candidate_probability",
        "candidate_local_fdr",
        "candidate_q_value",
        "probability_model_version",
        "probability_calibration_scope",
        "source_tsv",
    ]


def _suggestion_fields() -> list[str]:
    return [
        "pxd_accession",
        "unimod_accession",
        "unimod_name",
        "confident_runs",
        "successful_runs_with_scored_biological_candidates",
        "confident_run_fraction",
        "median_probability",
        "minimum_probability",
        "maximum_q_value",
        "median_abs_residual_da",
        "median_pair_support",
        "median_spectral_similarity",
        "sdrf_column",
        "sdrf_value",
        "sdrf_status",
        "sdrf_warning",
    ]


def calibrate(args: argparse.Namespace) -> None:
    clusters, candidates = read_mass_shift_tables(args.result_root.resolve())
    targets = _biological_targets(candidates)
    if args.catalog_source == "openms":
        catalog = catalog_from_openms()
    else:
        catalog = catalog_from_observed_candidates(candidates)
    biological_catalog = [item for item in catalog if item.category == "biological-ptm"]
    if not biological_catalog:
        raise SystemExit("No biological-PTM records are available for target/decoy calibration.")

    offsets = tuple(float(x) for x in args.decoy_offset)
    decoys = decoy_catalog_masses(
        biological_catalog,
        offsets=offsets,
        minimum_mass_da=args.minimum_mass_da,
        maximum_mass_da=args.maximum_mass_da,
    )
    decoy_scores = _candidate_decoy_scores(clusters, decoys)
    digest = catalog_sha256(biological_catalog)
    accessions = sorted({item.cluster.accession for item in targets})

    cross_validated = _score_candidates_cross_validated(
        targets,
        decoy_scores,
        catalog_hash=digest,
        biological_catalog_records=len(biological_catalog),
        decoy_offsets=offsets,
        minimum_targets_per_bin=args.minimum_targets_per_bin,
    )
    full_model = fit_probability_model(
        [item.score for item in targets],
        [score for _accession, score in decoy_scores],
        decoy_offsets=offsets,
        catalog_hash=digest,
        biological_catalog_records=len(biological_catalog),
        training_accessions=accessions,
        minimum_targets_per_bin=args.minimum_targets_per_bin,
    )

    successful_runs = successful_runs_by_accession(args.result_root.resolve())
    suggestions = build_sdrf_suggestions(
        cross_validated,
        successful_runs_by_accession=successful_runs,
        probability_threshold=args.probability_threshold,
        q_value_threshold=args.q_value_threshold,
        maximum_candidate_ambiguity=args.maximum_candidate_ambiguity,
        minimum_confident_runs=args.minimum_confident_runs,
        minimum_confident_run_fraction=args.minimum_confident_run_fraction,
        maximum_candidates_per_accession=args.maximum_candidates_per_accession,
    )

    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    model_path = output / "mass-shift-probability-v1.model.json"
    model_path.write_text(json.dumps(model_to_json(full_model), indent=2, sort_keys=True) + "\n")
    scored_path = output / "mass-shift-probability-v1.cross-validated-candidates.tsv"
    suggestions_path = output / "mass-shift-probability-v1.sdrf-review.tsv"
    _write_tsv(scored_path, cross_validated, _scored_fields())
    _write_tsv(suggestions_path, suggestions, _suggestion_fields())

    metadata = {
        "model_version": MODEL_VERSION,
        "result_root": str(args.result_root.resolve()),
        "catalog_source": args.catalog_source,
        "catalog_sha256": digest,
        "biological_catalog_records": len(biological_catalog),
        "reported_clusters": len(clusters),
        "biological_target_candidates": len(targets),
        "decoy_observations": len(decoy_scores),
        "decoy_offsets_da": offsets,
        "calibration_scope": "leave-one-accession-out",
        "probability_threshold": args.probability_threshold,
        "q_value_threshold": args.q_value_threshold,
        "maximum_candidate_ambiguity": args.maximum_candidate_ambiguity,
        "minimum_confident_runs": args.minimum_confident_runs,
        "minimum_confident_run_fraction": args.minimum_confident_run_fraction,
        "maximum_candidates_per_accession": args.maximum_candidates_per_accession,
        "sdrf_policy": (
            "review-only; never auto-write search modification parameters from RAW evidence"
        ),
    }
    (output / "mass-shift-probability-v1.metadata.json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    print(f"model={model_path}")
    print(f"cross_validated_candidates={scored_path}")
    print(f"sdrf_review={suggestions_path}")
    print(f"reported_clusters={len(clusters)}")
    print(f"biological_target_candidates={len(targets)}")
    print(f"decoy_observations={len(decoy_scores)}")
    print(f"sdrf_review_candidates={len(suggestions)}")


def apply_model(args: argparse.Namespace) -> None:
    payload = json.loads(args.model.read_text(encoding="utf-8"))
    model = model_from_json(payload)
    clusters, candidates = read_mass_shift_tables(args.result_root.resolve())
    targets = _biological_targets(candidates)
    scored = [
        scored_candidate_row(item, model, calibration_scope="fixed-model") for item in targets
    ]
    successful_runs = successful_runs_by_accession(args.result_root.resolve())
    suggestions = build_sdrf_suggestions(
        scored,
        successful_runs_by_accession=successful_runs,
        probability_threshold=args.probability_threshold,
        q_value_threshold=args.q_value_threshold,
        maximum_candidate_ambiguity=args.maximum_candidate_ambiguity,
        minimum_confident_runs=args.minimum_confident_runs,
        minimum_confident_run_fraction=args.minimum_confident_run_fraction,
        maximum_candidates_per_accession=args.maximum_candidates_per_accession,
    )
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    scored_path = output / "mass-shift-probability-v1.scored-candidates.tsv"
    suggestions_path = output / "mass-shift-probability-v1.sdrf-review.tsv"
    _write_tsv(scored_path, scored, _scored_fields())
    _write_tsv(suggestions_path, suggestions, _suggestion_fields())
    print(f"scored_candidates={scored_path}")
    print(f"sdrf_review={suggestions_path}")
    print(f"biological_target_candidates={len(targets)}")
    print(f"sdrf_review_candidates={len(suggestions)}")


def _common_filter_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--probability-threshold", type=float, default=0.95)
    parser.add_argument("--q-value-threshold", type=float, default=0.05)
    parser.add_argument("--maximum-candidate-ambiguity", type=int, default=3)
    parser.add_argument("--minimum-confident-runs", type=int, default=2)
    parser.add_argument("--minimum-confident-run-fraction", type=float, default=0.05)
    parser.add_argument("--maximum-candidates-per-accession", type=int, default=5)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Target/decoy calibration and conservative SDRF review shortlist for prideQC "
            "recurrent biological-PTM mass candidates."
        )
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    fit = subparsers.add_parser("calibrate", help="Fit and cross-validate a probability model.")
    fit.add_argument("--result-root", required=True, type=Path)
    fit.add_argument("--output-dir", required=True, type=Path)
    fit.add_argument("--catalog-source", choices=("openms", "observed"), default="openms")
    fit.add_argument(
        "--decoy-offset",
        action="append",
        type=float,
        default=None,
        help="Repeatable deterministic decoy mass shift. Defaults to four frozen offsets.",
    )
    fit.add_argument("--minimum-mass-da", type=float, default=DEFAULT_MINIMUM_MASS_DA)
    fit.add_argument("--maximum-mass-da", type=float, default=DEFAULT_MAXIMUM_MASS_DA)
    fit.add_argument("--minimum-targets-per-bin", type=int, default=25)
    _common_filter_arguments(fit)

    apply_parser = subparsers.add_parser("apply", help="Apply a frozen probability model.")
    apply_parser.add_argument("--result-root", required=True, type=Path)
    apply_parser.add_argument("--model", required=True, type=Path)
    apply_parser.add_argument("--output-dir", required=True, type=Path)
    _common_filter_arguments(apply_parser)
    return parser


def _validate_args(args: argparse.Namespace) -> None:
    for name in ("probability_threshold", "q_value_threshold"):
        value = float(getattr(args, name))
        if not 0.0 <= value <= 1.0:
            raise SystemExit(f"--{name.replace('_', '-')} must be in [0, 1].")
    if args.maximum_candidate_ambiguity < 1:
        raise SystemExit("--maximum-candidate-ambiguity must be >= 1.")
    if args.minimum_confident_runs < 1:
        raise SystemExit("--minimum-confident-runs must be >= 1.")
    if not 0.0 <= args.minimum_confident_run_fraction <= 1.0:
        raise SystemExit("--minimum-confident-run-fraction must be in [0, 1].")
    if args.maximum_candidates_per_accession < 1:
        raise SystemExit("--maximum-candidates-per-accession must be >= 1.")
    if args.command == "calibrate":
        if args.minimum_targets_per_bin < 1:
            raise SystemExit("--minimum-targets-per-bin must be >= 1.")
        if args.minimum_mass_da >= args.maximum_mass_da:
            raise SystemExit("--minimum-mass-da must be smaller than --maximum-mass-da.")
        if args.decoy_offset is None:
            args.decoy_offset = list(DEFAULT_DECOY_OFFSETS_DA)
        if not args.decoy_offset or any(abs(x) <= 0.02 for x in args.decoy_offset):
            raise SystemExit(
                "Decoy offsets must be nonzero and safely exceed the 0.02-Da match window."
            )


def main() -> None:
    args = build_parser().parse_args()
    _validate_args(args)
    if args.command == "calibrate":
        calibrate(args)
    else:
        apply_model(args)


if __name__ == "__main__":
    main()
