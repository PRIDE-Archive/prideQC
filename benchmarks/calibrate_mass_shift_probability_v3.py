#!/usr/bin/env python3
"""Accession-level PTM mass-family enrichment for SDRF review.

v3 deliberately changes the statistical question from the v1/v2 prototypes.
It does *not* ask whether an individual recurrent cluster or UniMod row is real.
Instead it:

1. aggregates reported ``putative-ptm`` clusters across RAW files into one
   accession-level mass family;
2. compares that family's run prevalence in the accession with prevalence of
   the same mass window across all other accessions;
3. estimates an exact Beta-posterior superiority probability;
4. annotates the robust cross-run family centroid against biological UniMod
   entries and quantifies mass-only candidate dominance; and
5. emits a conservative SDRF review table that still requires independent
   study/publication evidence before a modification can be accepted.

``study_enrichment_probability`` means the posterior probability that the
run-level prevalence of this mass family is higher in the accession than in
its calibration background. It is not a peptide-identification, localization,
search-parameter, or PTM-identity probability.
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
from collections import defaultdict
from collections.abc import Iterable, Sequence
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

MODEL_VERSION = "mass-shift-probability-v3.0"
DEFAULT_FAMILY_TOLERANCE_DA = 0.02
DEFAULT_ANNOTATION_WINDOW_DA = 0.02
DEFAULT_SYSTEMATIC_MASS_FLOOR_DA = 0.0015
DEFAULT_BETA_PRIOR_ALPHA = 1
DEFAULT_BETA_PRIOR_BETA = 1


@dataclass(frozen=True, slots=True)
class MassFamily:
    accession: str
    family_id: str
    median_mass_da: float
    mass_mad_da: float
    run_count: int
    successful_runs: int
    run_prevalence: float
    median_cluster_rank: float
    median_pair_support: float
    median_unique_spectrum_support: float
    median_spectral_similarity: float
    median_support_fraction: float
    median_cluster_sigma_da: float
    run_names: tuple[str, ...]
    source_cluster_keys: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class FamilyAnnotation:
    family: MassFamily
    unimod_accession: str
    unimod_name: str
    theoretical_mass_da: float
    residual_da: float
    candidate_count: int
    candidate_mass_dominance: float
    runner_up_residual_da: float | None
    residual_margin_da: float | None
    mass_likelihood_sigma_da: float


@dataclass(frozen=True, slots=True)
class StudyEnrichment:
    accession_hits: int
    accession_total: int
    background_hits: int
    background_total: int
    accession_prevalence: float
    background_prevalence: float
    prevalence_difference: float
    prevalence_ratio: float
    study_enrichment_probability: float
    study_enrichment_local_fdr: float


def _median(values: Iterable[float], default: float = 0.0) -> float:
    rows = [float(value) for value in values]
    return statistics.median(rows) if rows else default


def _mad(values: Sequence[float], center: float) -> float:
    return _median((abs(float(value) - center) for value in values), 0.0)


def _family_centroid(clusters: Sequence[Any]) -> float:
    return _median(cluster.delta_mass_da for cluster in clusters)


def _best_cluster_per_run(clusters: Sequence[Any], centroid: float) -> list[Any]:
    """Keep one representative cluster per RAW to avoid prevalence inflation."""

    by_run: dict[str, list[Any]] = defaultdict(list)
    for cluster in clusters:
        by_run[str(cluster.data_file)].append(cluster)
    selected = []
    for rows in by_run.values():
        selected.append(
            min(
                rows,
                key=lambda cluster: (
                    abs(float(cluster.delta_mass_da) - centroid),
                    -int(cluster.pair_support),
                    -int(cluster.unique_spectrum_support),
                    int(cluster.cluster_rank),
                    str(cluster.cluster_key),
                ),
            )
        )
    return selected


def aggregate_mass_families(
    clusters: Sequence[Any],
    *,
    successful_runs_by_accession: dict[str, set[str]],
    family_tolerance_da: float,
) -> list[MassFamily]:
    """Aggregate putative-PTM clusters into deterministic accession-level families."""

    if family_tolerance_da <= 0.0:
        raise ValueError("family_tolerance_da must be positive")

    by_accession: dict[str, list[Any]] = defaultdict(list)
    for cluster in clusters:
        if str(cluster.classification) == "putative-ptm":
            by_accession[str(cluster.accession)].append(cluster)

    output: list[MassFamily] = []
    for accession, rows in sorted(by_accession.items()):
        groups: list[list[Any]] = []
        centroids: list[float] = []
        for cluster in sorted(
            rows,
            key=lambda item: (
                float(item.delta_mass_da),
                str(item.data_file),
                int(item.cluster_rank),
            ),
        ):
            eligible = [
                (abs(float(cluster.delta_mass_da) - centroid), index)
                for index, centroid in enumerate(centroids)
                if abs(float(cluster.delta_mass_da) - centroid) <= family_tolerance_da
            ]
            if not eligible:
                groups.append([cluster])
                centroids.append(float(cluster.delta_mass_da))
                continue
            _distance, index = min(eligible)
            groups[index].append(cluster)
            centroids[index] = _family_centroid(groups[index])

        successful = len(successful_runs_by_accession.get(accession, set()))
        finalized: list[MassFamily] = []
        for group in groups:
            provisional_center = _family_centroid(group)
            selected = _best_cluster_per_run(group, provisional_center)
            center = _family_centroid(selected)
            # Re-select after the robust center is known in case one run had two
            # nearby reported clusters assigned to the same preliminary family.
            selected = _best_cluster_per_run(group, center)
            center = _family_centroid(selected)
            masses = [float(item.delta_mass_da) for item in selected]
            run_names = tuple(sorted({str(item.data_file) for item in selected}))
            run_count = len(run_names)
            prevalence = run_count / successful if successful else 0.0
            finalized.append(
                MassFamily(
                    accession=accession,
                    family_id="",
                    median_mass_da=center,
                    mass_mad_da=_mad(masses, center),
                    run_count=run_count,
                    successful_runs=successful,
                    run_prevalence=prevalence,
                    median_cluster_rank=_median(item.cluster_rank for item in selected),
                    median_pair_support=_median(item.pair_support for item in selected),
                    median_unique_spectrum_support=_median(
                        item.unique_spectrum_support for item in selected
                    ),
                    median_spectral_similarity=_median(
                        item.median_spectral_similarity for item in selected
                    ),
                    median_support_fraction=_median(
                        item.support_fraction for item in selected
                    ),
                    median_cluster_sigma_da=_median(
                        item.cluster_sigma_da for item in selected
                    ),
                    run_names=run_names,
                    source_cluster_keys=tuple(
                        sorted(str(item.cluster_key) for item in selected)
                    ),
                )
            )

        finalized.sort(
            key=lambda item: (
                item.median_mass_da,
                -item.run_count,
                item.median_cluster_rank,
            )
        )
        for index, family in enumerate(finalized, start=1):
            output.append(
                MassFamily(
                    **{
                        **asdict(family),
                        "family_id": f"{accession}:MF{index:04d}",
                    }
                )
            )
    return output


def _putative_run_masses(clusters: Sequence[Any]) -> dict[str, dict[str, list[float]]]:
    result: dict[str, dict[str, list[float]]] = defaultdict(lambda: defaultdict(list))
    for cluster in clusters:
        if str(cluster.classification) != "putative-ptm":
            continue
        result[str(cluster.accession)][str(cluster.data_file)].append(
            float(cluster.delta_mass_da)
        )
    return {
        accession: {run: sorted(masses) for run, masses in runs.items()}
        for accession, runs in result.items()
    }


def _has_mass_within(sorted_masses: Sequence[float], mass: float, tolerance: float) -> bool:
    left = bisect.bisect_left(sorted_masses, mass - tolerance)
    return left < len(sorted_masses) and sorted_masses[left] <= mass + tolerance


def prevalence_counts_for_family(
    family: MassFamily,
    *,
    run_masses_by_accession: dict[str, dict[str, list[float]]],
    successful_runs_by_accession: dict[str, set[str]],
    family_tolerance_da: float,
) -> tuple[int, int, int, int]:
    accession_hits = family.run_count
    accession_total = family.successful_runs
    background_hits = 0
    background_total = 0
    for accession, successful_runs in sorted(successful_runs_by_accession.items()):
        if accession == family.accession:
            continue
        background_total += len(successful_runs)
        run_masses = run_masses_by_accession.get(accession, {})
        for run in successful_runs:
            masses = run_masses.get(run, [])
            if _has_mass_within(
                masses,
                family.median_mass_da,
                family_tolerance_da,
            ):
                background_hits += 1
    return accession_hits, accession_total, background_hits, background_total


def _log_beta(a: int, b: int) -> float:
    return math.lgamma(a) + math.lgamma(b) - math.lgamma(a + b)


def beta_superiority_probability(
    accession_hits: int,
    accession_total: int,
    background_hits: int,
    background_total: int,
    *,
    prior_alpha: int = DEFAULT_BETA_PRIOR_ALPHA,
    prior_beta: int = DEFAULT_BETA_PRIOR_BETA,
) -> float:
    """Exact P(p_accession > p_background) for integer Beta posteriors.

    The formula is exact for positive integer shape parameters. v3 deliberately
    uses a Beta(1,1) prior so all posterior shapes remain positive integers.
    """

    if not 0 <= accession_hits <= accession_total:
        raise ValueError("invalid accession hit/total counts")
    if not 0 <= background_hits <= background_total:
        raise ValueError("invalid background hit/total counts")
    if prior_alpha < 1 or prior_beta < 1:
        raise ValueError("Beta prior shapes must be positive integers")

    a1 = accession_hits + prior_alpha
    b1 = accession_total - accession_hits + prior_beta
    a2 = background_hits + prior_alpha
    b2 = background_total - background_hits + prior_beta

    logs = []
    beta_a2_b2 = _log_beta(a2, b2)
    for index in range(a1):
        logs.append(
            _log_beta(a2 + index, b1 + b2)
            - math.log(b1 + index)
            - _log_beta(1 + index, b1)
            - beta_a2_b2
        )
    maximum = max(logs)
    probability = math.exp(maximum) * sum(math.exp(value - maximum) for value in logs)
    return min(max(probability, 0.0), 1.0)


def enrichment_for_family(
    family: MassFamily,
    *,
    run_masses_by_accession: dict[str, dict[str, list[float]]],
    successful_runs_by_accession: dict[str, set[str]],
    family_tolerance_da: float,
) -> StudyEnrichment:
    a_hits, a_total, b_hits, b_total = prevalence_counts_for_family(
        family,
        run_masses_by_accession=run_masses_by_accession,
        successful_runs_by_accession=successful_runs_by_accession,
        family_tolerance_da=family_tolerance_da,
    )
    accession_prevalence = a_hits / a_total if a_total else 0.0
    background_prevalence = b_hits / b_total if b_total else 0.0
    difference = accession_prevalence - background_prevalence
    if background_prevalence > 0.0:
        ratio = accession_prevalence / background_prevalence
    elif accession_prevalence > 0.0:
        ratio = math.inf
    else:
        ratio = 1.0
    probability = beta_superiority_probability(
        a_hits,
        a_total,
        b_hits,
        b_total,
    )
    return StudyEnrichment(
        accession_hits=a_hits,
        accession_total=a_total,
        background_hits=b_hits,
        background_total=b_total,
        accession_prevalence=accession_prevalence,
        background_prevalence=background_prevalence,
        prevalence_difference=difference,
        prevalence_ratio=ratio,
        study_enrichment_probability=probability,
        study_enrichment_local_fdr=1.0 - probability,
    )


def posterior_q_values(
    rows: Sequence[dict[str, Any]],
) -> dict[str, float]:
    """Within-accession posterior expected-FDR q-values by family id."""

    by_accession: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        by_accession[str(row["pxd_accession"])].append(row)

    result: dict[str, float] = {}
    for accession_rows in by_accession.values():
        ordered = sorted(
            accession_rows,
            key=lambda row: (
                -float(row["study_enrichment_probability"]),
                -float(row["prevalence_difference"]),
                str(row["family_id"]),
            ),
        )
        prefix_fdr = []
        cumulative = 0.0
        for index, row in enumerate(ordered, start=1):
            cumulative += float(row["study_enrichment_local_fdr"])
            prefix_fdr.append(cumulative / index)
        q_values = [1.0] * len(ordered)
        running = 1.0
        for index in range(len(ordered) - 1, -1, -1):
            running = min(running, prefix_fdr[index])
            q_values[index] = running
        for row, q_value in zip(ordered, q_values, strict=True):
            result[str(row["family_id"])] = q_value
    return result


def _catalog_index(catalog: Sequence[Any]) -> tuple[list[float], list[Any]]:
    biological = sorted(
        (item for item in catalog if item.category == "biological-ptm"),
        key=lambda item: (float(item.mass_da), str(item.accession), str(item.name)),
    )
    return [float(item.mass_da) for item in biological], biological


def annotate_family(
    family: MassFamily,
    catalog: Sequence[Any],
    *,
    annotation_window_da: float,
    systematic_mass_floor_da: float,
) -> FamilyAnnotation | None:
    masses, biological = _catalog_index(catalog)
    left = bisect.bisect_left(masses, family.median_mass_da - annotation_window_da)
    right = bisect.bisect_right(masses, family.median_mass_da + annotation_window_da)
    if left == right:
        return None

    # Deduplicate specificity-expanded catalogue rows with the same UniMod/name/mass.
    unique: dict[tuple[str, str, float], Any] = {}
    for record in biological[left:right]:
        key = (
            str(record.accession),
            str(record.name),
            round(float(record.mass_da), 9),
        )
        unique[key] = record
    records = sorted(
        unique.values(),
        key=lambda record: (
            abs(family.median_mass_da - float(record.mass_da)),
            str(record.accession).casefold(),
            str(record.name).casefold(),
        ),
    )
    if not records:
        return None

    robust_sem = (
        1.4826 * family.mass_mad_da / math.sqrt(max(family.run_count, 1))
        if family.run_count > 0
        else 0.0
    )
    sigma = max(systematic_mass_floor_da, robust_sem, 1e-9)
    residuals = [family.median_mass_da - float(record.mass_da) for record in records]
    log_weights = [-0.5 * (abs(residual) / sigma) ** 2 for residual in residuals]
    maximum = max(log_weights)
    weights = [math.exp(value - maximum) for value in log_weights]
    total = sum(weights)
    dominance = weights[0] / total if total > 0.0 else 1.0 / len(weights)
    runner = abs(residuals[1]) if len(residuals) > 1 else None
    best_residual = abs(residuals[0])
    margin = runner - best_residual if runner is not None else None
    best = records[0]
    return FamilyAnnotation(
        family=family,
        unimod_accession=str(best.accession),
        unimod_name=str(best.name),
        theoretical_mass_da=float(best.mass_da),
        residual_da=residuals[0],
        candidate_count=len(records),
        candidate_mass_dominance=dominance,
        runner_up_residual_da=runner,
        residual_margin_da=margin,
        mass_likelihood_sigma_da=sigma,
    )


def read_study_evidence(path: Path | None) -> dict[tuple[str, str], dict[str, str]]:
    if path is None:
        return {}
    allowed = {"supported", "conflicting", "not-found"}
    grouped: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"pxd_accession", "unimod_accession", "evidence_status"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise SystemExit(
                "study evidence TSV is missing required columns: "
                + ", ".join(sorted(missing))
            )
        for row in reader:
            status = (row.get("evidence_status") or "").strip().casefold()
            if status not in allowed:
                raise SystemExit(
                    "study evidence status must be supported, conflicting, or not-found"
                )
            key = (
                (row.get("pxd_accession") or "").strip().upper(),
                (row.get("unimod_accession") or "").strip(),
            )
            grouped[key].append(row)

    output: dict[tuple[str, str], dict[str, str]] = {}
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
        sources = sorted(
            {
                (row.get("evidence_source") or "").strip()
                for row in rows
                if (row.get("evidence_source") or "").strip()
            }
        )
        notes = sorted(
            {
                (row.get("evidence_note") or "").strip()
                for row in rows
                if (row.get("evidence_note") or "").strip()
            }
        )
        output[key] = {
            "semantic_evidence_status": status,
            "semantic_evidence_sources": " | ".join(sources),
            "semantic_evidence_notes": " | ".join(notes),
        }
    return output


def family_rows(
    families: Sequence[MassFamily],
    *,
    clusters: Sequence[Any],
    successful_runs_by_accession: dict[str, set[str]],
    catalog: Sequence[Any],
    family_tolerance_da: float,
    annotation_window_da: float,
    systematic_mass_floor_da: float,
) -> list[dict[str, Any]]:
    run_masses = _putative_run_masses(clusters)
    provisional: list[dict[str, Any]] = []
    for family in families:
        enrichment = enrichment_for_family(
            family,
            run_masses_by_accession=run_masses,
            successful_runs_by_accession=successful_runs_by_accession,
            family_tolerance_da=family_tolerance_da,
        )
        annotation = annotate_family(
            family,
            catalog,
            annotation_window_da=annotation_window_da,
            systematic_mass_floor_da=systematic_mass_floor_da,
        )
        provisional.append(
            {
                "pxd_accession": family.accession,
                "family_id": family.family_id,
                "median_mass_da": family.median_mass_da,
                "mass_mad_da": family.mass_mad_da,
                "family_runs": family.run_count,
                "successful_runs": family.successful_runs,
                "family_run_prevalence": family.run_prevalence,
                "background_hits": enrichment.background_hits,
                "background_runs": enrichment.background_total,
                "background_prevalence": enrichment.background_prevalence,
                "prevalence_difference": enrichment.prevalence_difference,
                "prevalence_ratio": enrichment.prevalence_ratio,
                "study_enrichment_probability": enrichment.study_enrichment_probability,
                "study_enrichment_local_fdr": enrichment.study_enrichment_local_fdr,
                "study_enrichment_q_value": 1.0,
                "median_cluster_rank": family.median_cluster_rank,
                "median_pair_support": family.median_pair_support,
                "median_unique_spectrum_support": family.median_unique_spectrum_support,
                "median_spectral_similarity": family.median_spectral_similarity,
                "median_support_fraction": family.median_support_fraction,
                "median_cluster_sigma_da": family.median_cluster_sigma_da,
                "unimod_accession": "" if annotation is None else annotation.unimod_accession,
                "unimod_name": "" if annotation is None else annotation.unimod_name,
                "unimod_theoretical_delta_mass_da": (
                    "" if annotation is None else annotation.theoretical_mass_da
                ),
                "unimod_residual_da": "" if annotation is None else annotation.residual_da,
                "biological_candidate_count": (
                    0 if annotation is None else annotation.candidate_count
                ),
                "candidate_mass_dominance": (
                    0.0 if annotation is None else annotation.candidate_mass_dominance
                ),
                "runner_up_residual_da": (
                    ""
                    if annotation is None or annotation.runner_up_residual_da is None
                    else annotation.runner_up_residual_da
                ),
                "residual_margin_da": (
                    ""
                    if annotation is None or annotation.residual_margin_da is None
                    else annotation.residual_margin_da
                ),
                "mass_likelihood_sigma_da": (
                    "" if annotation is None else annotation.mass_likelihood_sigma_da
                ),
                "run_names": " | ".join(family.run_names),
                "source_cluster_keys": " | ".join(family.source_cluster_keys),
            }
        )

    q_values = posterior_q_values(provisional)
    for row in provisional:
        row["study_enrichment_q_value"] = q_values[str(row["family_id"])]
    return sorted(
        provisional,
        key=lambda row: (
            str(row["pxd_accession"]),
            -float(row["study_enrichment_probability"]),
            -float(row["prevalence_difference"]),
            float(row["median_mass_da"]),
        ),
    )


def build_sdrf_review(
    rows: Sequence[dict[str, Any]],
    *,
    semantic_evidence: dict[tuple[str, str], dict[str, str]],
    enrichment_probability_threshold: float,
    enrichment_q_value_threshold: float,
    minimum_prevalence_difference: float,
    minimum_prevalence_ratio: float,
    minimum_accession_prevalence: float,
    minimum_family_runs: int,
    minimum_candidate_dominance: float,
    maximum_candidates_per_accession: int,
) -> list[dict[str, Any]]:
    by_accession: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        if not row["unimod_accession"]:
            continue
        if float(row["study_enrichment_probability"]) < enrichment_probability_threshold:
            continue
        if float(row["study_enrichment_q_value"]) > enrichment_q_value_threshold:
            continue
        if float(row["prevalence_difference"]) < minimum_prevalence_difference:
            continue
        if float(row["prevalence_ratio"]) < minimum_prevalence_ratio:
            continue
        if float(row["family_run_prevalence"]) < minimum_accession_prevalence:
            continue
        if int(row["family_runs"]) < minimum_family_runs:
            continue
        if float(row["candidate_mass_dominance"]) < minimum_candidate_dominance:
            continue

        accession = str(row["pxd_accession"])
        unimod_accession = str(row["unimod_accession"])
        evidence = semantic_evidence.get(
            (accession, unimod_accession),
            {
                "semantic_evidence_status": "not-evaluated",
                "semantic_evidence_sources": "",
                "semantic_evidence_notes": "",
            },
        )
        semantic_status = evidence["semantic_evidence_status"]
        if semantic_status == "supported":
            sdrf_status = "review-required"
        elif semantic_status in {"conflicting", "mixed"}:
            sdrf_status = "hold-conflicting-evidence"
        else:
            sdrf_status = "semantic-evidence-required"

        candidate = dict(row)
        candidate.update(evidence)
        candidate.update(
            {
                "sdrf_column": "comment[modification parameters]",
                "sdrf_value": (
                    f"NT={row['unimod_name']};AC={row['unimod_accession']}"
                ),
                "sdrf_status": sdrf_status,
                "sdrf_warning": (
                    "RAW-derived accession-enriched mass family; independent study/search "
                    "evidence is required before asserting an SDRF modification parameter. "
                    "No fixed/variable status, target residue, localization, or peptide "
                    "identity is inferred."
                ),
            }
        )
        by_accession[accession].append(candidate)

    output: list[dict[str, Any]] = []
    for _accession, accession_rows in sorted(by_accession.items()):
        accession_rows.sort(
            key=lambda row: (
                -float(row["study_enrichment_probability"]),
                float(row["study_enrichment_q_value"]),
                -float(row["prevalence_difference"]),
                -float(row["candidate_mass_dominance"]),
                -int(row["family_runs"]),
                abs(float(row["unimod_residual_da"])),
                str(row["unimod_accession"]),
            )
        )
        output.extend(accession_rows[:maximum_candidates_per_accession])
    return output


def _family_fields() -> list[str]:
    return [
        "pxd_accession",
        "family_id",
        "median_mass_da",
        "mass_mad_da",
        "family_runs",
        "successful_runs",
        "family_run_prevalence",
        "background_hits",
        "background_runs",
        "background_prevalence",
        "prevalence_difference",
        "prevalence_ratio",
        "study_enrichment_probability",
        "study_enrichment_local_fdr",
        "study_enrichment_q_value",
        "median_cluster_rank",
        "median_pair_support",
        "median_unique_spectrum_support",
        "median_spectral_similarity",
        "median_support_fraction",
        "median_cluster_sigma_da",
        "unimod_accession",
        "unimod_name",
        "unimod_theoretical_delta_mass_da",
        "unimod_residual_da",
        "biological_candidate_count",
        "candidate_mass_dominance",
        "runner_up_residual_da",
        "residual_margin_da",
        "mass_likelihood_sigma_da",
        "run_names",
        "source_cluster_keys",
    ]


def _review_fields() -> list[str]:
    return [
        "pxd_accession",
        "family_id",
        "median_mass_da",
        "family_runs",
        "successful_runs",
        "family_run_prevalence",
        "background_prevalence",
        "prevalence_difference",
        "prevalence_ratio",
        "study_enrichment_probability",
        "study_enrichment_q_value",
        "unimod_accession",
        "unimod_name",
        "unimod_residual_da",
        "candidate_mass_dominance",
        "biological_candidate_count",
        "median_pair_support",
        "median_unique_spectrum_support",
        "median_spectral_similarity",
        "semantic_evidence_status",
        "semantic_evidence_sources",
        "semantic_evidence_notes",
        "sdrf_column",
        "sdrf_value",
        "sdrf_status",
        "sdrf_warning",
    ]


def _write_tsv(path: Path, rows: Sequence[dict[str, Any]], fields: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def calibrate(args: argparse.Namespace) -> None:
    clusters, candidates = v1.read_mass_shift_tables(args.result_root.resolve())
    successful_runs = v1.successful_runs_by_accession(args.result_root.resolve())
    if len(successful_runs) < 2:
        raise SystemExit("v3 enrichment calibration requires at least two accessions.")

    families = aggregate_mass_families(
        clusters,
        successful_runs_by_accession=successful_runs,
        family_tolerance_da=args.family_tolerance_da,
    )
    if args.catalog_source == "openms":
        catalog = v1.catalog_from_openms()
    else:
        catalog = v1.catalog_from_observed_candidates(candidates)
    biological_catalog = [item for item in catalog if item.category == "biological-ptm"]
    if not biological_catalog:
        raise SystemExit("No biological-PTM catalogue entries are available.")

    rows = family_rows(
        families,
        clusters=clusters,
        successful_runs_by_accession=successful_runs,
        catalog=biological_catalog,
        family_tolerance_da=args.family_tolerance_da,
        annotation_window_da=args.annotation_window_da,
        systematic_mass_floor_da=args.systematic_mass_floor_da,
    )
    semantic = read_study_evidence(args.study_evidence_tsv)
    review = build_sdrf_review(
        rows,
        semantic_evidence=semantic,
        enrichment_probability_threshold=args.enrichment_probability_threshold,
        enrichment_q_value_threshold=args.enrichment_q_value_threshold,
        minimum_prevalence_difference=args.minimum_prevalence_difference,
        minimum_prevalence_ratio=args.minimum_prevalence_ratio,
        minimum_accession_prevalence=args.minimum_accession_prevalence,
        minimum_family_runs=args.minimum_family_runs,
        minimum_candidate_dominance=args.minimum_candidate_dominance,
        maximum_candidates_per_accession=args.maximum_candidates_per_accession,
    )

    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    families_path = output / "mass-shift-probability-v3.accession-families.tsv"
    review_path = output / "mass-shift-probability-v3.sdrf-review.tsv"
    metadata_path = output / "mass-shift-probability-v3.metadata.json"
    _write_tsv(families_path, rows, _family_fields())
    _write_tsv(review_path, review, _review_fields())

    catalog_digest = v1.catalog_sha256(biological_catalog)
    counts_by_status: dict[str, int] = defaultdict(int)
    for item in review:
        counts_by_status[str(item["sdrf_status"])] += 1
    metadata = {
        "model_version": MODEL_VERSION,
        "probability_semantics": (
            "Exact Beta-posterior probability that an accession-level recurrent mass "
            "family has higher run prevalence in the accession than across all other "
            "calibration accessions; not peptide/site/search-parameter/PTM-identity "
            "probability."
        ),
        "result_root": str(args.result_root.resolve()),
        "catalog_source": args.catalog_source,
        "catalog_sha256": catalog_digest,
        "biological_catalog_records": len(biological_catalog),
        "accessions": sorted(successful_runs),
        "reported_clusters": len(clusters),
        "putative_ptm_clusters": sum(
            str(cluster.classification) == "putative-ptm" for cluster in clusters
        ),
        "accession_mass_families": len(families),
        "family_tolerance_da": args.family_tolerance_da,
        "annotation_window_da": args.annotation_window_da,
        "systematic_mass_floor_da": args.systematic_mass_floor_da,
        "beta_prior": {
            "alpha": DEFAULT_BETA_PRIOR_ALPHA,
            "beta": DEFAULT_BETA_PRIOR_BETA,
        },
        "enrichment_probability_threshold": args.enrichment_probability_threshold,
        "enrichment_q_value_threshold": args.enrichment_q_value_threshold,
        "minimum_prevalence_difference": args.minimum_prevalence_difference,
        "minimum_prevalence_ratio": args.minimum_prevalence_ratio,
        "minimum_accession_prevalence": args.minimum_accession_prevalence,
        "minimum_family_runs": args.minimum_family_runs,
        "minimum_candidate_dominance": args.minimum_candidate_dominance,
        "maximum_candidates_per_accession": args.maximum_candidates_per_accession,
        "study_evidence_tsv": (
            "" if args.study_evidence_tsv is None else str(args.study_evidence_tsv.resolve())
        ),
        "sdrf_review_candidates": len(review),
        "sdrf_status_counts": dict(sorted(counts_by_status.items())),
        "sdrf_policy": (
            "RAW enrichment can nominate a review candidate but cannot assert a searched "
            "modification; independent semantic/search evidence is required before review "
            "acceptance."
        ),
        "v1_v2_findings_addressed": [
            "cluster/candidate reality separated from accession-level study relevance",
            "ubiquitous recurrent chemistry penalized by cross-accession prevalence",
            "minimum prevalence difference and ratio prevent tiny but certain effects",
            "cross-run family centroid used for mass identity and isobaric ambiguity",
            "semantic evidence kept as an independent SDRF confidence channel",
        ],
    }
    metadata_path.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")

    print(f"accession_families={families_path}")
    print(f"sdrf_review={review_path}")
    print(f"metadata={metadata_path}")
    print(f"reported_clusters={len(clusters)}")
    print(f"putative_ptm_clusters={metadata['putative_ptm_clusters']}")
    print(f"accession_mass_families={len(families)}")
    print(f"sdrf_review_candidates={len(review)}")
    for status, count in sorted(counts_by_status.items()):
        print(f"sdrf_status_{status.replace('-', '_')}={count}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Accession-level recurrent mass-family enrichment and mass-dominant UniMod "
            "annotation for conservative prideQC SDRF review."
        )
    )
    parser.add_argument("--result-root", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--catalog-source", choices=("openms", "observed"), default="openms")
    parser.add_argument("--study-evidence-tsv", type=Path)
    parser.add_argument("--family-tolerance-da", type=float, default=DEFAULT_FAMILY_TOLERANCE_DA)
    parser.add_argument("--annotation-window-da", type=float, default=DEFAULT_ANNOTATION_WINDOW_DA)
    parser.add_argument(
        "--systematic-mass-floor-da",
        type=float,
        default=DEFAULT_SYSTEMATIC_MASS_FLOOR_DA,
    )
    parser.add_argument("--enrichment-probability-threshold", type=float, default=0.95)
    parser.add_argument("--enrichment-q-value-threshold", type=float, default=0.05)
    parser.add_argument("--minimum-prevalence-difference", type=float, default=0.20)
    parser.add_argument("--minimum-prevalence-ratio", type=float, default=2.0)
    parser.add_argument("--minimum-accession-prevalence", type=float, default=0.10)
    parser.add_argument("--minimum-family-runs", type=int, default=2)
    parser.add_argument("--minimum-candidate-dominance", type=float, default=0.80)
    parser.add_argument("--maximum-candidates-per-accession", type=int, default=3)
    return parser


def _validate_args(args: argparse.Namespace) -> None:
    for name in (
        "enrichment_probability_threshold",
        "enrichment_q_value_threshold",
        "minimum_accession_prevalence",
        "minimum_candidate_dominance",
    ):
        value = float(getattr(args, name))
        if not 0.0 <= value <= 1.0:
            raise SystemExit(f"--{name.replace('_', '-')} must be in [0, 1].")
    if args.family_tolerance_da <= 0.0:
        raise SystemExit("--family-tolerance-da must be positive.")
    if args.annotation_window_da <= 0.0:
        raise SystemExit("--annotation-window-da must be positive.")
    if args.systematic_mass_floor_da <= 0.0:
        raise SystemExit("--systematic-mass-floor-da must be positive.")
    if args.minimum_prevalence_difference < 0.0:
        raise SystemExit("--minimum-prevalence-difference must be non-negative.")
    if args.minimum_prevalence_ratio < 1.0:
        raise SystemExit("--minimum-prevalence-ratio must be >= 1.")
    if args.minimum_family_runs < 1:
        raise SystemExit("--minimum-family-runs must be >= 1.")
    if args.maximum_candidates_per_accession < 1:
        raise SystemExit("--maximum-candidates-per-accession must be >= 1.")
    if args.study_evidence_tsv is not None and not args.study_evidence_tsv.is_file():
        raise SystemExit(f"study evidence TSV not found: {args.study_evidence_tsv}")


def main() -> None:
    args = build_parser().parse_args()
    _validate_args(args)
    calibrate(args)


if __name__ == "__main__":
    main()
