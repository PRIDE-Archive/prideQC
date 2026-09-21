#!/usr/bin/env python3
"""Semantic-first PTM evidence fusion for conservative SDRF review.

v4 is deliberately not another RAW-only ranking model. The v1-v3 experiments showed
that recurrent-cluster reality and accession enrichment do not establish which PTM
belongs in a study-level SDRF. v4 therefore requires an independent study-evidence
candidate first, then asks whether frozen prideQC mass-shift output contains a
compatible recurrent mass family.

The two confidence channels remain separate:

* ``raw_prevalence_probability`` is the Beta-posterior probability that the
  compatible mass family occurs in more than a configured fraction of successful
  RAW runs for the accession.
* ``semantic_evidence_status`` comes only from an external evidence TSV
  (publication, PRIDE metadata, deposited search parameters, etc.).

No single posterior is reported for PTM identity, peptide identity, localization, or
original search configuration.
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import math
import statistics
import sys
from collections import defaultdict
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

_V1_PATH = Path(__file__).with_name("calibrate_mass_shift_probability_v1.py")
_V1_SPEC = importlib.util.spec_from_file_location("prideqc_probability_v1", _V1_PATH)
if _V1_SPEC is None or _V1_SPEC.loader is None:  # pragma: no cover
    raise RuntimeError(f"Cannot load v1 probability helpers from {_V1_PATH}")
v1 = importlib.util.module_from_spec(_V1_SPEC)
sys.modules[_V1_SPEC.name] = v1
_V1_SPEC.loader.exec_module(v1)

_V3_PATH = Path(__file__).with_name("calibrate_mass_shift_probability_v3.py")
_V3_SPEC = importlib.util.spec_from_file_location("prideqc_probability_v3", _V3_PATH)
if _V3_SPEC is None or _V3_SPEC.loader is None:  # pragma: no cover
    raise RuntimeError(f"Cannot load v3 probability helpers from {_V3_PATH}")
v3 = importlib.util.module_from_spec(_V3_SPEC)
sys.modules[_V3_SPEC.name] = v3
_V3_SPEC.loader.exec_module(v3)

MODEL_VERSION = "mass-shift-probability-v4.0"
DEFAULT_FAMILY_TOLERANCE_DA = 0.02
DEFAULT_RAW_MATCH_WINDOW_DA = 0.02
DEFAULT_MINIMUM_RAW_PREVALENCE = 0.10
DEFAULT_RAW_PREVALENCE_PROBABILITY = 0.95
DEFAULT_MINIMUM_FAMILY_RUNS = 2

_ALLOWED_CLUSTER_CLASSIFICATIONS = {
    "putative-ptm",
    "sample-prep-modification",
    "putative-modification",
    "mass-compatible-other",
    "unknown",
}


@dataclass(frozen=True, slots=True)
class CatalogRecord:
    accession: str
    name: str
    mass_da: float
    source_classifications: tuple[str, ...]
    origins: tuple[str, ...]
    term_specificities: tuple[str, ...]

    @property
    def has_biological_specificity(self) -> bool:
        for value in self.source_classifications:
            folded = value.casefold()
            if any(
                token in folded
                for token in ("post", "co-translational", "pre-translational", "glycosyl")
            ):
                return True
        return False


@dataclass(frozen=True, slots=True)
class EvidenceRecord:
    accession: str
    unimod_accession: str
    status: str
    sources: tuple[str, ...]
    notes: tuple[str, ...]


def _text(value: Any) -> str:
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    return str(value or "")


def catalog_from_openms_specificity_aware(oms: Any | None = None) -> list[CatalogRecord]:
    """Load UniMod while preserving all per-specificity source classifications.

    The production v22.1 loader intentionally deduplicates OpenMS specificity rows,
    but stores one scalar source classification. That is insufficient for semantic
    evidence fusion because a UniMod accession can occur in multiple contexts.
    v4 accumulates every observed classification before deduplication.
    """

    if oms is None:
        import pyopenms as active_oms
    else:
        active_oms = oms

    db = active_oms.ModificationsDB()
    grouped: dict[tuple[str, str, float], dict[str, Any]] = {}
    for index in range(int(db.getNumberOfModifications())):
        mod = db.getModification(index)
        accession = _text(mod.getUniModAccession()).strip()
        if not accession or not accession.casefold().startswith("unimod:"):
            continue
        mass = abs(float(mod.getDiffMonoMass()))
        if not math.isfinite(mass) or mass <= 0.0:
            continue
        name = _text(mod.getFullName()).strip() or _text(mod.getId()).strip()
        key = (accession, name, round(mass, 9))
        row = grouped.setdefault(
            key,
            {
                "accession": accession,
                "name": name,
                "mass": mass,
                "classifications": set(),
                "origins": set(),
                "terms": set(),
            },
        )
        try:
            classification = mod.getSourceClassification()
            classification_name = _text(
                mod.getSourceClassificationName(classification)
            ).strip()
            if classification_name:
                row["classifications"].add(classification_name)
        except (AttributeError, RuntimeError, TypeError, ValueError):
            pass
        try:
            origin = _text(mod.getOrigin()).strip()
            if origin and origin not in {"X", "."}:
                row["origins"].add(origin)
        except (AttributeError, RuntimeError, TypeError, ValueError):
            pass
        try:
            term = mod.getTermSpecificity()
            term_name = _text(mod.getTermSpecificityName(term)).strip()
            if term_name:
                row["terms"].add(term_name)
        except (AttributeError, RuntimeError, TypeError, ValueError):
            pass

    output = [
        CatalogRecord(
            accession=str(row["accession"]),
            name=str(row["name"]),
            mass_da=float(row["mass"]),
            source_classifications=tuple(sorted(row["classifications"])),
            origins=tuple(sorted(row["origins"])),
            term_specificities=tuple(sorted(row["terms"])),
        )
        for row in grouped.values()
    ]
    return sorted(output, key=lambda item: (item.mass_da, item.accession, item.name))


def catalog_from_v1_records(records: Sequence[Any]) -> list[CatalogRecord]:
    """Fallback/testing adapter; does not recover classifications already discarded."""

    return [
        CatalogRecord(
            accession=str(item.accession),
            name=str(item.name),
            mass_da=float(item.mass_da),
            source_classifications=(str(item.source_classification),)
            if getattr(item, "source_classification", "")
            else (),
            origins=tuple(getattr(item, "origins", ())),
            term_specificities=tuple(getattr(item, "term_specificities", ())),
        )
        for item in records
    ]


def read_semantic_evidence(path: Path) -> dict[tuple[str, str], EvidenceRecord]:
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
            accession = (row.get("pxd_accession") or "").strip().upper()
            unimod = (row.get("unimod_accession") or "").strip()
            status = (row.get("evidence_status") or "").strip().casefold()
            if not accession or not unimod:
                raise SystemExit("study evidence rows require accession and UniMod accession")
            if status not in allowed:
                raise SystemExit(
                    "study evidence status must be supported, conflicting, or not-found"
                )
            grouped[(accession, unimod)].append(row)

    output: dict[tuple[str, str], EvidenceRecord] = {}
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
        output[key] = EvidenceRecord(
            accession=key[0],
            unimod_accession=key[1],
            status=status,
            sources=sources,
            notes=notes,
        )
    return output


def _median(values: Iterable[float], default: float = 0.0) -> float:
    rows = [float(value) for value in values]
    return statistics.median(rows) if rows else default


def _best_cluster_per_run(clusters: Sequence[Any], centroid: float) -> list[Any]:
    by_run: dict[str, list[Any]] = defaultdict(list)
    for cluster in clusters:
        by_run[str(cluster.data_file)].append(cluster)
    return [
        min(
            rows,
            key=lambda item: (
                abs(float(item.delta_mass_da) - centroid),
                -int(item.pair_support),
                -int(item.unique_spectrum_support),
                int(item.cluster_rank),
                str(item.cluster_key),
            ),
        )
        for rows in by_run.values()
    ]


def aggregate_recurrent_mass_families(
    clusters: Sequence[Any],
    *,
    successful_runs_by_accession: dict[str, set[str]],
    family_tolerance_da: float,
) -> list[v3.MassFamily]:
    if family_tolerance_da <= 0.0:
        raise ValueError("family_tolerance_da must be positive")
    by_accession: dict[str, list[Any]] = defaultdict(list)
    for cluster in clusters:
        if str(cluster.classification) in _ALLOWED_CLUSTER_CLASSIFICATIONS:
            by_accession[str(cluster.accession)].append(cluster)

    output: list[v3.MassFamily] = []
    for accession, rows in sorted(by_accession.items()):
        groups: list[list[Any]] = []
        centers: list[float] = []
        for cluster in sorted(
            rows,
            key=lambda item: (
                float(item.delta_mass_da),
                str(item.data_file),
                int(item.cluster_rank),
            ),
        ):
            eligible = [
                (abs(float(cluster.delta_mass_da) - center), index)
                for index, center in enumerate(centers)
                if abs(float(cluster.delta_mass_da) - center) <= family_tolerance_da
            ]
            if not eligible:
                groups.append([cluster])
                centers.append(float(cluster.delta_mass_da))
            else:
                _distance, index = min(eligible)
                groups[index].append(cluster)
                centers[index] = _median(item.delta_mass_da for item in groups[index])

        successful = len(successful_runs_by_accession.get(accession, set()))
        finalized: list[v3.MassFamily] = []
        for group in groups:
            center = _median(item.delta_mass_da for item in group)
            selected = _best_cluster_per_run(group, center)
            center = _median(item.delta_mass_da for item in selected)
            selected = _best_cluster_per_run(group, center)
            center = _median(item.delta_mass_da for item in selected)
            masses = [float(item.delta_mass_da) for item in selected]
            run_names = tuple(sorted({str(item.data_file) for item in selected}))
            run_count = len(run_names)
            finalized.append(
                v3.MassFamily(
                    accession=accession,
                    family_id="",
                    median_mass_da=center,
                    mass_mad_da=_median(abs(mass - center) for mass in masses),
                    run_count=run_count,
                    successful_runs=successful,
                    run_prevalence=run_count / successful if successful else 0.0,
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
        finalized.sort(key=lambda item: item.median_mass_da)
        for index, family in enumerate(finalized, start=1):
            output.append(
                v3.MassFamily(
                    accession=family.accession,
                    family_id=f"{accession}:RF{index:04d}",
                    median_mass_da=family.median_mass_da,
                    mass_mad_da=family.mass_mad_da,
                    run_count=family.run_count,
                    successful_runs=family.successful_runs,
                    run_prevalence=family.run_prevalence,
                    median_cluster_rank=family.median_cluster_rank,
                    median_pair_support=family.median_pair_support,
                    median_unique_spectrum_support=family.median_unique_spectrum_support,
                    median_spectral_similarity=family.median_spectral_similarity,
                    median_support_fraction=family.median_support_fraction,
                    median_cluster_sigma_da=family.median_cluster_sigma_da,
                    run_names=family.run_names,
                    source_cluster_keys=family.source_cluster_keys,
                )
            )
    return output


def beta_tail_probability_above(
    hits: int,
    total: int,
    threshold: float,
) -> float:
    """P(p > threshold | hits,total,Beta(1,1)) for integer posterior shapes."""

    if not 0.0 <= threshold <= 1.0:
        raise ValueError("threshold must be in [0, 1]")
    if hits < 0 or total < hits:
        raise ValueError("invalid hit/total counts")
    if threshold <= 0.0:
        return 1.0
    if threshold >= 1.0:
        return 0.0
    alpha = hits + 1
    beta = total - hits + 1
    n = alpha + beta - 1
    # Beta CDF identity: I_x(a,b)=P(Binomial(n,x)>=a).
    cdf = 0.0
    for k in range(alpha, n + 1):
        cdf += (
            math.comb(n, k)
            * (threshold**k)
            * ((1.0 - threshold) ** (n - k))
        )
    return min(max(1.0 - cdf, 0.0), 1.0)


def _catalog_by_accession(catalog: Sequence[CatalogRecord]) -> dict[str, list[CatalogRecord]]:
    output: dict[str, list[CatalogRecord]] = defaultdict(list)
    for record in catalog:
        output[record.accession].append(record)
    return output


def _find_best_family(
    families: Sequence[v3.MassFamily],
    target_mass: float,
    window_da: float,
) -> v3.MassFamily | None:
    compatible = [
        family
        for family in families
        if abs(float(family.median_mass_da) - target_mass) <= window_da
    ]
    if not compatible:
        return None
    return min(
        compatible,
        key=lambda family: (
            abs(float(family.median_mass_da) - target_mass),
            -int(family.run_count),
            -float(family.median_pair_support),
            family.family_id,
        ),
    )


def _mass_alternatives(
    catalog: Sequence[CatalogRecord],
    mass: float,
    window_da: float,
) -> list[CatalogRecord]:
    return sorted(
        [record for record in catalog if abs(record.mass_da - mass) <= window_da],
        key=lambda item: (abs(item.mass_da - mass), item.accession, item.name),
    )


def build_semantic_review(
    *,
    evidence: dict[tuple[str, str], EvidenceRecord],
    families: Sequence[v3.MassFamily],
    clusters: Sequence[Any],
    successful_runs_by_accession: dict[str, set[str]],
    catalog: Sequence[CatalogRecord],
    raw_match_window_da: float,
    minimum_raw_prevalence: float,
    raw_prevalence_probability_threshold: float,
    minimum_family_runs: int,
    family_tolerance_da: float,
) -> list[dict[str, Any]]:
    families_by_accession: dict[str, list[v3.MassFamily]] = defaultdict(list)
    for family in families:
        families_by_accession[family.accession].append(family)
    catalog_by_accession = _catalog_by_accession(catalog)
    run_masses = v3._putative_run_masses(clusters)

    output: list[dict[str, Any]] = []
    for key, item in sorted(evidence.items()):
        accession, unimod_accession = key
        records = catalog_by_accession.get(unimod_accession, [])
        if not records:
            output.append(
                {
                    "pxd_accession": accession,
                    "unimod_accession": unimod_accession,
                    "unimod_name": "",
                    "semantic_evidence_status": item.status,
                    "semantic_evidence_sources": " | ".join(item.sources),
                    "semantic_evidence_notes": " | ".join(item.notes),
                    "raw_confirmation_status": "catalog-accession-not-found",
                    "sdrf_status": "hold-catalog-unresolved",
                }
            )
            continue

        # UniMod accession should represent one chemical mass; choose the most
        # specificity-rich row if OpenMS exposes aliases with identical mass.
        record = min(
            records,
            key=lambda row: (
                -len(row.source_classifications),
                -len(row.origins),
                row.name.casefold(),
            ),
        )
        family = _find_best_family(
            families_by_accession.get(accession, []),
            record.mass_da,
            raw_match_window_da,
        )
        alternatives = _mass_alternatives(catalog, record.mass_da, raw_match_window_da)
        alt_accessions = tuple(
            sorted({alt.accession for alt in alternatives if alt.accession != record.accession})
        )
        alt_names = tuple(
            sorted({alt.name for alt in alternatives if alt.accession != record.accession})
        )
        raw_status = "not-detected"
        raw_probability = 0.0
        family_runs = 0
        successful_runs = len(successful_runs_by_accession.get(accession, set()))
        prevalence = 0.0
        enrichment_probability: float | str = ""
        enrichment_q: float | str = ""
        background_prevalence: float | str = ""
        family_mass: float | str = ""
        residual: float | str = ""
        median_pair_support: float | str = ""
        median_similarity: float | str = ""
        family_id = ""
        if family is not None:
            family_runs = int(family.run_count)
            successful_runs = int(family.successful_runs)
            prevalence = float(family.run_prevalence)
            raw_probability = beta_tail_probability_above(
                family_runs,
                successful_runs,
                minimum_raw_prevalence,
            )
            raw_status = (
                "confirmed"
                if family_runs >= minimum_family_runs
                and raw_probability >= raw_prevalence_probability_threshold
                else "weak-support"
            )
            enrichment = v3.enrichment_for_family(
                family,
                run_masses_by_accession=run_masses,
                successful_runs_by_accession=successful_runs_by_accession,
                family_tolerance_da=family_tolerance_da,
            )
            enrichment_probability = enrichment.study_enrichment_probability
            background_prevalence = enrichment.background_prevalence
            # q-value intentionally not recomputed across semantic candidates; it
            # is descriptive only in v4 and the independent semantic gate is primary.
            enrichment_q = ""
            family_mass = family.median_mass_da
            residual = family.median_mass_da - record.mass_da
            median_pair_support = family.median_pair_support
            median_similarity = family.median_spectral_similarity
            family_id = family.family_id

        if item.status in {"conflicting", "mixed"}:
            sdrf_status = "hold-conflicting-evidence"
        elif item.status == "not-found":
            sdrf_status = "hold-semantic-not-supported"
        elif raw_status == "confirmed":
            sdrf_status = "review-required"
        else:
            sdrf_status = "semantic-supported-raw-unconfirmed"

        output.append(
            {
                "pxd_accession": accession,
                "unimod_accession": record.accession,
                "unimod_name": record.name,
                "unimod_theoretical_delta_mass_da": record.mass_da,
                "catalog_source_classifications": " | ".join(record.source_classifications),
                "catalog_origins": " | ".join(record.origins),
                "catalog_term_specificities": " | ".join(record.term_specificities),
                "catalog_has_biological_specificity": int(record.has_biological_specificity),
                "mass_compatible_alternative_accessions": " | ".join(alt_accessions),
                "mass_compatible_alternative_names": " | ".join(alt_names),
                "mass_identity_ambiguous": int(bool(alt_accessions)),
                "semantic_evidence_status": item.status,
                "semantic_evidence_sources": " | ".join(item.sources),
                "semantic_evidence_notes": " | ".join(item.notes),
                "raw_confirmation_status": raw_status,
                "raw_family_id": family_id,
                "raw_family_mass_da": family_mass,
                "raw_mass_residual_da": residual,
                "raw_family_runs": family_runs,
                "successful_runs": successful_runs,
                "raw_family_prevalence": prevalence,
                "raw_prevalence_probability": raw_probability,
                "background_prevalence": background_prevalence,
                "study_enrichment_probability": enrichment_probability,
                "study_enrichment_q_value": enrichment_q,
                "median_pair_support": median_pair_support,
                "median_spectral_similarity": median_similarity,
                "sdrf_column": "comment[modification parameters]",
                "sdrf_value": f"NT={record.name};AC={record.accession}",
                "sdrf_status": sdrf_status,
                "sdrf_warning": (
                    "Independent study evidence nominates this UniMod accession; frozen "
                    "RAW evidence only confirms a compatible recurrent mass family. This "
                    "does not establish peptide/site localization or fixed/variable search "
                    "status. Exact/isobaric alternatives remain listed explicitly."
                ),
            }
        )
    return output


def _write_tsv(path: Path, rows: Sequence[dict[str, Any]], fields: Sequence[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(fields))
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def _review_fields() -> list[str]:
    return [
        "pxd_accession",
        "unimod_accession",
        "unimod_name",
        "unimod_theoretical_delta_mass_da",
        "catalog_source_classifications",
        "catalog_origins",
        "catalog_term_specificities",
        "catalog_has_biological_specificity",
        "mass_compatible_alternative_accessions",
        "mass_compatible_alternative_names",
        "mass_identity_ambiguous",
        "semantic_evidence_status",
        "semantic_evidence_sources",
        "semantic_evidence_notes",
        "raw_confirmation_status",
        "raw_family_id",
        "raw_family_mass_da",
        "raw_mass_residual_da",
        "raw_family_runs",
        "successful_runs",
        "raw_family_prevalence",
        "raw_prevalence_probability",
        "background_prevalence",
        "study_enrichment_probability",
        "study_enrichment_q_value",
        "median_pair_support",
        "median_spectral_similarity",
        "sdrf_column",
        "sdrf_value",
        "sdrf_status",
        "sdrf_warning",
    ]


def calibrate(args: argparse.Namespace) -> None:
    clusters, observed_candidates = v1.read_mass_shift_tables(args.result_root.resolve())
    successful_runs = v1.successful_runs_by_accession(args.result_root.resolve())
    evidence = read_semantic_evidence(args.study_evidence_tsv.resolve())
    families = aggregate_recurrent_mass_families(
        clusters,
        successful_runs_by_accession=successful_runs,
        family_tolerance_da=args.family_tolerance_da,
    )
    if args.catalog_source == "openms":
        catalog = catalog_from_openms_specificity_aware()
    else:
        catalog = catalog_from_v1_records(v1.catalog_from_observed_candidates(observed_candidates))
    review = build_semantic_review(
        evidence=evidence,
        families=families,
        clusters=clusters,
        successful_runs_by_accession=successful_runs,
        catalog=catalog,
        raw_match_window_da=args.raw_match_window_da,
        minimum_raw_prevalence=args.minimum_raw_prevalence,
        raw_prevalence_probability_threshold=args.raw_prevalence_probability_threshold,
        minimum_family_runs=args.minimum_family_runs,
        family_tolerance_da=args.family_tolerance_da,
    )

    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    review_path = output / "mass-shift-probability-v4.sdrf-review.tsv"
    metadata_path = output / "mass-shift-probability-v4.metadata.json"
    _write_tsv(review_path, review, _review_fields())

    counts: dict[str, int] = defaultdict(int)
    for row in review:
        counts[str(row["sdrf_status"])] += 1
    metadata = {
        "model_version": MODEL_VERSION,
        "architecture": "semantic-first-independent-evidence-plus-raw-mass-family-confirmation",
        "probability_semantics": (
            "raw_prevalence_probability is the Beta-posterior probability that a "
            "semantic-evidence-nominated mass family occurs in more than the configured "
            "fraction of successful RAW runs. It is not a PTM identity or localization "
            "probability."
        ),
        "result_root": str(args.result_root.resolve()),
        "study_evidence_tsv": str(args.study_evidence_tsv.resolve()),
        "catalog_source": args.catalog_source,
        "catalog_records": len(catalog),
        "reported_clusters": len(clusters),
        "eligible_recurrent_clusters": sum(
            str(item.classification) in _ALLOWED_CLUSTER_CLASSIFICATIONS for item in clusters
        ),
        "recurrent_mass_families": len(families),
        "semantic_evidence_candidates": len(evidence),
        "review_rows": len(review),
        "sdrf_status_counts": dict(sorted(counts.items())),
        "family_tolerance_da": args.family_tolerance_da,
        "raw_match_window_da": args.raw_match_window_da,
        "minimum_raw_prevalence": args.minimum_raw_prevalence,
        "raw_prevalence_probability_threshold": args.raw_prevalence_probability_threshold,
        "minimum_family_runs": args.minimum_family_runs,
        "catalog_policy": (
            "Preserve all OpenMS source classifications and specificities per UniMod "
            "accession/name/mass; semantic evidence may nominate mixed-context entries."
        ),
        "sdrf_policy": (
            "No RAW-only top-N selection. Only independently supported study candidates "
            "can become review-required, and only when a compatible recurrent RAW mass "
            "family is confirmed."
        ),
    }
    metadata_path.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")

    print(f"sdrf_review={review_path}")
    print(f"metadata={metadata_path}")
    print(f"reported_clusters={len(clusters)}")
    print(f"eligible_recurrent_clusters={metadata['eligible_recurrent_clusters']}")
    print(f"recurrent_mass_families={len(families)}")
    print(f"semantic_evidence_candidates={len(evidence)}")
    print(f"review_rows={len(review)}")
    for status, count in sorted(counts.items()):
        print(f"sdrf_status_{status.replace('-', '_')}={count}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Fuse independent study PTM evidence with frozen prideQC recurrent mass "
            "families for conservative SDRF review."
        )
    )
    parser.add_argument("--result-root", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--study-evidence-tsv", required=True, type=Path)
    parser.add_argument("--catalog-source", choices=("openms", "observed"), default="openms")
    parser.add_argument("--family-tolerance-da", type=float, default=DEFAULT_FAMILY_TOLERANCE_DA)
    parser.add_argument("--raw-match-window-da", type=float, default=DEFAULT_RAW_MATCH_WINDOW_DA)
    parser.add_argument(
        "--minimum-raw-prevalence",
        type=float,
        default=DEFAULT_MINIMUM_RAW_PREVALENCE,
    )
    parser.add_argument(
        "--raw-prevalence-probability-threshold",
        type=float,
        default=DEFAULT_RAW_PREVALENCE_PROBABILITY,
    )
    parser.add_argument("--minimum-family-runs", type=int, default=DEFAULT_MINIMUM_FAMILY_RUNS)
    return parser


def _validate_args(args: argparse.Namespace) -> None:
    if not args.study_evidence_tsv.is_file():
        raise SystemExit(f"study evidence TSV not found: {args.study_evidence_tsv}")
    if args.family_tolerance_da <= 0.0:
        raise SystemExit("--family-tolerance-da must be positive")
    if args.raw_match_window_da <= 0.0:
        raise SystemExit("--raw-match-window-da must be positive")
    if not 0.0 <= args.minimum_raw_prevalence <= 1.0:
        raise SystemExit("--minimum-raw-prevalence must be in [0, 1]")
    if not 0.0 <= args.raw_prevalence_probability_threshold <= 1.0:
        raise SystemExit("--raw-prevalence-probability-threshold must be in [0, 1]")
    if args.minimum_family_runs < 1:
        raise SystemExit("--minimum-family-runs must be >= 1")


def main() -> None:
    args = build_parser().parse_args()
    _validate_args(args)
    calibrate(args)


if __name__ == "__main__":
    main()
