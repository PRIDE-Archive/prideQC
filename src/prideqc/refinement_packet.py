"""Deterministic accession-level evidence packet for later LLM SDRF adjudication.

The packet is deliberately model-provider agnostic. prideQC computes measurements,
cohort groups and candidate values; a later adjudication layer may only accept,
reject or abstain on the candidates encoded here.
"""

from __future__ import annotations

import hashlib
import math
from collections import defaultdict
from collections.abc import Mapping, Sequence
from math import comb
from pathlib import Path
from statistics import median
from typing import Any

from prideqc.cohort import (
    COHORT_FRAGMENT_FIELD,
    COHORT_PRECURSOR_FIELD,
    CohortSynthesis,
    SemanticEvidence,
)
from prideqc.models import AnalysisResult, Annotation, CVTerm
from prideqc.sdrf import SDRFDocument, file_name

PACKET_SCHEMA_VERSION = "prideqc-llm-refinement-packet-v1"
PACKET_SCOPE = "accession-sdrf"
ALLOWED_DECISIONS = ("accept", "reject", "abstain")
PRECURSOR_COLUMN = "comment[precursor mass tolerance]"
FRAGMENT_COLUMN = "comment[fragment mass tolerance]"
MODIFICATION_COLUMN = "comment[modification parameters]"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _annotation(result: AnalysisResult, field: str) -> Annotation | None:
    # Cohort synthesis appends annotations. Prefer the most recent one if a
    # summary already contained an earlier cohort annotation.
    return next((item for item in reversed(result.annotations) if item.field == field), None)


def _payload(result: AnalysisResult, field: str) -> tuple[Annotation, dict[str, Any]] | None:
    item = _annotation(result, field)
    if item is None or not isinstance(item.value, dict):
        return None
    return item, item.value


def _tolerance_evidence(result: AnalysisResult, *, fragment: bool) -> dict[str, Any]:
    fields = (
        ("suggested_fragment_search_tolerance_ppm", "ppm"),
        ("suggested_fragment_search_tolerance_da", "Da"),
    ) if fragment else (("suggested_precursor_search_tolerance_ppm", "ppm"),)
    for field, fallback_unit in fields:
        found = _payload(result, field)
        if found is None:
            continue
        annotation, payload = found
        value = payload.get("suggested_tolerance")
        if not isinstance(value, (int, float)) or not math.isfinite(float(value)):
            continue
        unit = str(payload.get("unit") or fallback_unit)
        return {
            "status": "available",
            "source_field": field,
            "value": float(value),
            "unit": unit,
            "confidence": payload.get("confidence"),
            "resolution_regime": payload.get("resolution_regime"),
            "support": annotation.support,
            "total": annotation.total,
            "method": annotation.method,
            "detail": annotation.detail,
        }
    return {"status": "unavailable"}


def _candidate_rows(record: dict[str, Any]) -> list[dict[str, Any]]:
    raw = record.get("diagnostic_unimod_candidates") or record.get("unimod_candidates") or []
    if not isinstance(raw, list):
        return []
    output: list[dict[str, Any]] = []
    for candidate in raw:
        if not isinstance(candidate, dict):
            continue
        category = str(candidate.get("candidate_category") or "")
        if category in {"decoy", "amino-acid-substitution"}:
            continue
        accession = str(candidate.get("unimod_accession") or "").strip()
        name = str(candidate.get("name") or "").strip()
        if not accession or not name:
            continue
        output.append(
            {
                "accession": accession,
                "name": name,
                "category": category or None,
                "theoretical_delta_mass_da": candidate.get("theoretical_delta_mass_da"),
                "residual_da": candidate.get("residual_da"),
                "origins": list(candidate.get("origins") or []),
                "term_specificities": list(candidate.get("term_specificities") or []),
            }
        )
    output.sort(key=lambda item: (str(item["accession"]), str(item["name"])))
    return output


def _mass_shift_records(result: AnalysisResult) -> list[dict[str, Any]]:
    item = _annotation(result, "putative_modification_mass_shifts")
    if item is None or not isinstance(item.value, list):
        return []
    return [record for record in item.value if isinstance(record, dict)]


def _matching_ptm_support(
    result: AnalysisResult,
    *,
    accession: str,
    family_mass_da: float,
    tolerance_da: float = 0.02,
) -> list[dict[str, Any]]:
    matches: list[dict[str, Any]] = []
    for record in _mass_shift_records(result):
        value = record.get("delta_mass_da")
        if not isinstance(value, (int, float)) or not math.isfinite(float(value)):
            continue
        if abs(float(value) - family_mass_da) > tolerance_da:
            continue
        candidates = _candidate_rows(record)
        if accession.upper() not in {
            str(candidate["accession"]).upper() for candidate in candidates
        }:
            continue
        matches.append(
            {
                "delta_mass_da": float(value),
                "cluster_sigma_da": record.get("cluster_sigma_da"),
                "pair_support": record.get("pair_support"),
                "unique_spectrum_support": record.get("unique_spectrum_support"),
                "median_spectral_similarity": record.get("median_spectral_similarity"),
                "classification": record.get("classification"),
                "confidence": record.get("confidence"),
                "match_tolerance_da": record.get("match_tolerance_da"),
                "candidate_options": candidates,
            }
        )
    matches.sort(
        key=lambda item: (
            -int(item.get("pair_support") or 0),
            float(item["delta_mass_da"]),
        )
    )
    return matches


def _column_context(
    document: SDRFDocument,
    row_numbers: Sequence[int],
    column: str,
) -> dict[str, Any]:
    indices = document.indices(column)
    rows: list[dict[str, Any]] = []
    for row_number in row_numbers:
        row = document.rows[row_number - 2]
        rows.append(
            {
                "row": row_number,
                "data_file": row[document.file_column],
                "values": [row[index] for index in indices],
            }
        )
    return {
        "column_present": bool(indices),
        "column_count": len(indices),
        "rows": rows,
    }


def _row_mapping(
    document: SDRFDocument,
    results: Sequence[AnalysisResult],
    aliases: dict[str, str] | None,
) -> tuple[dict[int, list[int]], dict[int, list[str]]]:
    lookup = SDRFDocument._result_lookup(results, aliases)
    rows_by_result: dict[int, list[int]] = {id(result): [] for result in results}
    data_files_by_result: dict[int, list[str]] = {id(result): [] for result in results}
    for row_number, row in enumerate(document.rows, start=2):
        data_file = row[document.file_column]
        result = lookup.get(file_name(data_file))
        if result is None:
            continue
        rows_by_result[id(result)].append(row_number)
        data_files_by_result[id(result)].append(data_file)
    return rows_by_result, data_files_by_result


def _group_rows(
    members: Sequence[str],
    by_name: Mapping[str, AnalysisResult],
    rows_by_result: Mapping[int, Sequence[int]],
) -> list[int]:
    rows: set[int] = set()
    for name in members:
        result = by_name.get(name)
        if result is not None:
            rows.update(rows_by_result.get(id(result), ()))
    return sorted(rows)


def _group_tolerance_decision(
    *,
    label: str,
    members: Sequence[str],
    field: str,
    target_column: str,
    decision_suffix: str,
    by_name: Mapping[str, AnalysisResult],
    rows_by_result: Mapping[int, Sequence[int]],
    document: SDRFDocument,
) -> dict[str, Any] | None:
    annotations: list[Annotation] = []
    for name in members:
        result = by_name.get(name)
        if result is None:
            continue
        item = _annotation(result, field)
        if item is not None and item.sdrf_value:
            annotations.append(item)
    if not annotations:
        return None
    first = annotations[0]
    candidate_value = str(first.sdrf_value)
    payload = first.value if isinstance(first.value, dict) else {}
    per_run: list[dict[str, Any]] = []
    is_fragment = field == COHORT_FRAGMENT_FIELD
    for name in sorted(members, key=str.casefold):
        result = by_name.get(name)
        if result is None:
            continue
        evidence = _tolerance_evidence(result, fragment=is_fragment)
        per_run.append({"run_id": name, **evidence})
    target_rows = _group_rows(members, by_name, rows_by_result)
    return {
        "decision_id": f"{label}:{decision_suffix}",
        "decision_type": "mass_tolerance",
        "experiment_group": label,
        "target_field": target_column,
        "write_semantics": "fill_or_replace_canonical_value",
        "target_rows": target_rows,
        "target_runs": sorted(members, key=str.casefold),
        "original": _column_context(document, target_rows, target_column),
        "candidate_values": [candidate_value],
        "allowed_values": [candidate_value],
        "allowed_decisions": list(ALLOWED_DECISIONS),
        "evidence": {
            "unit": payload.get("unit"),
            "group_median": payload.get("median"),
            "group_common_max": payload.get("common_max"),
            "supporting_runs": payload.get("files_with_estimate"),
            "group_runs": payload.get("group_files"),
            "source_measurement_support": payload.get("source_measurement_support"),
            "selection_policy": "maximum-supported-per-run-estimate-rounded-up",
            "minimum_compatible_run_fraction": 0.80,
            "per_run_estimates": per_run,
            "method": first.method,
            "detail": first.detail,
        },
    }


def _beta_prevalence_probability(hits: int, runs: int, threshold: float = 0.10) -> float:
    if runs <= 0 or hits < 0 or hits > runs:
        return 0.0
    total = runs + 1
    probability = sum(
        comb(total, index) * threshold**index * (1.0 - threshold) ** (total - index)
        for index in range(hits + 1)
    )
    return min(1.0, max(0.0, probability))


def _cluster_packet_mass_shifts(
    observations: list[tuple[float, str, dict[str, Any]]],
    tolerance_da: float = 0.02,
) -> list[list[tuple[float, str, dict[str, Any]]]]:
    ordered = sorted(observations, key=lambda item: item[0])
    families: list[list[tuple[float, str, dict[str, Any]]]] = []
    for observation in ordered:
        if not families:
            families.append([observation])
            continue
        center = median(item[0] for item in families[-1])
        if abs(observation[0] - center) <= tolerance_da:
            families[-1].append(observation)
        else:
            families.append([observation])
    return families


def _ptm_observation_sort_key(item: Mapping[str, Any]) -> tuple[int, float]:
    return (
        -int(item.get("pair_support") or 0),
        float(item["delta_mass_da"]),
    )


def _compact_ptm_observation(
    value: float,
    record: Mapping[str, Any],
    allowed_accessions: set[str],
) -> dict[str, Any]:
    candidates = [
        candidate
        for candidate in _candidate_rows(dict(record))
        if str(candidate["accession"]).upper() in allowed_accessions
    ]
    return {
        "delta_mass_da": value,
        "cluster_sigma_da": record.get("cluster_sigma_da"),
        "pair_support": record.get("pair_support"),
        "unique_spectrum_support": record.get("unique_spectrum_support"),
        "median_spectral_similarity": record.get("median_spectral_similarity"),
        "classification": record.get("classification"),
        "confidence": record.get("confidence"),
        "match_tolerance_da": record.get("match_tolerance_da"),
        "candidate_options": candidates,
    }


def _packet_ptm_families(
    results: Sequence[AnalysisResult],
    *,
    label: str,
    project_accession: str | None,
    semantic_evidence: Mapping[tuple[str, str], SemanticEvidence] | None,
) -> list[dict[str, Any]]:
    observations: list[tuple[float, str, dict[str, Any]]] = []
    for result in results:
        for record in _mass_shift_records(result):
            value = record.get("delta_mass_da")
            if isinstance(value, (int, float)) and math.isfinite(float(value)):
                observations.append((float(value), result.input_path.name, record))
    if len(results) < 3:
        return []

    accession_key = (project_accession or "").strip().upper()
    evidence_lookup = semantic_evidence or {}
    output: list[dict[str, Any]] = []
    for family in _cluster_packet_mass_shifts(observations):
        run_hits = {item[1] for item in family}
        if len(run_hits) < 3:
            continue
        prevalence = len(run_hits) / len(results)
        probability = _beta_prevalence_probability(len(run_hits), len(results))
        high_support_runs = {
            run_name
            for _, run_name, record in family
            if str(record.get("confidence") or "") == "high-support"
        }
        high_support_fraction = len(high_support_runs) / len(run_hits)
        if prevalence < 0.90 or high_support_fraction < 0.80 or probability < 0.99:
            continue

        support: dict[str, dict[str, Any]] = {}
        candidate_runs: defaultdict[str, set[str]] = defaultdict(set)
        residuals: defaultdict[str, list[float]] = defaultdict(list)
        for _, run_name, record in family:
            for candidate in _candidate_rows(record):
                accession = str(candidate["accession"])
                candidate_runs[accession].add(run_name)
                residual = candidate.get("residual_da")
                if isinstance(residual, (int, float)) and math.isfinite(float(residual)):
                    residuals[accession].append(float(residual))
                support.setdefault(accession, candidate)

        options: list[dict[str, Any]] = []
        for accession in sorted(support, key=str.casefold):
            run_fraction = len(candidate_runs[accession]) / len(run_hits)
            if run_fraction < 0.80:
                continue
            candidate = support[accession]
            semantic = evidence_lookup.get((accession_key, accession.upper()))
            options.append(
                {
                    **candidate,
                    "supporting_runs": len(candidate_runs[accession]),
                    "family_runs": len(run_hits),
                    "candidate_run_fraction": run_fraction,
                    "median_residual_da": (
                        median(residuals[accession]) if residuals[accession] else None
                    ),
                    "semantic_evidence": {
                        "status": semantic.status if semantic else "not-evaluated",
                        "sources": list(semantic.sources) if semantic else [],
                        "notes": list(semantic.notes) if semantic else [],
                    },
                }
            )
        if not options:
            continue

        allowed_accessions = {str(option["accession"]).upper() for option in options}
        per_run: defaultdict[str, list[dict[str, Any]]] = defaultdict(list)
        for value, run_name, record in family:
            compact = _compact_ptm_observation(value, record, allowed_accessions)
            if compact["candidate_options"]:
                per_run[run_name].append(compact)
        per_run_support: list[dict[str, Any]] = []
        for run_name in sorted(per_run, key=str.casefold):
            run_observations: list[Mapping[str, Any]] = list(per_run[run_name])
            per_run_support.append(
                {
                    "run_id": run_name,
                    "observations": sorted(
                        run_observations,
                        key=_ptm_observation_sort_key,
                    ),
                }
            )
        output.append(
            {
                "experiment_group": label,
                "median_mass_da": median(item[0] for item in family),
                "family_runs": len(run_hits),
                "group_runs": len(results),
                "run_prevalence": prevalence,
                "high_support_run_fraction": high_support_fraction,
                "recurrent_family_probability": probability,
                "median_pair_support": median(
                    int(item[2].get("pair_support", 0) or 0) for item in family
                ),
                "candidate_options": options,
                "mass_identity_ambiguous": len(options) != 1,
                "per_run_support": per_run_support,
            }
        )
    output.sort(key=lambda item: float(item["median_mass_da"]))
    return output


def _ptm_decision(
    *,
    family: Mapping[str, Any],
    members: Sequence[str],
    by_name: Mapping[str, AnalysisResult],
    rows_by_result: Mapping[int, Sequence[int]],
    document: SDRFDocument,
) -> dict[str, Any]:
    label = str(family.get("experiment_group") or "")
    family_mass = float(family.get("median_mass_da") or 0.0)
    options = list(family.get("candidate_options") or [])
    allowed_values = [
        CVTerm(str(option["accession"]), str(option["name"])).sdrf_value()
        for option in options
    ]
    target_rows = _group_rows(members, by_name, rows_by_result)
    return {
        "decision_id": f"{label}:modification-family:{family_mass:.6f}",
        "decision_type": "modification",
        "experiment_group": label,
        "target_field": MODIFICATION_COLUMN,
        "write_semantics": "append_one_allowed_canonical_value_if_accepted_and_missing",
        "target_rows": target_rows,
        "target_runs": sorted(members, key=str.casefold),
        "original": _column_context(document, target_rows, MODIFICATION_COLUMN),
        "candidate_values": allowed_values,
        "allowed_values": allowed_values,
        "allowed_decisions": list(ALLOWED_DECISIONS),
        "evidence": {
            "observed_delta_mass_da": family.get("median_mass_da"),
            "supporting_runs": family.get("family_runs"),
            "group_runs": family.get("group_runs"),
            "run_prevalence": family.get("run_prevalence"),
            "high_support_run_fraction": family.get("high_support_run_fraction"),
            "recurrent_family_probability": family.get("recurrent_family_probability"),
            "median_pair_support": family.get("median_pair_support"),
            "mass_identity_ambiguous": family.get("mass_identity_ambiguous"),
            "candidate_options": options,
            "per_run_support": list(family.get("per_run_support") or []),
        },
        "evidence_semantics": {
            "recurrent_family_probability_is_identity_probability": False,
            "candidate_is_peptide_or_site_localized": False,
            "candidate_was_searched_in_original_analysis": False,
        },
    }


def build_llm_refinement_packet(
    document: SDRFDocument,
    results: Sequence[AnalysisResult],
    synthesis: CohortSynthesis,
    *,
    project_accession: str | None,
    sdrf_path: Path,
    prideqc_version: str,
    aliases: dict[str, str] | None = None,
    source_artifacts: Mapping[str, Mapping[str, str]] | None = None,
    semantic_evidence: Mapping[tuple[str, str], SemanticEvidence] | None = None,
) -> dict[str, Any]:
    """Build one deterministic evidence packet for one accession-level SDRF."""
    ordered_results = sorted(results, key=lambda item: item.input_path.name.casefold())
    rows_by_result, data_files_by_result = _row_mapping(document, ordered_results, aliases)
    by_name = {result.input_path.name: result for result in ordered_results}
    artifact_lookup = {
        name.casefold(): dict(values) for name, values in (source_artifacts or {}).items()
    }

    run_entries: list[dict[str, Any]] = []
    for result in ordered_results:
        run_name = result.input_path.name
        rows = list(rows_by_result.get(id(result), ()))
        run_entries.append(
            {
                "run_id": run_name,
                "experiment_group": synthesis.assignments.get(run_name),
                "sdrf_rows": rows,
                "sdrf_data_files": list(data_files_by_result.get(id(result), ())),
                "original_values": {
                    PRECURSOR_COLUMN: _column_context(document, rows, PRECURSOR_COLUMN),
                    FRAGMENT_COLUMN: _column_context(document, rows, FRAGMENT_COLUMN),
                    MODIFICATION_COLUMN: _column_context(document, rows, MODIFICATION_COLUMN),
                },
                "tolerance_evidence": {
                    "precursor": _tolerance_evidence(result, fragment=False),
                    "fragment": _tolerance_evidence(result, fragment=True),
                },
                "source_artifacts": artifact_lookup.get(run_name.casefold(), {}),
            }
        )

    decisions: list[dict[str, Any]] = []
    group_entries: list[dict[str, Any]] = []
    for label in sorted(synthesis.groups, key=str.casefold):
        summary = synthesis.groups[label]
        members = sorted([str(item) for item in summary.get("members", [])], key=str.casefold)
        decision_ids: list[str] = []
        precursor = _group_tolerance_decision(
            label=label,
            members=members,
            field=COHORT_PRECURSOR_FIELD,
            target_column=PRECURSOR_COLUMN,
            decision_suffix="precursor-mass-tolerance",
            by_name=by_name,
            rows_by_result=rows_by_result,
            document=document,
        )
        if precursor is not None:
            decisions.append(precursor)
            decision_ids.append(str(precursor["decision_id"]))
        fragment = _group_tolerance_decision(
            label=label,
            members=members,
            field=COHORT_FRAGMENT_FIELD,
            target_column=FRAGMENT_COLUMN,
            decision_suffix="fragment-mass-tolerance",
            by_name=by_name,
            rows_by_result=rows_by_result,
            document=document,
        )
        if fragment is not None:
            decisions.append(fragment)
            decision_ids.append(str(fragment["decision_id"]))

        member_results = [by_name[name] for name in members if name in by_name]
        group_ptms = _packet_ptm_families(
            member_results,
            label=label,
            project_accession=project_accession,
            semantic_evidence=semantic_evidence,
        )
        for family in group_ptms:
            ptm = _ptm_decision(
                family=family,
                members=members,
                by_name=by_name,
                rows_by_result=rows_by_result,
                document=document,
            )
            decisions.append(ptm)
            decision_ids.append(str(ptm["decision_id"]))

        group_entries.append(
            {
                "group_id": label,
                "member_runs": members,
                "target_sdrf_rows": _group_rows(members, by_name, rows_by_result),
                "summary": {
                    key: value
                    for key, value in summary.items()
                    if key not in {"ptm_review_families", "sdrf_eligible_ptms", "members"}
                },
                "decision_ids": decision_ids,
            }
        )

    decisions.sort(key=lambda item: str(item["decision_id"]).casefold())
    packet = {
        "schema_version": PACKET_SCHEMA_VERSION,
        "packet_scope": PACKET_SCOPE,
        "project_accession": project_accession,
        "provenance": {
            "prideqc_version": prideqc_version,
            "cohort_grouping_evidence": list(synthesis.evidence),
            "decision_policy": "candidate-generation-only; no LLM adjudication performed",
        },
        "sdrf": {
            "source_name": sdrf_path.name,
            "sha256": _sha256(sdrf_path),
            "row_count": len(document.rows),
            "target_columns": {
                PRECURSOR_COLUMN: {
                    "present": bool(document.indices(PRECURSOR_COLUMN)),
                    "count": len(document.indices(PRECURSOR_COLUMN)),
                },
                FRAGMENT_COLUMN: {
                    "present": bool(document.indices(FRAGMENT_COLUMN)),
                    "count": len(document.indices(FRAGMENT_COLUMN)),
                },
                MODIFICATION_COLUMN: {
                    "present": bool(document.indices(MODIFICATION_COLUMN)),
                    "count": len(document.indices(MODIFICATION_COLUMN)),
                },
            },
        },
        "runs": run_entries,
        "experiment_groups": group_entries,
        "decision_candidates": decisions,
        "llm_contract": {
            "allowed_decisions": list(ALLOWED_DECISIONS),
            "candidate_values_must_come_from_packet": True,
            "llm_must_not_invent_measurements_or_ontology_terms": True,
            "uncertainty_should_abstain": True,
            "packet_does_not_modify_sdrf": True,
        },
        "limitations": {
            "ptm_candidate_scope": (
                "strict recurrent mass families with one or more supplied non-decoy "
                "ontology candidates; ambiguity is retained for LLM/human adjudication"
            ),
            "full_mzqc_documents_embedded": False,
        },
    }
    return packet


def validate_llm_refinement_packet(packet: Mapping[str, Any]) -> None:
    """Validate invariants required by a later constrained adjudication layer."""
    if packet.get("schema_version") != PACKET_SCHEMA_VERSION:
        raise ValueError("Unexpected LLM refinement packet schema version")
    if packet.get("packet_scope") != PACKET_SCOPE:
        raise ValueError("LLM refinement packet must be accession scoped")
    runs = packet.get("runs")
    groups = packet.get("experiment_groups")
    decisions = packet.get("decision_candidates")
    if (
        not isinstance(runs, list)
        or not isinstance(groups, list)
        or not isinstance(decisions, list)
    ):
        raise ValueError("LLM refinement packet runs/groups/decisions must be arrays")
    run_ids = [str(item.get("run_id") or "") for item in runs if isinstance(item, dict)]
    if not run_ids or any(not value for value in run_ids) or len(set(run_ids)) != len(run_ids):
        raise ValueError("LLM refinement packet must contain unique non-empty run IDs")
    decision_ids: set[str] = set()
    for item in decisions:
        if not isinstance(item, dict):
            raise ValueError("LLM refinement decision candidate must be an object")
        decision_id = str(item.get("decision_id") or "")
        if not decision_id or decision_id in decision_ids:
            raise ValueError("LLM refinement decision IDs must be unique and non-empty")
        decision_ids.add(decision_id)
        allowed = item.get("allowed_decisions")
        if allowed != list(ALLOWED_DECISIONS):
            raise ValueError(f"Unexpected allowed decisions for {decision_id}")
        values = item.get("allowed_values")
        candidates = item.get("candidate_values")
        if not isinstance(values, list) or not values:
            raise ValueError(f"Decision {decision_id} must provide allowed candidate values")
        if not isinstance(candidates, list) or candidates != values:
            raise ValueError(
                f"Decision {decision_id} candidate values must equal allowed values"
            )
    for group in groups:
        if not isinstance(group, dict):
            raise ValueError("LLM refinement experiment group must be an object")
        for decision_id in group.get("decision_ids") or []:
            if str(decision_id) not in decision_ids:
                raise ValueError(f"Experiment group references unknown decision {decision_id}")
