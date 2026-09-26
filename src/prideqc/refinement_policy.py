"""Deterministic policy gates for confidence-gated SDRF reconstruction.

RAW-derived evidence may reconstruct missing canonical metadata only when both the
scientific confidence gate and the canonical parameter-scope gate are satisfied.
Existing reported values remain authoritative unless an exact no-op match is observed.
High-confidence PTM identity is model eligibility, not automatic canonical writeback.
"""

from __future__ import annotations

import math
import re
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from statistics import median
from typing import Any

PRE_ADJUDICATION_POLICY_VERSION = "prideqc-pre-adjudication-policy-v3"

_MISSING_VALUES = {
    "",
    "n/a",
    "na",
    "none",
    "not applicable",
    "not available",
    "null",
    "unknown",
}
_TOLERANCE_RE = re.compile(r"^\s*([0-9]+(?:\.[0-9]+)?)\s*(ppm|da)\s*$", re.IGNORECASE)

# These thresholds deliberately sit above the broad packet/review admission thresholds.
# They are an identity-confidence gate, not a second discovery algorithm.
_PTM_MIN_SUPPORTING_RUNS = 3
_PTM_MIN_RUN_PREVALENCE = 0.95
_PTM_MIN_HIGH_SUPPORT_RUN_FRACTION = 0.90
_PTM_MIN_CANDIDATE_RUN_FRACTION = 0.95
_PTM_MIN_RECURRENCE_PROBABILITY = 0.99
_PTM_MIN_STRONG_RUN_FRACTION = 0.95
_PTM_MAX_RESIDUAL_FRACTION_OF_MATCH_TOLERANCE = 0.25

_TOLERANCE_MIN_SUPPORTING_RUNS = 3
_TOLERANCE_MIN_COVERAGE = 0.80
_TOLERANCE_MIN_HIGH_CONFIDENCE_FRACTION = 0.50
_TOLERANCE_MAX_RELATIVE_MAD = 0.10
_TOLERANCE_MAX_COMMON_MAX_TO_MEDIAN_RATIO = 1.50
_VERIFIED_PARAMETER_SCOPES = frozenset({"accession-complete", "search-config-confirmed"})


def pre_adjudication_policy_metadata() -> dict[str, Any]:
    """Return the frozen thresholds recorded in every adjudication audit."""
    return {
        "version": PRE_ADJUDICATION_POLICY_VERSION,
        "tolerance_reconstruction": {
            "minimum_supporting_runs": _TOLERANCE_MIN_SUPPORTING_RUNS,
            "minimum_coverage": _TOLERANCE_MIN_COVERAGE,
            "minimum_high_confidence_fraction": _TOLERANCE_MIN_HIGH_CONFIDENCE_FRACTION,
            "maximum_relative_mad": _TOLERANCE_MAX_RELATIVE_MAD,
            "maximum_common_max_to_median_ratio": (
                _TOLERANCE_MAX_COMMON_MAX_TO_MEDIAN_RATIO
            ),
            "verified_parameter_scopes": sorted(_VERIFIED_PARAMETER_SCOPES),
        },
        "ptm_high_confidence_raw_identity": {
            "minimum_supporting_runs": _PTM_MIN_SUPPORTING_RUNS,
            "minimum_run_prevalence": _PTM_MIN_RUN_PREVALENCE,
            "minimum_high_support_run_fraction": _PTM_MIN_HIGH_SUPPORT_RUN_FRACTION,
            "minimum_candidate_run_fraction": _PTM_MIN_CANDIDATE_RUN_FRACTION,
            "minimum_recurrent_family_probability": _PTM_MIN_RECURRENCE_PROBABILITY,
            "minimum_strong_identity_run_fraction": _PTM_MIN_STRONG_RUN_FRACTION,
            "maximum_residual_fraction_of_match_tolerance": (
                _PTM_MAX_RESIDUAL_FRACTION_OF_MATCH_TOLERANCE
            ),
            "requires_non_ambiguous_identity": True,
            "requires_biological_ptm_classification": True,
            "verified_parameter_scopes": sorted(_VERIFIED_PARAMETER_SCOPES),
        },
        "missing_original_modification_is_negative_evidence": False,
        "raw_precision_can_fill_missing_reported_tolerance_when_confident": True,
        "raw_ptm_identity_is_model_eligibility_not_automatic_acceptance": True,
        "repository_summary_metadata_is_hard_negative_evidence": False,
    }


@dataclass(frozen=True)
class PolicyResolution:
    """One deterministic pre-LLM decision and its audit metadata."""

    decision_id: str
    decision: str
    selected_value: str | None
    reason: str
    rule: str
    metrics: dict[str, Any]

    def response_item(self) -> dict[str, Any]:
        return {
            "decision_id": self.decision_id,
            "decision": self.decision,
            "selected_value": self.selected_value,
            "reason": self.reason,
        }

    def audit_item(self) -> dict[str, Any]:
        return {
            "decision_id": self.decision_id,
            "decision": self.decision,
            "selected_value": self.selected_value,
            "rule": self.rule,
            "metrics": self.metrics,
        }


def _as_mapping(value: Any) -> Mapping[str, Any]:
    return value if isinstance(value, Mapping) else {}


def _finite_float(value: Any) -> float | None:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return None
    converted = float(value)
    return converted if math.isfinite(converted) else None


def _flatten_original_values(decision: Mapping[str, Any]) -> list[str]:
    original = _as_mapping(decision.get("original"))
    rows = original.get("rows")
    if not isinstance(rows, Sequence) or isinstance(rows, (str, bytes)):
        return []
    values: list[str] = []
    for row in rows:
        if not isinstance(row, Mapping):
            continue
        raw_values = row.get("values")
        if not isinstance(raw_values, Sequence) or isinstance(raw_values, (str, bytes)):
            continue
        values.extend(str(value).strip() for value in raw_values)
    return values


def _parse_tolerance(value: str) -> tuple[float, str] | None:
    match = _TOLERANCE_RE.match(value)
    if match is None:
        return None
    amount = float(match.group(1))
    if not math.isfinite(amount) or amount <= 0:
        return None
    unit = match.group(2).lower()
    return amount, unit


def _same_tolerance(left: tuple[float, str], right: tuple[float, str]) -> bool:
    if left[1] != right[1]:
        return False
    return math.isclose(left[0], right[0], rel_tol=1e-9, abs_tol=1e-12)


def _tolerance_confidence(decision: Mapping[str, Any]) -> tuple[bool, dict[str, Any]]:
    evidence = _as_mapping(decision.get("evidence"))
    raw_runs = evidence.get("per_run_estimates")
    if not isinstance(raw_runs, Sequence) or isinstance(raw_runs, (str, bytes)):
        raw_runs = []
    estimates: list[float] = []
    confidences: list[str] = []
    units: set[str] = set()
    regimes: set[str] = set()
    for item in raw_runs:
        if not isinstance(item, Mapping) or item.get("status") != "available":
            continue
        value = _finite_float(item.get("value"))
        if value is None or value <= 0:
            continue
        estimates.append(value)
        confidence = str(item.get("confidence") or "").strip().casefold()
        if confidence:
            confidences.append(confidence)
        unit = str(item.get("unit") or "").strip().casefold()
        if unit:
            units.add(unit)
        regime = str(item.get("resolution_regime") or "").strip().casefold()
        if regime:
            regimes.add(regime)
    supporting_runs = len(estimates)
    group_runs = _finite_float(evidence.get("group_runs"))
    coverage = (supporting_runs / group_runs) if group_runs and group_runs > 0 else 0.0
    high_fraction = (
        sum(1 for item in confidences if item == "high") / supporting_runs
        if supporting_runs
        else 0.0
    )
    center = median(estimates) if estimates else None
    relative_mad: float | None = None
    if center is not None and center > 0:
        relative_mad = median(abs(value - center) for value in estimates) / center
    common_max = _finite_float(evidence.get("group_common_max"))
    max_to_median: float | None = None
    if common_max is not None and center is not None and center > 0:
        max_to_median = common_max / center
    units_consistent = len(units) <= 1
    regimes_consistent = len(regimes) <= 1
    metrics = {
        "supporting_runs": supporting_runs,
        "group_runs": int(group_runs) if group_runs is not None else None,
        "coverage": coverage,
        "high_confidence_fraction": high_fraction,
        "relative_mad": relative_mad,
        "common_max_to_median_ratio": max_to_median,
        "units_consistent": units_consistent,
        "resolution_regimes_consistent": regimes_consistent,
        "parameter_scope_status": str(decision.get("parameter_scope_status") or "legacy-unverified"),
    }
    passed = (
        supporting_runs >= _TOLERANCE_MIN_SUPPORTING_RUNS
        and coverage >= _TOLERANCE_MIN_COVERAGE
        and high_fraction >= _TOLERANCE_MIN_HIGH_CONFIDENCE_FRACTION
        and relative_mad is not None
        and relative_mad <= _TOLERANCE_MAX_RELATIVE_MAD
        and max_to_median is not None
        and max_to_median <= _TOLERANCE_MAX_COMMON_MAX_TO_MEDIAN_RATIO
        and units_consistent
        and regimes_consistent
        and metrics["parameter_scope_status"] in _VERIFIED_PARAMETER_SCOPES
    )
    return passed, metrics


def _tolerance_resolution(decision: Mapping[str, Any]) -> PolicyResolution | None:
    decision_id = str(decision.get("decision_id") or "")
    candidates = decision.get("candidate_values")
    if not isinstance(candidates, list) or len(candidates) != 1:
        return None
    candidate_value = str(candidates[0])
    candidate = _parse_tolerance(candidate_value)
    if candidate is None:
        return None

    original_values = _flatten_original_values(decision)
    meaningful = [value for value in original_values if value.casefold() not in _MISSING_VALUES]
    parsed = [(value, _parse_tolerance(value)) for value in meaningful]
    valid = [(value, item) for value, item in parsed if item is not None]

    if not meaningful:
        confident, metrics = _tolerance_confidence(decision)
        metrics.update({"original_values": original_values, "candidate_value": candidate_value})
        if confident:
            return PolicyResolution(
                decision_id=decision_id,
                decision="accept",
                selected_value=candidate_value,
                reason=(
                    "The original SDRF tolerance is missing and prideQC's cohort estimate "
                    "passes the frozen stability and parameter-scope confidence gate."
                ),
                rule="tolerance-missing-high-confidence-reconstruction",
                metrics=metrics,
            )
        return PolicyResolution(
            decision_id=decision_id,
            decision="abstain",
            selected_value=None,
            reason=(
                "The original SDRF tolerance is missing, but the prideQC estimate does not "
                "meet the frozen confidence and verified-scope requirements for writeback."
            ),
            rule="tolerance-missing-insufficient-confidence",
            metrics=metrics,
        )
    if len(valid) != len(meaningful):
        return PolicyResolution(
            decision_id=decision_id,
            decision="abstain",
            selected_value=None,
            reason=(
                "The original SDRF contains a non-empty tolerance value that prideQC cannot "
                "parse safely, so the candidate is not applied automatically."
            ),
            rule="tolerance-original-malformed-abstain",
            metrics={"original_values": meaningful, "candidate_value": candidate_value},
        )

    parsed_values = [item for _, item in valid]
    assert all(item is not None for item in parsed_values)
    canonical = [item for item in parsed_values if item is not None]
    if all(_same_tolerance(item, candidate) for item in canonical):
        return PolicyResolution(
            decision_id=decision_id,
            decision="accept",
            selected_value=candidate_value,
            reason="The supplied candidate already matches the reported SDRF tolerance.",
            rule="tolerance-original-already-matches",
            metrics={"original_values": meaningful, "candidate_value": candidate_value},
        )

    unique_originals = {(round(item[0], 12), item[1]) for item in canonical}
    if len(unique_originals) > 1:
        return PolicyResolution(
            decision_id=decision_id,
            decision="abstain",
            selected_value=None,
            reason=(
                "The target experiment group contains conflicting reported tolerance values; "
                "measured precision is insufficient to choose a replacement automatically."
            ),
            rule="tolerance-original-conflict-abstain",
            metrics={"original_values": meaningful, "candidate_value": candidate_value},
        )

    return PolicyResolution(
        decision_id=decision_id,
        decision="reject",
        selected_value=None,
        reason=(
            "A parseable reported SDRF tolerance is already present and differs from the "
            "RAW-derived precision estimate. Measured precision alone does not establish "
            "that the original search tolerance was wrong."
        ),
        rule="tolerance-preserve-reported-value",
        metrics={"original_values": meaningful, "candidate_value": candidate_value},
    )


def _candidate_identity(decision: Mapping[str, Any]) -> tuple[str, str] | None:
    evidence = _as_mapping(decision.get("evidence"))
    options = evidence.get("candidate_options")
    if not isinstance(options, list) or len(options) != 1 or not isinstance(options[0], Mapping):
        return None
    accession = str(options[0].get("accession") or "").strip()
    name = str(options[0].get("name") or "").strip()
    if not accession or not name:
        return None
    return accession, name


def _candidate_is_already_reported(decision: Mapping[str, Any]) -> bool:
    candidates = decision.get("candidate_values")
    if not isinstance(candidates, list) or len(candidates) != 1:
        return False
    candidate = str(candidates[0]).strip().casefold()
    return any(
        value.strip().casefold() == candidate
        for value in _flatten_original_values(decision)
    )


def _observation_supports_identity(
    observation: Mapping[str, Any],
    *,
    accession: str,
    name: str,
) -> tuple[bool, float | None]:
    if str(observation.get("confidence") or "") != "high-support":
        return False, None
    match_tolerance = _finite_float(observation.get("match_tolerance_da"))
    if match_tolerance is None or match_tolerance <= 0:
        return False, None
    options = observation.get("candidate_options")
    if not isinstance(options, list):
        return False, None
    matching: list[Mapping[str, Any]] = []
    for option in options:
        if not isinstance(option, Mapping):
            continue
        if str(option.get("accession") or "").casefold() != accession.casefold():
            continue
        if str(option.get("name") or "").casefold() != name.casefold():
            continue
        matching.append(option)
    if len(matching) != 1:
        return False, None
    option = matching[0]
    if str(option.get("category") or "") != "biological-ptm":
        return False, None
    residual = _finite_float(option.get("residual_da"))
    if residual is None:
        return False, None
    residual_fraction = abs(residual) / match_tolerance
    if residual_fraction > _PTM_MAX_RESIDUAL_FRACTION_OF_MATCH_TOLERANCE:
        return False, residual_fraction
    return True, residual_fraction


def _strong_ptm_run_fraction(
    decision: Mapping[str, Any],
    *,
    accession: str,
    name: str,
) -> tuple[float, list[float]]:
    evidence = _as_mapping(decision.get("evidence"))
    per_run = evidence.get("per_run_support")
    if not isinstance(per_run, list) or not per_run:
        return 0.0, []
    strong_runs = 0
    residual_fractions: list[float] = []
    considered_runs = 0
    for run in per_run:
        if not isinstance(run, Mapping):
            continue
        observations = run.get("observations")
        if not isinstance(observations, list) or not observations:
            continue
        considered_runs += 1
        supported = False
        best_fraction: float | None = None
        for observation in observations:
            if not isinstance(observation, Mapping):
                continue
            matches, fraction = _observation_supports_identity(
                observation,
                accession=accession,
                name=name,
            )
            if fraction is not None and (best_fraction is None or fraction < best_fraction):
                best_fraction = fraction
            if matches:
                supported = True
        if supported:
            strong_runs += 1
            if best_fraction is not None:
                residual_fractions.append(best_fraction)
    if considered_runs == 0:
        return 0.0, []
    return strong_runs / considered_runs, residual_fractions


def _ptm_resolution(decision: Mapping[str, Any]) -> PolicyResolution | None:
    decision_id = str(decision.get("decision_id") or "")
    candidates = decision.get("candidate_values")
    if not isinstance(candidates, list) or len(candidates) != 1:
        return None
    candidate_value = str(candidates[0])

    if _candidate_is_already_reported(decision):
        return PolicyResolution(
            decision_id=decision_id,
            decision="accept",
            selected_value=candidate_value,
            reason="The candidate PTM is already present in the original SDRF metadata.",
            rule="ptm-original-already-reported",
            metrics={"candidate_value": candidate_value},
        )

    evidence = _as_mapping(decision.get("evidence"))
    if evidence.get("mass_identity_ambiguous") is True:
        return PolicyResolution(
            decision_id=decision_id,
            decision="abstain",
            selected_value=None,
            reason=(
                "The recurrent mass family maps to more than one chemical identity, so prideQC "
                "will not promote a PTM into canonical SDRF metadata."
            ),
            rule="ptm-mass-identity-ambiguous",
            metrics={"candidate_value": candidate_value},
        )

    identity = _candidate_identity(decision)
    if identity is None:
        return PolicyResolution(
            decision_id=decision_id,
            decision="abstain",
            selected_value=None,
            reason="The PTM candidate does not have one unambiguous UniMod identity.",
            rule="ptm-identity-not-unique",
            metrics={"candidate_value": candidate_value},
        )
    accession, name = identity

    semantic_status = str(evidence.get("semantic_evidence_status") or "not-evaluated")
    # Repository/project-level semantic metadata is contextual evidence, not a
    # deterministic veto. Direct deposited search evidence is surfaced to the constrained
    # model and can support reject/abstain there.

    supporting_runs = _finite_float(evidence.get("supporting_runs"))
    run_prevalence = _finite_float(evidence.get("run_prevalence"))
    high_support_fraction = _finite_float(evidence.get("high_support_run_fraction"))
    candidate_run_fraction = _finite_float(evidence.get("candidate_run_fraction"))
    recurrence_probability = _finite_float(evidence.get("recurrent_family_probability"))
    strong_run_fraction, residual_fractions = _strong_ptm_run_fraction(
        decision,
        accession=accession,
        name=name,
    )
    metrics = {
        "supporting_runs": supporting_runs,
        "run_prevalence": run_prevalence,
        "high_support_run_fraction": high_support_fraction,
        "candidate_run_fraction": candidate_run_fraction,
        "recurrent_family_probability": recurrence_probability,
        "strong_identity_run_fraction": strong_run_fraction,
        "max_mass_residual_fraction_of_match_tolerance": (
            max(residual_fractions) if residual_fractions else None
        ),
        "semantic_evidence_status": semantic_status,
    }

    high_confidence_raw_identity = (
        supporting_runs is not None
        and supporting_runs >= _PTM_MIN_SUPPORTING_RUNS
        and run_prevalence is not None
        and run_prevalence >= _PTM_MIN_RUN_PREVALENCE
        and high_support_fraction is not None
        and high_support_fraction >= _PTM_MIN_HIGH_SUPPORT_RUN_FRACTION
        and candidate_run_fraction is not None
        and candidate_run_fraction >= _PTM_MIN_CANDIDATE_RUN_FRACTION
        and recurrence_probability is not None
        and recurrence_probability >= _PTM_MIN_RECURRENCE_PROBABILITY
        and strong_run_fraction >= _PTM_MIN_STRONG_RUN_FRACTION
    )
    metrics["raw_identity_gate_met"] = high_confidence_raw_identity
    scope_status = str(decision.get("parameter_scope_status") or "legacy-unverified")
    metrics["parameter_scope_status"] = scope_status
    scope_verified = scope_status in _VERIFIED_PARAMETER_SCOPES

    if not scope_verified:
        return PolicyResolution(
            decision_id=decision_id,
            decision="abstain",
            selected_value=None,
            reason=(
                "The PTM evidence may be strong, but the canonical search-parameter scope "
                "is not verified; prideQC will not write a subset-of-rows search parameter."
            ),
            rule="ptm-parameter-scope-unverified",
            metrics=metrics,
        )

    if high_confidence_raw_identity or semantic_status == "supported":
        return None

    return PolicyResolution(
        decision_id=decision_id,
        decision="abstain",
        selected_value=None,
        reason=(
            "The PTM candidate lacks either the frozen high-confidence RAW identity gate "
            "or independent supporting evidence required for constrained adjudication."
        ),
        rule="ptm-insufficient-identity-confidence",
        metrics=metrics,
    )


def resolve_pre_adjudication_policy(
    decision: Mapping[str, Any],
) -> PolicyResolution | None:
    """Resolve one decision when deterministic scientific policy is sufficient."""
    decision_type = decision.get("decision_type")
    if decision_type == "mass_tolerance":
        return _tolerance_resolution(decision)
    if decision_type == "modification":
        return _ptm_resolution(decision)
    return None
