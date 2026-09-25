"""Deterministic scientific policy gates before LLM SDRF adjudication.

The gates resolve cases where prideQC can enforce metadata semantics directly and leave
only genuinely semantic/borderline cases for the model. Missing SDRF modification
parameters are not evidence against a PTM: a sufficiently strong, unambiguous RAW-derived
PTM identity may be accepted and later appended by the deterministic SDRF apply stage.
"""

from __future__ import annotations

import math
import re
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Any

PRE_ADJUDICATION_POLICY_VERSION = "prideqc-pre-adjudication-policy-v2"

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


def pre_adjudication_policy_metadata() -> dict[str, Any]:
    """Return the frozen thresholds recorded in every adjudication audit."""
    return {
        "version": PRE_ADJUDICATION_POLICY_VERSION,
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
        },
        "missing_original_modification_is_negative_evidence": False,
        "raw_precision_alone_can_replace_reported_tolerance": False,
        "raw_precision_alone_can_fill_missing_reported_tolerance": False,
        "raw_ptm_identity_alone_can_fill_missing_modification": False,
        "canonical_ptm_write_requires_original_or_semantic_support": True,
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
        return PolicyResolution(
            decision_id=decision_id,
            decision="abstain",
            selected_value=None,
            reason=(
                "The original SDRF does not contain a usable reported search tolerance. "
                "The RAW-derived precision estimate is retained as QC/reanalysis evidence "
                "but cannot establish the historical search setting for canonical SDRF metadata."
            ),
            rule="tolerance-original-missing-abstain",
            metrics={
                "original_values": original_values,
                "candidate_value": candidate_value,
            },
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
    if semantic_status in {"conflicting", "mixed"}:
        return PolicyResolution(
            decision_id=decision_id,
            decision="abstain",
            selected_value=None,
            reason=(
                "Independent study evidence is conflicting for this PTM identity, so prideQC "
                "keeps the candidate out of canonical SDRF metadata."
            ),
            rule="ptm-semantic-evidence-conflicting",
            metrics={"semantic_evidence_status": semantic_status},
        )

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

    if semantic_status == "supported":
        return None

    return PolicyResolution(
        decision_id=decision_id,
        decision="abstain",
        selected_value=None,
        reason=(
            "RAW-derived recurrent mass evidence is retained as QC/reanalysis evidence but "
            "cannot establish that this PTM was searched or reported in the experiment. "
            "Independent semantic/search support is required before a missing PTM may be "
            "promoted into canonical SDRF metadata."
        ),
        rule="ptm-raw-evidence-only-abstain",
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
