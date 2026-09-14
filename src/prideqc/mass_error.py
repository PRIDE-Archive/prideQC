"""Streaming repeat-spectrum mass-error evidence with bounded memory.

This module estimates *measurement precision* from likely repeated MS2 spectra.
It deliberately does not claim to recover the historical database-search
settings. Search tolerances remain separate SDRF/search provenance.
"""

from __future__ import annotations

import math
from collections import defaultdict, deque
from dataclasses import dataclass, field
from typing import Any

import numpy as np

from prideqc.models import Annotation, EvidenceKind, Spectrum


@dataclass(frozen=True, slots=True)
class _Fingerprint:
    rt: float
    precursor_mz: float
    charge: int
    mz: np.ndarray
    intensity: np.ndarray


def _top_peaks(spectrum: Spectrum, limit: int) -> tuple[np.ndarray, np.ndarray]:
    mz = np.asarray(spectrum.mz, dtype=float)
    intensity = np.asarray(spectrum.intensity, dtype=float)
    valid = np.isfinite(mz) & (mz > 0) & np.isfinite(intensity) & (intensity > 0)
    mz, intensity = mz[valid], intensity[valid]
    if not mz.size:
        return np.array([], dtype=float), np.array([], dtype=float)
    if mz.size > limit:
        indices = np.argpartition(intensity, -limit)[-limit:]
        mz, intensity = mz[indices], intensity[indices]
    order = np.argsort(mz, kind="stable")
    mz, intensity = mz[order], intensity[order]
    norm = float(np.linalg.norm(intensity))
    if norm > 0:
        intensity = intensity / norm
    return mz, intensity


def _matched_deltas(
    left_mz: np.ndarray,
    right_mz: np.ndarray,
    tolerance_da: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Greedy one-to-one nearest matching for two sorted peak lists."""

    i = j = 0
    deltas: list[float] = []
    means: list[float] = []
    while i < left_mz.size and j < right_mz.size:
        delta = float(right_mz[j] - left_mz[i])
        if abs(delta) <= tolerance_da:
            # If the next peak is closer, advance that side before committing.
            if j + 1 < right_mz.size:
                next_delta = float(right_mz[j + 1] - left_mz[i])
                if abs(next_delta) < abs(delta) and abs(next_delta) <= tolerance_da:
                    j += 1
                    continue
            if i + 1 < left_mz.size:
                next_delta = float(right_mz[j] - left_mz[i + 1])
                if abs(next_delta) < abs(delta) and abs(next_delta) <= tolerance_da:
                    i += 1
                    continue
            deltas.append(delta)
            means.append((float(left_mz[i]) + float(right_mz[j])) / 2.0)
            i += 1
            j += 1
        elif delta < -tolerance_da:
            j += 1
        else:
            i += 1
    return np.asarray(deltas, dtype=float), np.asarray(means, dtype=float)


def _robust_error(values: list[float]) -> dict[str, float] | None:
    """Robust symmetric error summary for pairwise repeated measurements.

    Pairwise differences contain variance from two measurements. The reported
    single-measurement sigma therefore divides the robust pairwise scale by
    sqrt(2). The median is retained to expose residual pairwise centering.
    """

    array = np.asarray(values, dtype=float)
    array = array[np.isfinite(array)]
    if array.size < 2:
        return None
    center = float(np.median(array))
    absolute = np.abs(array - center)
    mad = float(np.median(absolute))
    pair_sigma = 1.4826 * mad
    if pair_sigma == 0.0 and array.size >= 4:
        q25, q75 = np.quantile(array, [0.25, 0.75])
        pair_sigma = float((q75 - q25) / 1.349)
    single_sigma = pair_sigma / math.sqrt(2.0)
    return {
        "pairwise_median": center,
        "single_measurement_sigma": single_sigma,
        "pairwise_p95_abs": float(np.quantile(np.abs(array), 0.95)),
    }


@dataclass(slots=True)
class RepeatSpectrumMassErrorCollector:
    """Estimate precursor/fragment precision from likely repeated MS2 spectra.

    Candidate spectra must have the same positive precursor charge, nearby
    precursor m/z and retention time, and strong overlap among their most
    intense centroid peaks. Only compact top-peak fingerprints are retained.
    """

    rt_window_seconds: float = 120.0
    precursor_candidate_ppm: float = 100.0
    fragment_match_da: float = 0.2
    top_peaks: int = 50
    min_matched_peaks: int = 8
    min_overlap_fraction: float = 0.25
    max_candidates_per_bin: int = 32
    min_spectrum_pairs: int = 25
    min_fragment_pairs: int = 200
    _bins: defaultdict[tuple[int, int], deque[_Fingerprint]] = field(
        default_factory=lambda: defaultdict(deque),
    )
    precursor_errors_ppm: list[float] = field(default_factory=list)
    precursor_errors_da: list[float] = field(default_factory=list)
    fragment_errors_da: list[float] = field(default_factory=list)
    fragment_errors_ppm: list[float] = field(default_factory=list)
    eligible_ms2: int = 0
    paired_spectra: int = 0
    excluded_profile_or_unknown: int = 0
    excluded_missing_precursor: int = 0

    def __post_init__(self) -> None:
        if not math.isfinite(self.rt_window_seconds) or self.rt_window_seconds <= 0:
            raise ValueError("rt_window_seconds must be finite and positive")
        if not math.isfinite(self.precursor_candidate_ppm) or self.precursor_candidate_ppm <= 0:
            raise ValueError("precursor_candidate_ppm must be finite and positive")
        if not math.isfinite(self.fragment_match_da) or self.fragment_match_da <= 0:
            raise ValueError("fragment_match_da must be finite and positive")
        if self.top_peaks < 8 or self.min_matched_peaks < 3:
            raise ValueError("top_peaks/min_matched_peaks are too small")
        if not 0 < self.min_overlap_fraction <= 1:
            raise ValueError("min_overlap_fraction must be in (0, 1]")
        if self.max_candidates_per_bin < 1:
            raise ValueError("max_candidates_per_bin must be positive")

    def consume_spectrum(self, spectrum: Spectrum) -> None:
        if spectrum.ms_level != 2:
            return
        if spectrum.representation != "centroid":
            self.excluded_profile_or_unknown += 1
            return
        first = spectrum.precursors[0] if spectrum.precursors else None
        if (
            first is None
            or first.charge <= 0
            or not math.isfinite(first.mz)
            or first.mz <= 0
        ):
            self.excluded_missing_precursor += 1
            return
        mz, intensity = _top_peaks(spectrum, self.top_peaks)
        if mz.size < self.min_matched_peaks:
            return
        self.eligible_ms2 += 1
        fingerprint = _Fingerprint(spectrum.rt, first.mz, first.charge, mz, intensity)
        # A 1 Th bin plus immediate neighbors keeps candidate lookup bounded,
        # while the actual acceptance test is in ppm below.
        coarse = int(math.floor(first.mz))
        candidates: list[_Fingerprint] = []
        for bin_id in (coarse - 1, coarse, coarse + 1):
            key = (first.charge, bin_id)
            queue = self._bins[key]
            while queue and spectrum.rt - queue[0].rt > self.rt_window_seconds:
                queue.popleft()
            candidates.extend(queue)

        best: tuple[float, _Fingerprint, np.ndarray, np.ndarray] | None = None
        for candidate in candidates:
            denominator = (first.mz + candidate.precursor_mz) / 2.0
            precursor_ppm = (first.mz - candidate.precursor_mz) / denominator * 1e6
            if abs(precursor_ppm) > self.precursor_candidate_ppm:
                continue
            deltas, means = _matched_deltas(candidate.mz, mz, self.fragment_match_da)
            required = max(
                self.min_matched_peaks,
                int(math.ceil(min(candidate.mz.size, mz.size) * self.min_overlap_fraction)),
            )
            if deltas.size < required:
                continue
            score = float(deltas.size / max(candidate.mz.size, mz.size))
            if best is None or score > best[0]:
                best = (score, candidate, deltas, means)

        if best is not None:
            _, candidate, deltas, means = best
            denominator = (first.mz + candidate.precursor_mz) / 2.0
            precursor_da = first.mz - candidate.precursor_mz
            self.precursor_errors_da.append(float(precursor_da))
            self.precursor_errors_ppm.append(float(precursor_da / denominator * 1e6))
            self.fragment_errors_da.extend(float(x) for x in deltas)
            self.fragment_errors_ppm.extend(
                float(delta / mean_mz * 1e6)
                for delta, mean_mz in zip(deltas, means, strict=True)
                if mean_mz > 0
            )
            self.paired_spectra += 1

        key = (first.charge, coarse)
        queue = self._bins[key]
        queue.append(fingerprint)
        while len(queue) > self.max_candidates_per_bin:
            queue.popleft()

    def annotations(self) -> list[Annotation]:
        precursor_ppm = _robust_error(self.precursor_errors_ppm)
        precursor_da = _robust_error(self.precursor_errors_da)
        fragment_da = _robust_error(self.fragment_errors_da)
        fragment_ppm = _robust_error(self.fragment_errors_ppm)
        enough_precursor = self.paired_spectra >= self.min_spectrum_pairs
        enough_fragment = len(self.fragment_errors_da) >= self.min_fragment_pairs
        method = (
            "repeat-spectrum mass-error estimator v1; centroid MS2; same charge; <=120 s; "
            "<=100 ppm precursor candidate window; top-50 peak overlap; <=0.2 Da fragment matching"
        )
        detail = (
            "Measurement-precision evidence from likely repeated spectra. Pairwise robust scale is "
            "converted to an approximate single-measurement sigma by dividing by sqrt(2). This is "
            "not the historical database-search tolerance and is never written into SDRF "
            "automatically."
        )
        common_total = self.eligible_ms2
        result: list[Annotation] = []
        for field_name, value, enough, support, unit in (
            (
                "estimated_precursor_mass_error_ppm",
                precursor_ppm,
                enough_precursor,
                self.paired_spectra,
                "ppm",
            ),
            (
                "estimated_precursor_mass_error_da",
                precursor_da,
                enough_precursor,
                self.paired_spectra,
                "Da",
            ),
            (
                "estimated_fragment_mass_error_da",
                fragment_da,
                enough_fragment,
                len(self.fragment_errors_da),
                "Da",
            ),
            (
                "estimated_fragment_mass_error_ppm",
                fragment_ppm,
                enough_fragment,
                len(self.fragment_errors_ppm),
                "ppm",
            ),
        ):
            payload: dict[str, Any] | None = None
            if enough and value is not None and value["single_measurement_sigma"] > 0:
                payload = {"unit": unit, **value}
            result.append(Annotation(
                field_name,
                payload,
                EvidenceKind.INFERRED if payload is not None else EvidenceKind.UNAVAILABLE,
                method,
                detail,
                support=support,
                total=common_total,
            ))
        result.append(Annotation(
            "mass_error_estimator_diagnostics",
            {
                "eligible_ms2": self.eligible_ms2,
                "paired_spectra": self.paired_spectra,
                "fragment_pairs": len(self.fragment_errors_da),
                "excluded_profile_or_unknown": self.excluded_profile_or_unknown,
                "excluded_missing_precursor": self.excluded_missing_precursor,
            },
            EvidenceKind.INFERRED if self.eligible_ms2 else EvidenceKind.UNAVAILABLE,
            method,
            "Estimator support/accounting only; not an SDRF annotation.",
            support=self.paired_spectra,
            total=common_total,
        ))
        return result
