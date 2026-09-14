"""Streaming mass-error evidence with bounded memory.

This module estimates *measurement precision* from repeated precursor
observations and, when centroid fragment peaks are available, likely repeated
MS2 spectra. It deliberately does not claim to recover the historical
database-search settings. Search tolerances remain separate SDRF/search
provenance.
"""

from __future__ import annotations

import math
from collections import defaultdict, deque
from dataclasses import dataclass, field
from typing import Any

import numpy as np

from prideqc.models import Annotation, EvidenceKind, Spectrum


@dataclass(slots=True)
class _PrecursorCluster:
    last_rt: float
    charge: int
    mean_mz: float
    count: int = 1
    last_mz: float = 0.0
    pending_da: list[float] = field(default_factory=list)
    pending_ppm: list[float] = field(default_factory=list)
    committed: bool = False


@dataclass(frozen=True, slots=True)
class _FragmentFingerprint:
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
    """Robust symmetric error summary for pairwise repeated measurements."""

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
    """Estimate precursor and fragment precision without altering spectra.

    Precursor precision uses repeated MS2 precursor observations with the same
    positive charge, nearby m/z, and nearby retention time. This path does not
    depend on fragment peak representation and therefore remains usable for
    profile MS2 data.

    Fragment precision is stricter: candidate spectra must additionally have
    centroid fragment peaks and strong overlap among their most intense peaks.
    Explicit native profile spectra are never treated as centroid. Native
    ``unknown`` spectra may use OpenMS peak-type evidence when it independently
    classifies the peak array as centroid.
    """

    rt_window_seconds: float = 120.0
    precursor_candidate_ppm: float = 20.0
    fragment_match_da: float = 0.2
    top_peaks: int = 50
    min_matched_peaks: int = 8
    min_overlap_fraction: float = 0.25
    max_candidates_per_bin: int = 32
    min_spectrum_pairs: int = 25
    min_precursor_clusters: int = 10
    min_precursor_cluster_size: int = 3
    min_fragment_pairs: int = 200
    _precursor_bins: defaultdict[tuple[int, int], list[_PrecursorCluster]] = field(
        default_factory=lambda: defaultdict(list),
    )
    _fragment_bins: defaultdict[tuple[int, int], deque[_FragmentFingerprint]] = field(
        default_factory=lambda: defaultdict(deque),
    )
    precursor_errors_ppm: list[float] = field(default_factory=list)
    precursor_errors_da: list[float] = field(default_factory=list)
    fragment_errors_da: list[float] = field(default_factory=list)
    fragment_errors_ppm: list[float] = field(default_factory=list)
    precursor_eligible_ms2: int = 0
    fragment_eligible_ms2: int = 0
    precursor_paired_spectra: int = 0
    precursor_clusters_used: int = 0
    fragment_paired_spectra: int = 0
    excluded_fragment_profile_or_unknown: int = 0
    excluded_missing_precursor: int = 0
    excluded_low_fragment_peaks: int = 0

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
        if self.min_precursor_clusters < 1 or self.min_precursor_cluster_size < 3:
            raise ValueError("precursor cluster support thresholds are too small")

    def _candidate_precursor_clusters(
        self,
        spectrum: Spectrum,
        precursor_mz: float,
        charge: int,
    ) -> list[_PrecursorCluster]:
        coarse = int(math.floor(precursor_mz))
        candidates: list[_PrecursorCluster] = []
        for bin_id in (coarse - 1, coarse, coarse + 1):
            key = (charge, bin_id)
            active = [
                cluster
                for cluster in self._precursor_bins[key]
                if spectrum.rt - cluster.last_rt <= self.rt_window_seconds
            ]
            self._precursor_bins[key] = active
            candidates.extend(active)
        return candidates

    def _commit_precursor_cluster(self, cluster: _PrecursorCluster) -> None:
        if cluster.committed or cluster.count < self.min_precursor_cluster_size:
            return
        self.precursor_errors_da.extend(cluster.pending_da)
        self.precursor_errors_ppm.extend(cluster.pending_ppm)
        self.precursor_paired_spectra += len(cluster.pending_da)
        self.precursor_clusters_used += 1
        cluster.pending_da.clear()
        cluster.pending_ppm.clear()
        cluster.committed = True

    def _consume_precursor(self, spectrum: Spectrum, precursor_mz: float, charge: int) -> None:
        self.precursor_eligible_ms2 += 1
        best: tuple[float, _PrecursorCluster] | None = None
        for cluster in self._candidate_precursor_clusters(spectrum, precursor_mz, charge):
            denominator = (precursor_mz + cluster.mean_mz) / 2.0
            ppm = (precursor_mz - cluster.mean_mz) / denominator * 1e6
            if abs(ppm) > self.precursor_candidate_ppm:
                continue
            score = abs(ppm)
            if best is None or score < best[0]:
                best = (score, cluster)

        if best is None:
            coarse = int(math.floor(precursor_mz))
            clusters = self._precursor_bins[(charge, coarse)]
            clusters.append(_PrecursorCluster(
                spectrum.rt,
                charge,
                precursor_mz,
                last_mz=precursor_mz,
            ))
            if len(clusters) > self.max_candidates_per_bin:
                del clusters[: len(clusters) - self.max_candidates_per_bin]
            return

        _, cluster = best
        denominator = (precursor_mz + cluster.last_mz) / 2.0
        delta_da = precursor_mz - cluster.last_mz
        delta_ppm = delta_da / denominator * 1e6
        if cluster.committed:
            self.precursor_errors_da.append(float(delta_da))
            self.precursor_errors_ppm.append(float(delta_ppm))
            self.precursor_paired_spectra += 1
        else:
            cluster.pending_da.append(float(delta_da))
            cluster.pending_ppm.append(float(delta_ppm))

        cluster.count += 1
        cluster.mean_mz += (precursor_mz - cluster.mean_mz) / cluster.count
        cluster.last_mz = precursor_mz
        cluster.last_rt = spectrum.rt
        self._commit_precursor_cluster(cluster)

    def _fragment_is_centroid(self, spectrum: Spectrum) -> bool:
        if spectrum.representation == "centroid":
            return True
        return (
            spectrum.representation == "unknown"
            and spectrum.estimated_representation == "centroid"
        )

    def _consume_fragments(self, spectrum: Spectrum, precursor_mz: float, charge: int) -> None:
        if not self._fragment_is_centroid(spectrum):
            self.excluded_fragment_profile_or_unknown += 1
            return
        mz, intensity = _top_peaks(spectrum, self.top_peaks)
        if mz.size < self.min_matched_peaks:
            self.excluded_low_fragment_peaks += 1
            return
        self.fragment_eligible_ms2 += 1
        fingerprint = _FragmentFingerprint(spectrum.rt, precursor_mz, charge, mz, intensity)
        coarse = int(math.floor(precursor_mz))
        candidates: list[_FragmentFingerprint] = []
        for bin_id in (coarse - 1, coarse, coarse + 1):
            queue = self._fragment_bins[(charge, bin_id)]
            while queue and spectrum.rt - queue[0].rt > self.rt_window_seconds:
                queue.popleft()
            candidates.extend(queue)

        best: tuple[float, _FragmentFingerprint, np.ndarray, np.ndarray] | None = None
        for candidate in candidates:
            denominator = (precursor_mz + candidate.precursor_mz) / 2.0
            precursor_ppm = (precursor_mz - candidate.precursor_mz) / denominator * 1e6
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
            _, _, deltas, means = best
            self.fragment_errors_da.extend(float(x) for x in deltas)
            self.fragment_errors_ppm.extend(
                float(delta / mean_mz * 1e6)
                for delta, mean_mz in zip(deltas, means, strict=True)
                if mean_mz > 0
            )
            self.fragment_paired_spectra += 1

        queue = self._fragment_bins[(charge, coarse)]
        queue.append(fingerprint)
        while len(queue) > self.max_candidates_per_bin:
            queue.popleft()

    def consume_spectrum(self, spectrum: Spectrum) -> None:
        if spectrum.ms_level != 2:
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
        self._consume_precursor(spectrum, first.mz, first.charge)
        self._consume_fragments(spectrum, first.mz, first.charge)

    def annotations(self) -> list[Annotation]:
        precursor_ppm = _robust_error(self.precursor_errors_ppm)
        precursor_da = _robust_error(self.precursor_errors_da)
        fragment_da = _robust_error(self.fragment_errors_da)
        fragment_ppm = _robust_error(self.fragment_errors_ppm)
        enough_precursor = (
            self.precursor_paired_spectra >= self.min_spectrum_pairs
            and self.precursor_clusters_used >= self.min_precursor_clusters
        )
        enough_fragment = len(self.fragment_errors_da) >= self.min_fragment_pairs
        method = (
            "repeat-observation mass-error estimator v2; precursor: same positive charge, <=120 s, "
            "<=20 ppm repeat clusters with >=3 observations; fragment: centroid MS2, top-50 peak "
            "overlap, <=0.2 Da matching"
        )
        detail = (
            "Measurement-precision evidence from repeated precursor observations and, where centroid "
            "fragment peaks are available, likely repeated spectra. Pairwise robust scale is converted "
            "to an approximate single-measurement sigma by dividing by sqrt(2). Explicit profile "
            "fragment spectra are never centroided or peak-picked. This is not the historical "
            "database-search tolerance and is never written into SDRF automatically."
        )
        result: list[Annotation] = []
        for field_name, value, enough, support, total, unit in (
            (
                "estimated_precursor_mass_error_ppm",
                precursor_ppm,
                enough_precursor,
                self.precursor_paired_spectra,
                self.precursor_eligible_ms2,
                "ppm",
            ),
            (
                "estimated_precursor_mass_error_da",
                precursor_da,
                enough_precursor,
                self.precursor_paired_spectra,
                self.precursor_eligible_ms2,
                "Da",
            ),
            (
                "estimated_fragment_mass_error_da",
                fragment_da,
                enough_fragment,
                len(self.fragment_errors_da),
                self.fragment_eligible_ms2,
                "Da",
            ),
            (
                "estimated_fragment_mass_error_ppm",
                fragment_ppm,
                enough_fragment,
                len(self.fragment_errors_ppm),
                self.fragment_eligible_ms2,
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
                total=total,
            ))
        result.append(Annotation(
            "mass_error_estimator_diagnostics",
            {
                "precursor_eligible_ms2": self.precursor_eligible_ms2,
                "precursor_paired_spectra": self.precursor_paired_spectra,
                "precursor_clusters_used": self.precursor_clusters_used,
                "fragment_eligible_ms2": self.fragment_eligible_ms2,
                "fragment_paired_spectra": self.fragment_paired_spectra,
                "fragment_pairs": len(self.fragment_errors_da),
                "excluded_fragment_profile_or_unknown": self.excluded_fragment_profile_or_unknown,
                "excluded_low_fragment_peaks": self.excluded_low_fragment_peaks,
                "excluded_missing_precursor": self.excluded_missing_precursor,
            },
            EvidenceKind.INFERRED if self.precursor_eligible_ms2 else EvidenceKind.UNAVAILABLE,
            method,
            "Estimator support/accounting only; not an SDRF annotation.",
            support=self.precursor_paired_spectra,
            total=self.precursor_eligible_ms2,
        ))
        return result
