"""Streaming mass-error evidence with bounded memory.

This module estimates *measurement precision* from repeated precursor
observations and likely repeated MS2 spectra. Centroid fragments use their native
peak lists; profile fragments use ephemeral local peak-center estimates. It deliberately does not claim to recover the historical
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



def _profile_peak_centers(
    spectrum: Spectrum,
    limit: int,
    *,
    min_neighbor_fraction: float = 0.01,
) -> tuple[np.ndarray, np.ndarray]:
    """Estimate strong peak centers from a profile spectrum without mutating it.

    Candidate local maxima are ranked by apex intensity. For each candidate,
    the center is estimated from the apex and its immediate neighbours by a
    three-point quadratic fit to log(intensity), the local form of a Gaussian
    profile. Ill-conditioned, non-concave, or out-of-bracket fits are rejected
    rather than falling back to the sampled apex position.

    The calculation is vectorized because profile MS2 scans can contain many
    samples. Returned peak centers are ephemeral evidence only; the input
    spectrum is never centroided in-place or written back to disk.
    """

    mz = np.asarray(spectrum.mz, dtype=float)
    intensity = np.asarray(spectrum.intensity, dtype=float)
    valid = np.isfinite(mz) & (mz > 0) & np.isfinite(intensity) & (intensity >= 0)
    mz, intensity = mz[valid], intensity[valid]
    if mz.size < 3:
        return np.array([], dtype=float), np.array([], dtype=float)

    if np.any(np.diff(mz) < 0):
        order = np.argsort(mz, kind="stable")
        mz, intensity = mz[order], intensity[order]
    maxima = np.flatnonzero(
        (intensity[1:-1] > intensity[:-2])
        & (intensity[1:-1] >= intensity[2:])
        & (intensity[1:-1] > 0)
    ) + 1
    if not maxima.size:
        return np.array([], dtype=float), np.array([], dtype=float)

    if maxima.size > limit:
        strongest = np.argpartition(intensity[maxima], -limit)[-limit:]
        maxima = maxima[strongest]

    apex = intensity[maxima]
    left_i = intensity[maxima - 1]
    right_i = intensity[maxima + 1]
    usable = (
        (left_i > 0)
        & (right_i > 0)
        & (left_i >= apex * min_neighbor_fraction)
        & (right_i >= apex * min_neighbor_fraction)
    )
    maxima = maxima[usable]
    apex = apex[usable]
    left_i = left_i[usable]
    right_i = right_i[usable]
    if not maxima.size:
        return np.array([], dtype=float), np.array([], dtype=float)

    center_sample = mz[maxima]
    x_left = mz[maxima - 1] - center_sample
    x_right = mz[maxima + 1] - center_sample
    y_left = np.log(left_i) - np.log(apex)
    y_right = np.log(right_i) - np.log(apex)
    determinant = x_left * x_right * (x_left - x_right)
    with np.errstate(divide="ignore", invalid="ignore"):
        a = (y_left * x_right - y_right * x_left) / determinant
        b = (x_left * x_left * y_right - x_right * x_right * y_left) / determinant
        offset = -b / (2.0 * a)
    centers = center_sample + offset
    usable = (
        np.isfinite(a)
        & np.isfinite(b)
        & np.isfinite(centers)
        & (a < 0)
        & (centers >= mz[maxima - 1])
        & (centers <= mz[maxima + 1])
    )
    centers = centers[usable]
    heights = apex[usable]
    if not centers.size:
        return np.array([], dtype=float), np.array([], dtype=float)

    order = np.argsort(centers, kind="stable")
    centers, heights = centers[order], heights[order]
    norm = float(np.linalg.norm(heights))
    if norm > 0:
        heights = heights / norm
    return centers.astype(float, copy=False), heights.astype(float, copy=False)

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


def _robust_error(values: list[float]) -> dict[str, float | int] | None:
    """Robust symmetric error summary for pairwise repeated measurements.

    The central scale is MAD-based. A diagnostic ``robust_inlier_fraction``
    reports the fraction of pairwise deltas within three robust pairwise sigmas
    of the median. This intentionally does not remove observations from the
    stored distribution; it quantifies the broad-match outlier component so
    downstream confidence rules can be validated without a mixture-model
    dependency.
    """

    array = np.asarray(values, dtype=float)
    array = array[np.isfinite(array)]
    if array.size < 2:
        return None
    center = float(np.median(array))
    centered_absolute = np.abs(array - center)
    mad = float(np.median(centered_absolute))
    pair_sigma = 1.4826 * mad
    if pair_sigma == 0.0 and array.size >= 4:
        q25, q75 = np.quantile(array, [0.25, 0.75])
        pair_sigma = float((q75 - q25) / 1.349)
    single_sigma = pair_sigma / math.sqrt(2.0)
    inlier_threshold = 3.0 * pair_sigma
    if pair_sigma > 0:
        inlier_mask = centered_absolute <= inlier_threshold
        inlier_count = int(np.count_nonzero(inlier_mask))
        inlier_fraction = float(inlier_count / array.size)
        inlier_p95 = float(np.quantile(centered_absolute[inlier_mask], 0.95))
    else:
        inlier_count = int(np.count_nonzero(centered_absolute == 0))
        inlier_fraction = float(inlier_count / array.size)
        inlier_p95 = 0.0 if inlier_count else float("nan")
    return {
        "pairwise_median": center,
        "pairwise_sigma": pair_sigma,
        "single_measurement_sigma": single_sigma,
        "pairwise_p95_abs": float(np.quantile(np.abs(array), 0.95)),
        "robust_inlier_threshold_3sigma": inlier_threshold,
        "robust_inlier_count": inlier_count,
        "robust_outlier_count": int(array.size - inlier_count),
        "robust_inlier_fraction": inlier_fraction,
        "robust_inlier_p95_abs_centered": inlier_p95,
    }


@dataclass(slots=True)
class RepeatSpectrumMassErrorCollector:
    """Estimate precursor and fragment precision without altering spectra.

    Precursor precision uses repeated MS2 precursor observations with the same
    positive charge, nearby m/z, and nearby retention time. This path does not
    depend on fragment peak representation and therefore remains usable for
    profile MS2 data.

    Fragment precision uses strong peak centers from likely repeated MS2 spectra.
    Native centroid spectra use their existing peak lists. Native profile spectra
    use an ephemeral local Gaussian-apex estimate derived from the original profile
    samples; the spectrum is never modified or written back. Native ``unknown``
    spectra use OpenMS peak-type evidence to select the centroid/profile path and
    otherwise abstain.
    """

    rt_window_seconds: float = 120.0
    precursor_candidate_ppm: float = 20.0
    fragment_match_da: float = 0.2
    fragment_low_res_match_da: float = 0.5
    top_peaks: int = 50
    min_matched_peaks: int = 8
    min_overlap_fraction: float = 0.25
    max_candidates_per_bin: int = 32
    min_spectrum_pairs: int = 25
    min_precursor_clusters: int = 10
    min_precursor_cluster_size: int = 3
    min_fragment_pairs: int = 200
    min_tolerance_pairs: int = 200
    min_tolerance_clusters: int = 100
    tolerance_sigma_multiplier: float = 6.0
    min_fragment_tolerance_pairs: int = 1000
    min_fragment_tolerance_spectra: int = 50
    fragment_tolerance_sigma_multiplier: float = 6.0
    fragment_high_res_max_sigma_da: float = 0.01
    fragment_high_res_max_sigma_ppm: float = 10.0
    fragment_low_res_min_sigma_da: float = 0.01
    fragment_low_res_min_sigma_ppm: float = 20.0
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
    fragment_low_res_errors_da: list[float] = field(default_factory=list)
    fragment_low_res_errors_ppm: list[float] = field(default_factory=list)
    precursor_eligible_ms2: int = 0
    fragment_eligible_ms2: int = 0
    precursor_paired_spectra: int = 0
    precursor_clusters_used: int = 0
    fragment_paired_spectra: int = 0
    fragment_centroid_spectra: int = 0
    fragment_profile_spectra: int = 0
    fragment_profile_centroids: int = 0
    excluded_fragment_profile_or_unknown: int = 0
    excluded_fragment_unknown: int = 0
    excluded_profile_peak_pick_failure: int = 0
    excluded_missing_precursor: int = 0
    excluded_low_fragment_peaks: int = 0

    def __post_init__(self) -> None:
        if not math.isfinite(self.rt_window_seconds) or self.rt_window_seconds <= 0:
            raise ValueError("rt_window_seconds must be finite and positive")
        if not math.isfinite(self.precursor_candidate_ppm) or self.precursor_candidate_ppm <= 0:
            raise ValueError("precursor_candidate_ppm must be finite and positive")
        if not math.isfinite(self.fragment_match_da) or self.fragment_match_da <= 0:
            raise ValueError("fragment_match_da must be finite and positive")
        if (
            not math.isfinite(self.fragment_low_res_match_da)
            or self.fragment_low_res_match_da <= self.fragment_match_da
        ):
            raise ValueError(
                "fragment_low_res_match_da must be finite and greater than fragment_match_da"
            )
        if self.top_peaks < 8 or self.min_matched_peaks < 3:
            raise ValueError("top_peaks/min_matched_peaks are too small")
        if not 0 < self.min_overlap_fraction <= 1:
            raise ValueError("min_overlap_fraction must be in (0, 1]")
        if self.max_candidates_per_bin < 1:
            raise ValueError("max_candidates_per_bin must be positive")
        if self.min_precursor_clusters < 1 or self.min_precursor_cluster_size < 3:
            raise ValueError("precursor cluster support thresholds are too small")
        if self.min_tolerance_pairs < self.min_spectrum_pairs:
            raise ValueError("min_tolerance_pairs must be >= min_spectrum_pairs")
        if self.min_tolerance_clusters < self.min_precursor_clusters:
            raise ValueError("min_tolerance_clusters must be >= min_precursor_clusters")
        if not math.isfinite(self.tolerance_sigma_multiplier) or self.tolerance_sigma_multiplier <= 0:
            raise ValueError("tolerance_sigma_multiplier must be finite and positive")
        if self.min_fragment_tolerance_pairs < self.min_fragment_pairs:
            raise ValueError("min_fragment_tolerance_pairs must be >= min_fragment_pairs")
        if self.min_fragment_tolerance_spectra < 1:
            raise ValueError("min_fragment_tolerance_spectra must be positive")
        if (
            not math.isfinite(self.fragment_tolerance_sigma_multiplier)
            or self.fragment_tolerance_sigma_multiplier <= 0
        ):
            raise ValueError("fragment_tolerance_sigma_multiplier must be finite and positive")
        for name, value in (
            ("fragment_high_res_max_sigma_da", self.fragment_high_res_max_sigma_da),
            ("fragment_high_res_max_sigma_ppm", self.fragment_high_res_max_sigma_ppm),
            ("fragment_low_res_min_sigma_da", self.fragment_low_res_min_sigma_da),
            ("fragment_low_res_min_sigma_ppm", self.fragment_low_res_min_sigma_ppm),
        ):
            if not math.isfinite(value) or value <= 0:
                raise ValueError(f"{name} must be finite and positive")

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

    def _fragment_peaks(self, spectrum: Spectrum) -> tuple[np.ndarray, np.ndarray]:
        representation = spectrum.representation
        if representation == "unknown":
            representation = spectrum.estimated_representation or "unknown"

        if representation == "centroid":
            self.fragment_centroid_spectra += 1
            return _top_peaks(spectrum, self.top_peaks)

        if representation == "profile":
            self.fragment_profile_spectra += 1
            mz, intensity = _profile_peak_centers(spectrum, self.top_peaks)
            self.fragment_profile_centroids += int(mz.size)
            if mz.size < self.min_matched_peaks:
                self.excluded_profile_peak_pick_failure += 1
            return mz, intensity

        self.excluded_fragment_unknown += 1
        self.excluded_fragment_profile_or_unknown += 1
        return np.array([], dtype=float), np.array([], dtype=float)

    def _consume_fragments(self, spectrum: Spectrum, precursor_mz: float, charge: int) -> None:
        mz, intensity = self._fragment_peaks(spectrum)
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
            _, candidate, deltas, means = best
            self.fragment_errors_da.extend(float(x) for x in deltas)
            self.fragment_errors_ppm.extend(
                float(delta / mean_mz * 1e6)
                for delta, mean_mz in zip(deltas, means, strict=True)
                if mean_mz > 0
            )
            wide_deltas, wide_means = _matched_deltas(
                candidate.mz,
                mz,
                self.fragment_low_res_match_da,
            )
            self.fragment_low_res_errors_da.extend(float(x) for x in wide_deltas)
            self.fragment_low_res_errors_ppm.extend(
                float(delta / mean_mz * 1e6)
                for delta, mean_mz in zip(wide_deltas, wide_means, strict=True)
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
        fragment_da_narrow = _robust_error(self.fragment_errors_da)
        fragment_ppm_narrow = _robust_error(self.fragment_errors_ppm)
        fragment_da_wide = _robust_error(self.fragment_low_res_errors_da)
        fragment_ppm_wide = _robust_error(self.fragment_low_res_errors_ppm)

        fragment_resolution_regime = "unavailable"
        fragment_match_window_da = self.fragment_match_da
        fragment_da = fragment_da_narrow
        fragment_ppm = fragment_ppm_narrow
        fragment_window_censored = False
        if (
            fragment_da_narrow is not None
            and fragment_ppm_narrow is not None
            and fragment_da_narrow["single_measurement_sigma"] > 0
            and fragment_ppm_narrow["single_measurement_sigma"] > 0
        ):
            narrow_sigma_da = float(fragment_da_narrow["single_measurement_sigma"])
            narrow_sigma_ppm = float(fragment_ppm_narrow["single_measurement_sigma"])
            if (
                narrow_sigma_da <= self.fragment_high_res_max_sigma_da
                and narrow_sigma_ppm <= self.fragment_high_res_max_sigma_ppm
            ):
                fragment_resolution_regime = "high-resolution"
            elif (
                narrow_sigma_da >= self.fragment_low_res_min_sigma_da
                and narrow_sigma_ppm >= self.fragment_low_res_min_sigma_ppm
            ):
                fragment_resolution_regime = "low-resolution"
                if fragment_da_wide is not None and fragment_ppm_wide is not None:
                    fragment_da = fragment_da_wide
                    fragment_ppm = fragment_ppm_wide
                    fragment_match_window_da = self.fragment_low_res_match_da
                    fragment_window_censored = (
                        float(fragment_da_wide["robust_inlier_threshold_3sigma"])
                        >= 0.9 * self.fragment_low_res_match_da
                    )
                else:
                    fragment_window_censored = True

        enough_precursor = (
            self.precursor_paired_spectra >= self.min_spectrum_pairs
            and self.precursor_clusters_used >= self.min_precursor_clusters
        )
        enough_fragment = len(self.fragment_errors_da) >= self.min_fragment_pairs
        method = (
            "repeat-observation mass-error estimator v5; precursor: same positive charge, <=120 s, "
            "<=20 ppm repeat clusters with >=3 observations; fragment: top-50 peak centers, "
            "centroid directly or profile local Gaussian-apex estimate; repeated-spectrum pairing "
            f"uses <={self.fragment_match_da:g} Da and low-resolution precision may use "
            f"<={self.fragment_low_res_match_da:g} Da after regime classification"
        )
        detail = (
            "Measurement-precision evidence from repeated precursor observations and, where centroid "
            "fragment peak centers are available, likely repeated spectra. Pairwise robust scale is "
            "converted to an approximate single-measurement sigma by dividing by sqrt(2). Profile "
            "fragment centers are ephemeral local estimates from the original profile samples; input "
            "spectra are never modified or written back. This is not the historical database-search "
            "tolerance and is never written into SDRF automatically."
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
                (
                    len(self.fragment_low_res_errors_da)
                    if (
                        fragment_resolution_regime == "low-resolution"
                        and fragment_da is fragment_da_wide
                    )
                    else len(self.fragment_errors_da)
                ),
                self.fragment_eligible_ms2,
                "Da",
            ),
            (
                "estimated_fragment_mass_error_ppm",
                fragment_ppm,
                enough_fragment,
                (
                    len(self.fragment_low_res_errors_ppm)
                    if (
                        fragment_resolution_regime == "low-resolution"
                        and fragment_ppm is fragment_ppm_wide
                    )
                    else len(self.fragment_errors_ppm)
                ),
                self.fragment_eligible_ms2,
                "ppm",
            ),
        ):
            payload: dict[str, Any] | None = None
            if enough and value is not None and value["single_measurement_sigma"] > 0:
                payload = {"unit": unit, **value}
                if field_name.startswith("estimated_fragment_mass_error_"):
                    payload.update(
                        {
                            "resolution_regime": fragment_resolution_regime,
                            "fragment_match_window_da": fragment_match_window_da,
                            "window_censored": fragment_window_censored,
                        }
                    )
            result.append(Annotation(
                field_name,
                payload,
                EvidenceKind.INFERRED if payload is not None else EvidenceKind.UNAVAILABLE,
                method,
                detail,
                support=support,
                total=total,
            ))
        tolerance_payload: dict[str, Any] | None = None
        enough_tolerance_support = (
            enough_precursor
            and self.precursor_paired_spectra >= self.min_tolerance_pairs
            and self.precursor_clusters_used >= self.min_tolerance_clusters
        )
        if (
            enough_tolerance_support
            and precursor_ppm is not None
            and precursor_ppm["single_measurement_sigma"] > 0
        ):
            sigma = float(precursor_ppm["single_measurement_sigma"])
            suggested = self.tolerance_sigma_multiplier * sigma
            confidence = (
                "high"
                if self.precursor_paired_spectra >= 1000
                and self.precursor_clusters_used >= 250
                else "moderate"
            )
            tolerance_payload = {
                "unit": "ppm",
                "suggested_tolerance": suggested,
                "single_measurement_sigma": sigma,
                "sigma_multiplier": self.tolerance_sigma_multiplier,
                "confidence": confidence,
                "precursor_clusters": self.precursor_clusters_used,
            }
        result.append(Annotation(
            "suggested_precursor_search_tolerance_ppm",
            tolerance_payload,
            EvidenceKind.INFERRED if tolerance_payload is not None else EvidenceKind.UNAVAILABLE,
            (
                "precision-derived precursor search-tolerance heuristic v1; "
                f"{self.tolerance_sigma_multiplier:g} x robust single-measurement sigma; "
                f">={self.min_tolerance_pairs} repeat differences; "
                f">={self.min_tolerance_clusters} precursor clusters"
            ),
            (
                "Suggested starting tolerance for diverse repeated precursor observations. "
                "The repeat-observation method measures random precision but cannot observe a fixed "
                "calibration offset or choose isotope-error handling. The suggestion is therefore "
                "not a recovered historical search setting, is not guaranteed search-optimal, and "
                "is never written into SDRF automatically. Small fixed-target runs abstain via the "
                "precursor-cluster support requirement."
            ),
            support=self.precursor_paired_spectra,
            total=self.precursor_eligible_ms2,
        ))
        fragment_tolerance_ppm: dict[str, Any] | None = None
        fragment_tolerance_da: dict[str, Any] | None = None
        selected_fragment_pairs = (
            len(self.fragment_low_res_errors_da)
            if (
                fragment_resolution_regime == "low-resolution"
                and fragment_da is fragment_da_wide
            )
            else len(self.fragment_errors_da)
        )
        enough_fragment_tolerance_support = (
            enough_fragment
            and selected_fragment_pairs >= self.min_fragment_tolerance_pairs
            and self.fragment_paired_spectra >= self.min_fragment_tolerance_spectra
        )
        if (
            enough_fragment_tolerance_support
            and fragment_da is not None
            and fragment_ppm is not None
            and fragment_da["single_measurement_sigma"] > 0
            and fragment_ppm["single_measurement_sigma"] > 0
        ):
            sigma_da = float(fragment_da["single_measurement_sigma"])
            sigma_ppm = float(fragment_ppm["single_measurement_sigma"])
            confidence = (
                "high"
                if selected_fragment_pairs >= 10_000
                and self.fragment_paired_spectra >= 500
                else "moderate"
            )
            common_payload = {
                "sigma_multiplier": self.fragment_tolerance_sigma_multiplier,
                "confidence": confidence,
                "fragment_pairs": selected_fragment_pairs,
                "paired_spectra": self.fragment_paired_spectra,
                "fragment_match_window_da": fragment_match_window_da,
                "window_censored": fragment_window_censored,
            }
            if fragment_resolution_regime == "high-resolution":
                fragment_tolerance_ppm = {
                    "unit": "ppm",
                    "suggested_tolerance": self.fragment_tolerance_sigma_multiplier * sigma_ppm,
                    "single_measurement_sigma": sigma_ppm,
                    "resolution_regime": "high-resolution",
                    "robust_inlier_fraction": fragment_ppm["robust_inlier_fraction"],
                    "robust_inlier_count": fragment_ppm["robust_inlier_count"],
                    "robust_outlier_count": fragment_ppm["robust_outlier_count"],
                    **common_payload,
                }
            elif fragment_resolution_regime == "low-resolution" and not fragment_window_censored:
                fragment_tolerance_da = {
                    "unit": "Da",
                    "suggested_tolerance": self.fragment_tolerance_sigma_multiplier * sigma_da,
                    "single_measurement_sigma": sigma_da,
                    "resolution_regime": "low-resolution",
                    "robust_inlier_fraction": fragment_da["robust_inlier_fraction"],
                    "robust_inlier_count": fragment_da["robust_inlier_count"],
                    "robust_outlier_count": fragment_da["robust_outlier_count"],
                    **common_payload,
                }

        fragment_tolerance_method = (
            "precision-derived fragment search-tolerance heuristic v3; "
            f"{self.fragment_tolerance_sigma_multiplier:g} x robust single-measurement sigma; "
            f">={self.min_fragment_tolerance_pairs} matched fragment differences; "
            f">={self.min_fragment_tolerance_spectra} paired spectra; unit selected only for "
            "clearly separated high- or low-resolution precision regimes; low-resolution "
            f"precision uses a {self.fragment_low_res_match_da:g} Da window only after repeated-spectrum "
            f"pairing with {self.fragment_match_da:g} Da"
        )
        fragment_tolerance_detail = (
            "Suggested starting fragment tolerance from repeated-spectrum measurement precision. "
            "High-resolution evidence emits ppm only; low-resolution evidence emits Da only; "
            "intermediate or discordant regimes abstain. Low-resolution recommendations also abstain "
            "when the robust 3-sigma pairwise core approaches the wider matching-window boundary, "
            "preventing a window-truncated precision estimate from becoming a tolerance suggestion. "
            "The payload reports the fraction of broad-match fragment deltas within three robust "
            "pairwise sigmas as a diagnostic only. Profile-mode evidence uses ephemeral local "
            "peak-center estimates and never modifies the input spectra. The suggestion cannot "
            "recover historical search settings or guarantee search-optimal parameters and is never "
            "written into SDRF automatically."
        )
        for field_name, payload in (
            ("suggested_fragment_search_tolerance_ppm", fragment_tolerance_ppm),
            ("suggested_fragment_search_tolerance_da", fragment_tolerance_da),
        ):
            result.append(Annotation(
                field_name,
                payload,
                EvidenceKind.INFERRED if payload is not None else EvidenceKind.UNAVAILABLE,
                fragment_tolerance_method,
                fragment_tolerance_detail,
                support=selected_fragment_pairs,
                total=self.fragment_eligible_ms2,
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
                "fragment_low_res_pairs": len(self.fragment_low_res_errors_da),
                "fragment_resolution_regime": fragment_resolution_regime,
                "fragment_match_window_da": fragment_match_window_da,
                "fragment_window_censored": fragment_window_censored,
                "fragment_centroid_spectra": self.fragment_centroid_spectra,
                "fragment_profile_spectra": self.fragment_profile_spectra,
                "fragment_profile_centroids": self.fragment_profile_centroids,
                "excluded_fragment_profile_or_unknown": self.excluded_fragment_profile_or_unknown,
                "excluded_fragment_unknown": self.excluded_fragment_unknown,
                "excluded_profile_peak_pick_failure": self.excluded_profile_peak_pick_failure,
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
