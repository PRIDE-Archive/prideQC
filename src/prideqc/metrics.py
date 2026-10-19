"""ID-free metrics from compact per-scan summaries, inspired by rawQC.

Peak arrays are reduced once. Exact quantiles retain O(number of spectra)
scalar values, not O(total peaks) data. RT calculations use stable time order.
"""

from __future__ import annotations

import math
import struct
from array import array
from collections import Counter
from dataclasses import dataclass, field
from typing import Any

import numpy as np

from prideqc.models import FloatArray, Metric, Spectrum


def _doubles() -> array:
    return array("d")


def _finite(values: Any) -> FloatArray:
    result = np.asarray(values, dtype=float)
    return result[np.isfinite(result)]


def _range(values: FloatArray) -> list[float] | None:
    return [float(values.min()), float(values.max())] if values.size else None


def _quantiles(values: FloatArray) -> list[float] | None:
    return np.quantile(values, [0.25, 0.5, 0.75]).tolist() if values.size else None


def _mean(values: FloatArray) -> float | None:
    return float(values.mean()) if values.size else None


def _fastest_frequency(rt: FloatArray) -> float | None:
    if not rt.size:
        return None
    upper = np.searchsorted(rt, rt + 60, side="right")
    lower = np.searchsorted(rt, rt, side="left")
    return float(np.max(upper - lower) / 60)


def _interval_fractions(rt: FloatArray, span: float) -> list[float] | None:
    if span <= 0:
        return None
    return (np.diff(np.quantile(rt, [0, 0.25, 0.5, 0.75, 1])) / span).tolist()


def _mode_and_fraction(values: list[int]) -> tuple[int | None, float | None]:
    if not values:
        return None, None
    counts = Counter(values)
    # Deterministic tie rule: lower count value wins.
    mode, frequency = max(counts.items(), key=lambda item: (item[1], -item[0]))
    return mode, frequency / len(values)


def _integer_signature(values: list[int]) -> bytes:
    if not values:
        return b""
    return struct.pack(f"<{len(values)}i", *values)


@dataclass(slots=True)
class LevelSummary:
    rt: array = field(default_factory=_doubles)
    tic: array = field(default_factory=_doubles)
    peak_count: array = field(default_factory=_doubles)
    charges: array = field(default_factory=_doubles)
    isolation_widths: array = field(default_factory=_doubles)
    isolation_precursor_mz: array = field(default_factory=_doubles)
    activation: Counter[tuple[str, str]] = field(default_factory=Counter)
    collision_energy: Counter[tuple[float, str]] = field(default_factory=Counter)
    invalid_peak_count: int = 0
    multiple_precursors: int = 0
    spectra_with_precursor: int = 0
    missing_activation: int = 0
    missing_collision_energy: int = 0
    observed_mz_low: float = math.inf
    observed_mz_high: float = -math.inf

    @property
    def count(self) -> int:
        return len(self.rt)

    def consume(self, spectrum: Spectrum) -> None:
        mz, intensity = spectrum.mz, spectrum.intensity
        if mz.ndim != 1 or intensity.ndim != 1 or mz.size != intensity.size:
            raise ValueError(f"Invalid peak-array dimensions in spectrum {spectrum.native_id!r}")
        valid = np.isfinite(mz) & (mz > 0) & np.isfinite(intensity) & (intensity >= 0)
        self.invalid_peak_count += int(mz.size - np.count_nonzero(valid))
        mz, intensity = mz[valid], intensity[valid]
        tic = float(intensity.sum(dtype=np.float64))
        self.rt.append(spectrum.rt)
        self.tic.append(tic)
        self.peak_count.append(len(mz))
        if mz.size:
            self.observed_mz_low = min(self.observed_mz_low, float(mz.min()))
            self.observed_mz_high = max(self.observed_mz_high, float(mz.max()))
        precursors = spectrum.precursors
        self.multiple_precursors += int(len(precursors) > 1)
        self.spectra_with_precursor += int(bool(precursors))
        # Per-scan precursor summaries use the first precursor, as in rawQC.
        # All precursor activation methods and energies still contribute metadata.
        first = precursors[0] if precursors else None
        self.charges.append(first.charge if first else 0)
        if first and first.isolation_width is not None:
            if math.isfinite(first.isolation_width) and first.isolation_width > 0:
                self.isolation_widths.append(first.isolation_width)
                if math.isfinite(first.mz) and first.mz > 0:
                    self.isolation_precursor_mz.append(first.mz)
        methods = {(t.accession, t.name) for p in precursors for t in p.activation}
        self.activation.update(methods)
        self.missing_activation += int(not methods)
        energies = {p.collision_energy for p in precursors if p.collision_energy is not None and math.isfinite(
            p.collision_energy[0],
        )}
        self.collision_energy.update(energies)
        self.missing_collision_energy += int(not energies)


@dataclass(slots=True)
class AcquisitionCycleTracker:
    """Compact MS1-delimited MS2-cycle evidence; no peak arrays are retained.

    A cycle starts at an MS1 spectrum and includes subsequent MS2 spectra up to
    the next MS1. MS2 spectra before the first MS1 are ignored. The final
    trailing cycle is included at calculation time when it contains MS2 scans.
    Target-set/order evidence is restricted to cycles where every MS2 scan has
    a finite positive first-precursor m/z paired with a positive isolation width.
    Targets are quantized to 0.1 Th, matching the run-level repetition metrics.
    """

    ms2_counts: list[int] = field(default_factory=list)
    unique_target_counts: list[int] = field(default_factory=list)
    target_set_signatures: Counter[bytes] = field(default_factory=Counter)
    target_order_signatures: Counter[bytes] = field(default_factory=Counter)
    target_observations: int = 0
    started: bool = False
    current_ms2_count: int = 0
    current_targets: list[int] = field(default_factory=list)

    def consume(self, spectrum: Spectrum) -> None:
        if spectrum.ms_level == 1:
            if self.started:
                self._finalize_current()
            self.started = True
            self.current_ms2_count = 0
            self.current_targets.clear()
            return
        if spectrum.ms_level != 2 or not self.started:
            return
        self.current_ms2_count += 1
        first = spectrum.precursors[0] if spectrum.precursors else None
        if (
            first
            and first.isolation_width is not None
            and math.isfinite(first.isolation_width)
            and first.isolation_width > 0
            and math.isfinite(first.mz)
            and first.mz > 0
        ):
            self.current_targets.append(int(np.rint(first.mz * 10.0)))

    def _finalize_current(self) -> None:
        if self.current_ms2_count <= 0:
            return
        self.ms2_counts.append(self.current_ms2_count)
        self.target_observations += len(self.current_targets)
        if len(self.current_targets) != self.current_ms2_count:
            return
        unique_targets = sorted(set(self.current_targets))
        self.unique_target_counts.append(len(unique_targets))
        self.target_set_signatures[_integer_signature(unique_targets)] += 1
        self.target_order_signatures[_integer_signature(self.current_targets)] += 1

    def metrics(self) -> dict[str, Any]:
        ms2_counts = list(self.ms2_counts)
        unique_counts = list(self.unique_target_counts)
        target_observations = self.target_observations

        if self.started and self.current_ms2_count > 0:
            ms2_counts.append(self.current_ms2_count)
            target_observations += len(self.current_targets)
            if len(self.current_targets) == self.current_ms2_count:
                unique_targets = sorted(set(self.current_targets))
                unique_counts.append(len(unique_targets))

        cycle_count = len(ms2_counts)
        eligible_count = len(unique_counts)
        total_ms2 = sum(ms2_counts)
        mode, _ = _mode_and_fraction(ms2_counts)
        ms2_array = np.asarray(ms2_counts, dtype=float)
        unique_array = np.asarray(unique_counts, dtype=float)

        return {
            "AcquisitionCycle_Count": cycle_count,
            "AcquisitionCycle_MS2Count_Quantiles": _quantiles(ms2_array),
            "AcquisitionCycle_MS2Count_Mode": mode,
            "AcquisitionCycle_TargetCoverageFraction": (
                target_observations / total_ms2 if total_ms2 else None
            ),
            "AcquisitionCycle_TargetEligibleFraction": (
                eligible_count / cycle_count if cycle_count else None
            ),
            "AcquisitionCycle_UniqueTargetCount_Quantiles": _quantiles(unique_array),
        }

    def classifier_consensus(self) -> tuple[float | None, float | None, float | None]:
        """Private evidence for DDA/DIA inference; not reported as QC metrics."""
        counts = list(self.ms2_counts)
        set_counts = self.target_set_signatures.copy()
        order_counts = self.target_order_signatures.copy()
        if self.started and self.current_ms2_count > 0:
            counts.append(self.current_ms2_count)
            if len(self.current_targets) == self.current_ms2_count:
                set_counts[_integer_signature(sorted(set(self.current_targets)))] += 1
                order_counts[_integer_signature(self.current_targets)] += 1
        _, modal_fraction = _mode_and_fraction(counts)
        eligible = sum(set_counts.values())
        return (
            modal_fraction,
            max(set_counts.values()) / eligible if eligible else None,
            max(order_counts.values()) / eligible if eligible else None,
        )


class RunSummary:
    """Accumulate compact scan summaries; spectra and peak arrays are not retained."""

    def __init__(self) -> None:
        self.levels: dict[int, LevelSummary] = {}
        self.acquisition_cycles = AcquisitionCycleTracker()
        self.faims: set[float] = set()
        # Boundaries, FAIMS CV and window ion mobility define a DIA window.
        # Count distinct geometries, never rounded precursor centres.
        self.dia_windows: set[tuple[float | None, ...]] = set()
        self.dia_window_incomplete = False
        self.dia_faims_present: set[bool] = set()
        self.dia_im_present: set[bool] = set()

    def _record_dia_window(self, spectrum: Spectrum) -> None:
        if spectrum.ms_level != 2:
            return
        self.dia_faims_present.add(spectrum.faims_cv is not None)
        if len(spectrum.precursors) != 1 or spectrum.unresolved_ion_mobility:
            self.dia_window_incomplete = True
            return
        precursor = spectrum.precursors[0]
        target, low, high = (
            precursor.mz, precursor.isolation_lower_offset,
            precursor.isolation_upper_offset,
        )
        im_low, im_high = precursor.ion_mobility_lower, precursor.ion_mobility_upper
        self.dia_im_present.add(im_low is not None or im_high is not None)
        if low is None or high is None:
            self.dia_window_incomplete = True
            return
        if (not all(math.isfinite(v) for v in (target, low, high))
            or target <= 0 or low < 0 or high < 0 or low + high <= 0
            or ((im_low is None) != (im_high is None))
            or (im_low is not None and im_high is not None
                and (not math.isfinite(im_low) or not math.isfinite(im_high)
                     or im_low >= im_high))
            or (spectrum.faims_cv is not None and not math.isfinite(spectrum.faims_cv))
        ):
            self.dia_window_incomplete = True
            return
        self.dia_windows.add(tuple(
            round(v, 6) if v is not None else None
            for v in (target - low, target + high, spectrum.faims_cv, im_low, im_high)
        ))

    def dia_window_metrics(self) -> dict[str, Any]:
        """Only complete, unambiguous MS2 geometry can produce PSI-MS DIA metrics."""
        if (not self.dia_windows or self.dia_window_incomplete
            or len(self.dia_faims_present) > 1 or len(self.dia_im_present) > 1):
            return {}
        widths: list[float] = []
        for window in self.dia_windows:
            if window[0] is None or window[1] is None:
                return {}
            widths.append(window[1] - window[0])
        if not widths:
            return {}
        return {
            "DIAIsolationWindowCount": len(self.dia_windows),
            "DIAIsolationWindowWidths": [min(widths), max(widths)],
        }

    def consume_spectrum(self, spectrum: Spectrum) -> None:
        if spectrum.ms_level < 1:
            raise ValueError(f"Invalid MS level: {spectrum.ms_level}")
        self.acquisition_cycles.consume(spectrum)
        self._record_dia_window(spectrum)
        if spectrum.ms_level not in self.levels:
            self.levels[spectrum.ms_level] = LevelSummary()
        self.levels[spectrum.ms_level].consume(spectrum)
        if spectrum.faims_cv is not None and math.isfinite(spectrum.faims_cv):
            self.faims.add(spectrum.faims_cv)

    def warnings(self) -> list[str]:
        messages = []
        for level, summary in sorted(self.levels.items()):
            invalid_rt = summary.count - _finite(summary.rt).size
            if invalid_rt:
                messages.append(
                    f"MS{level}: {invalid_rt} non-finite retention times excluded from RT metrics.",
                )
            if summary.invalid_peak_count:
                messages.append(
                    f"MS{level}: {summary.invalid_peak_count} invalid/negative peaks excluded.",
                )
            if summary.multiple_precursors:
                messages.append(
                    f"MS{level}: {summary.multiple_precursors} scans have multiple precursors; precursor metrics use the first.",
                )
        if not self.levels:
            messages.append("No spectra; spectrum-derived quantities are unavailable.")
        return messages


class QCMetricCalculator:
    """Pure calculations; no file access, CV serialization, or SDRF mutation."""

    def calculate(self, run: RunSummary, *, dia: bool = False) -> list[Metric]:
        values: dict[str, Any] = {}
        summaries = list(run.levels.values())
        ranges = [_range(_finite(s.rt)) for s in summaries]
        finite_ranges = [r for r in ranges if r is not None]
        values["NumberOfMSLevels"] = len(summaries)
        values["ChromatographyDuration"] = (
            max(r[1] for r in finite_ranges) - min(r[0] for r in finite_ranges)
            if finite_ranges else None
        )
        values["FAIMS_CV_Values"] = sorted(run.faims)
        values["FAIMS_CV_Count"] = len(run.faims)
        values["FAIMS_CV_Range"] = [min(run.faims), max(run.faims)] if run.faims else None
        for level in sorted(set(run.levels) | {1, 2}):
            values.update(self._level(level, run.levels.get(level, LevelSummary())))
        n1, n2 = values["NumberOfSpectra_MS1"], values["NumberOfSpectra_MS2"]
        values["MS1_to_MS2_Ratio"] = n1 / n2 if n2 else None
        values.update(run.acquisition_cycles.metrics())
        values.update(self._charge_metrics(run.levels.get(2, LevelSummary())))
        if dia:
            values.update(run.dia_window_metrics())
        return [Metric(key, value) for key, value in values.items()]

    def _level(self, level: int, summary: LevelSummary) -> dict[str, Any]:
        label = f"MS{level}"
        rt_all = np.asarray(summary.rt)
        tic_all = np.asarray(summary.tic)
        finite_rt = _finite(rt_all)
        mask = np.isfinite(rt_all) & np.isfinite(tic_all)
        rt, tic = rt_all[mask], tic_all[mask]
        if rt.size > 1 and np.any(np.diff(rt) < 0):
            order = np.argsort(rt, kind="stable")
            rt, tic = rt[order], tic[order]
        rts = np.sort(finite_rt) if np.any(np.diff(finite_rt) < 0) else finite_rt
        span = float(rts[-1] - rts[0]) if rts.size else 0.0
        result: dict[str, Any] = {
            f"NumberOfSpectra_{label}": summary.count,
            f"PeakDensity_{label}_Quantiles": _quantiles(np.asarray(summary.peak_count)),
            f"FastestFrequency_{label}": _fastest_frequency(rts),
            f"RT_{label}_Quantiles": _interval_fractions(rts, span),
            f"ObservedMzRange_{label}": (
                [summary.observed_mz_low, summary.observed_mz_high]
                if math.isfinite(summary.observed_mz_low)
                else None
            ),
        }
        if level == 1:
            differences = np.diff(rts)
            result["AvgCycleTime_MS1"] = _mean(differences[differences > 0])
            if rt.size > 1 and rt[-1] > rt[0] and tic.sum() > 0:
                cumulative = np.cumsum(tic)
                indices = np.searchsorted(cumulative, cumulative[-1] * np.array([0.25, 0.5, 0.75]))
                boundaries = np.concatenate(([rt[0]], rt[indices], [rt[-1]]))
                result["RT_TIC_Quantiles"] = (np.diff(boundaries) / (rt[-1] - rt[0])).tolist()
            else:
                result["RT_TIC_Quantiles"] = None
        return result

    def _charge_metrics(self, summary: LevelSummary) -> dict[str, Any]:
        charges = np.asarray(summary.charges)
        known = charges[charges > 0]
        bins = [int(np.count_nonzero(charges == k)) for k in range(1, 6)]
        bins.extend([int(np.count_nonzero(charges >= 6)), int(np.count_nonzero(charges <= 0))])
        return {
            "ChargeMean": _mean(known),
            "ChargeMedian": float(np.median(known)) if known.size else None,
            "ChargeRatio_3over2": bins[2] / bins[1] if bins[1] and bins[2] else None,
            "ChargeRatio_4over2": bins[3] / bins[1] if bins[1] and bins[3] else None,
        }
