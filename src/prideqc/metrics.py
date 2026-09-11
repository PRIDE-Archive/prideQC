"""ID-free metrics from compact per-scan summaries, inspired by rawQC.

Peak arrays are reduced once. Exact quantiles retain O(number of spectra)
scalar values, not O(total peaks) data. RT calculations use stable time order.
"""

from __future__ import annotations

import math
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


def _cv(values: FloatArray) -> float | None:
    return float(values.std(ddof=1) / values.mean()) if values.size > 1 and values.mean() else None


def _counter_table(counts: Counter[str], column: str) -> dict[str, list[Any]]:
    keys = sorted(counts)
    return {column: keys, "count": [counts[key] for key in keys]}


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


def _integral(rt: FloatArray, tic: FloatArray) -> float | None:
    if rt.size < 2:
        return None
    return float(np.sum(np.diff(rt) * (tic[1:] + tic[:-1]) * 0.5))


def _quartile_areas(rt: FloatArray, tic: FloatArray) -> list[float] | None:
    """Integrate piecewise-linear TIC, including the partial boundary segments.

    Unlike linear interpolation of cumulative area, this integrates the
    interpolated intensity exactly when a quartile cuts through a trapezoid.
    Duplicate RTs have zero area; interpolation uses the last value at that RT.
    """
    if rt.size < 2:
        return None
    bounds = np.quantile(rt, [0, 0.25, 0.5, 0.75, 1])
    dt = np.diff(rt)
    segments = (tic[1:] + tic[:-1]) * dt * 0.5
    cumulative = np.concatenate(([0.0], np.cumsum(segments)))
    areas = []
    for bound in bounds:
        index = min(int(np.searchsorted(rt, bound, side="right")) - 1, rt.size - 1)
        value = cumulative[index]
        if index < rt.size - 1 and dt[index] > 0:
            distance = bound - rt[index]
            slope = (tic[index + 1] - tic[index]) / dt[index]
            value += tic[index] * distance + 0.5 * slope * distance * distance
        areas.append(value)
    return np.diff(areas).tolist()


def _log_ratios(values: FloatArray) -> list[float | None] | None:
    if not values.size:
        return None
    q = np.quantile(values, [0.25, 0.5, 0.75, 1])
    return [float(
        np.log(b / a),
    ) if a > 0 and b > 0 else None for a, b in zip(
        q[:-1],
        q[1:],
        strict=True,
    )]


@dataclass(slots=True)
class LevelSummary:
    rt: array = field(default_factory=_doubles)
    tic: array = field(default_factory=_doubles)
    peak_count: array = field(default_factory=_doubles)
    base_peak: array = field(default_factory=_doubles)
    precursor_mz: array = field(default_factory=_doubles)
    precursor_intensity: array = field(default_factory=_doubles)
    charges: array = field(default_factory=_doubles)
    isolation_widths: array = field(default_factory=_doubles)
    polarity: Counter[str] = field(default_factory=Counter)
    representation: Counter[str] = field(default_factory=Counter)
    estimated_representation: Counter[str] = field(default_factory=Counter)
    activation: Counter[tuple[str, str]] = field(default_factory=Counter)
    collision_energy: Counter[tuple[float, str]] = field(default_factory=Counter)
    total_peaks: int = 0
    empty_scans: int = 0
    invalid_peak_count: int = 0
    multiple_precursors: int = 0
    spectra_with_precursor: int = 0
    missing_activation: int = 0
    missing_collision_energy: int = 0
    scan_low: float = math.inf
    scan_high: float = -math.inf
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
        self.base_peak.append(float(intensity.max()) if intensity.size else math.nan)
        self.total_peaks += len(mz)
        self.empty_scans += int(tic == 0)
        self.polarity[spectrum.polarity] += 1
        self.representation[spectrum.representation] += 1
        if spectrum.estimated_representation is not None:
            self.estimated_representation[spectrum.estimated_representation] += 1
        if mz.size:
            self.observed_mz_low = min(self.observed_mz_low, float(mz.min()))
            self.observed_mz_high = max(self.observed_mz_high, float(mz.max()))
        for lower, upper in spectrum.scan_windows:
            if math.isfinite(lower) and math.isfinite(upper) and 0 <= lower < upper:
                self.scan_low = min(self.scan_low, lower)
                self.scan_high = max(self.scan_high, upper)
        precursors = spectrum.precursors
        self.multiple_precursors += int(len(precursors) > 1)
        self.spectra_with_precursor += int(bool(precursors))
        # Per-scan precursor summaries use the first precursor, as in rawQC.
        # All precursor activation methods and energies still contribute metadata.
        first = precursors[0] if precursors else None
        self.charges.append(first.charge if first else 0)
        if first and math.isfinite(first.mz) and first.mz > 0:
            self.precursor_mz.append(first.mz)
        if first and math.isfinite(first.intensity) and first.intensity > 0:
            self.precursor_intensity.append(first.intensity)
        if first and first.isolation_width is not None:
            if math.isfinite(first.isolation_width) and first.isolation_width > 0:
                self.isolation_widths.append(first.isolation_width)
        methods = {(t.accession, t.name) for p in precursors for t in p.activation}
        self.activation.update(methods)
        self.missing_activation += int(not methods)
        energies = {p.collision_energy for p in precursors if p.collision_energy is not None and math.isfinite(
            p.collision_energy[0],
        )}
        self.collision_energy.update(energies)
        self.missing_collision_energy += int(not energies)


class RunSummary:
    """Accumulate scalars only; neither spectra nor chromatograms are retained."""

    def __init__(self) -> None:
        self.levels: dict[int, LevelSummary] = {}
        self.chromatograms: Counter[int] = Counter()
        self.chromatogram_points = 0
        self.chromatogram_low = math.inf
        self.chromatogram_high = -math.inf
        self.faims: set[float] = set()

    def consume_spectrum(self, spectrum: Spectrum) -> None:
        if spectrum.ms_level < 1:
            raise ValueError(f"Invalid MS level: {spectrum.ms_level}")
        if spectrum.ms_level not in self.levels:
            self.levels[spectrum.ms_level] = LevelSummary()
        self.levels[spectrum.ms_level].consume(spectrum)
        if spectrum.faims_cv is not None and math.isfinite(spectrum.faims_cv):
            self.faims.add(spectrum.faims_cv)

    def consume_chromatogram(self, rt: FloatArray, kind: int) -> None:
        self.chromatograms[kind] += 1
        self.chromatogram_points += len(rt)
        finite = _finite(rt)
        if finite.size:
            self.chromatogram_low = min(self.chromatogram_low, float(finite.min()))
            self.chromatogram_high = max(self.chromatogram_high, float(finite.max()))

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

    def calculate(self, run: RunSummary) -> list[Metric]:
        values: dict[str, Any] = {}
        summaries = list(run.levels.values())
        ranges = [_range(_finite(s.rt)) for s in summaries]
        finite_ranges = [r for r in ranges if r is not None]
        values["NumberOfMSLevels"] = len(summaries)
        values["NumberOfSpectra"] = sum(s.count for s in summaries)
        values["NumberOfSpectralPeaks"] = sum(s.total_peaks for s in summaries)
        values["ChromatographyDuration"] = (
            max(r[1] for r in finite_ranges) - min(r[0] for r in finite_ranges)
            if finite_ranges else None
        )
        values["NumberOfChromatograms"] = sum(run.chromatograms.values())
        values["NumberOfChromatogramDataPoints"] = run.chromatogram_points
        values["ChromatogramTypes"] = {"openms_type": sorted(
            run.chromatograms,
        ), "count": [run.chromatograms[k] for k in sorted(
            run.chromatograms,
        )]}
        values["ChromatogramRTRange"] = (
            [run.chromatogram_low, run.chromatogram_high]
            if math.isfinite(run.chromatogram_low) else None
        )
        maxima = [float(x.max()) for s in summaries if (x := _finite(s.base_peak)).size]
        values["BasePeak_All_Max"] = max(maxima) if maxima else None
        values["FAIMS_CV_Values"] = sorted(run.faims)
        values["FAIMS_CV_Count"] = len(run.faims)
        values["FAIMS_CV_Range"] = [min(run.faims), max(run.faims)] if run.faims else None
        for level in sorted(set(run.levels) | {1, 2}):
            values.update(self._level(level, run.levels.get(level, LevelSummary())))
        n1, n2 = values["NumberOfSpectra_MS1"], values["NumberOfSpectra_MS2"]
        values["MS1_to_MS2_Ratio"] = n1 / n2 if n2 else None
        values.update(self._charge_metrics(run.levels.get(2, LevelSummary())))
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
        finite_tic = _finite(tic_all)
        result: dict[str, Any] = {
            f"NumberOfSpectra_{label}": summary.count,
            f"RtRange_{label}": _range(rts),
            f"EmptyScans_{label}": summary.empty_scans,
            f"PeakDensity_{label}_Quantiles": _quantiles(np.asarray(summary.peak_count)),
            f"BasePeak_{label}_Mean": _mean(_finite(summary.base_peak)),
            f"TIC_{label}_Median": float(np.median(finite_tic)) if finite_tic.size else None,
            f"TIC_{label}_CV": _cv(finite_tic),
            f"TIC_{label}_Area": _integral(rt, tic),
            f"TIC_{label}_Area_RTQuantiles": _quartile_areas(rt, tic),
            f"TIC_{label}_QuartileLogRatios": _log_ratios(tic),
            f"TIC_{label}_ChangeQuartileLogRatios": _log_ratios(np.abs(np.diff(tic))),
            f"ScanRate_{label}": float(rts.size / span * 60) if span > 0 else None,
            f"FastestFrequency_{label}": _fastest_frequency(rts),
            f"RT_{label}_Quantiles": _interval_fractions(rts, span),
            f"ObservedMzRange_{label}": (
                [summary.observed_mz_low, summary.observed_mz_high]
                if math.isfinite(summary.observed_mz_low)
                else None
            ),
            f"ScanWindow_{label}": (
                [summary.scan_low, summary.scan_high]
                if math.isfinite(summary.scan_low)
                else None
            ),
            f"PeakTypes_{label}": _counter_table(summary.representation, "type"),
            f"EstimatedPeakTypes_{label}": _counter_table(
                summary.estimated_representation, "type"
            ),
            f"Polarity_{label}": _counter_table(summary.polarity, "polarity"),
            f"InvalidPeakCount_{label}": summary.invalid_peak_count,
        }
        previous, following = tic[:-1], tic[1:]
        valid = previous > 0
        ratios = following[valid] / previous[valid]
        result[f"TIC_{label}_SignalJump10x_Count"] = int(np.count_nonzero(ratios >= 10))
        result[f"TIC_{label}_SignalFall10x_Count"] = int(np.count_nonzero(ratios <= 0.1))
        result[f"RT_{label}_IQR"] = float(
            np.diff(np.quantile(rts, [0.25, 0.75]))[0],
        ) if rts.size else None
        if rts.size:
            low, high = np.quantile(rts, [0.25, 0.75])
            result[f"RT_{label}_IQRRate"] = float(
                np.count_nonzero((rts >= low) & (rts <= high)) / (high - low),
            ) if high > low else None
        else:
            result[f"RT_{label}_IQRRate"] = None
        if level == 1:
            differences = np.diff(rts)
            result["AvgCycleTime_MS1"] = _mean(differences[differences > 0])
            groups = np.sort(np.resize(np.arange(1, 5), tic.size))
            middle = tic[(groups == 2) | (groups == 3)]
            result["MedianTIC_in_RT_MS1_IQR"] = float(np.median(middle)) if middle.size else None
            if rt.size:
                half = (rt.size + 1) // 2
                spans = rt[half - 1:] - rt[:rt.size - half + 1]
                start = int(np.argmin(spans))
                result["TIC_MS1_MedianInHalfRange"] = float(np.median(tic[start:start + half]))
            else:
                result["TIC_MS1_MedianInHalfRange"] = None
            if rt.size > 1 and rt[-1] > rt[0] and tic.sum() > 0:
                cumulative = np.cumsum(tic)
                indices = np.searchsorted(cumulative, cumulative[-1] * np.array([0.25, 0.5, 0.75]))
                boundaries = np.concatenate(([rt[0]], rt[indices], [rt[-1]]))
                result["RT_TIC_Quantiles"] = (np.diff(boundaries) / (rt[-1] - rt[0])).tolist()
            else:
                result["RT_TIC_Quantiles"] = None
        if level > 1:
            mz = _finite(summary.precursor_mz)
            intensity = _finite(summary.precursor_intensity)
            widths = _finite(summary.isolation_widths)
            result.update({
                f"MzRange_{label}": _range(mz),
                f"PrecursorMz_{label}_Median": float(np.median(mz)) if mz.size else None,
                f"PrecursorIntensity_{label}_Quantiles": _quantiles(intensity),
                f"PrecursorIntensity_{label}_Mean": _mean(intensity),
                f"PrecursorIntensity_{label}_Sd": float(
                    intensity.std(ddof=1),
                ) if intensity.size > 1 else None,
                f"PrecursorIntensity_{label}_MissingCount": summary.count - intensity.size,
                f"IsolationWidth_{label}_Median": float(np.median(widths)) if widths.size else None,
                f"IsolationWidth_{label}_Quantiles": _quantiles(widths),
                f"IsolationWidth_{label}_Min": float(widths.min()) if widths.size else None,
                f"IsolationWidth_{label}_Max": float(widths.max()) if widths.size else None,
                f"IsolationWidth_{label}_Count": int(widths.size),
                f"IsolationWidth_{label}_FractionLe15": (
                    float(np.count_nonzero(widths <= 15) / widths.size)
                    if widths.size else None
                ),
                f"IsolationWidth_{label}_FractionGe15": (
                    float(np.count_nonzero(widths >= 15) / widths.size)
                    if widths.size else None
                ),
                f"MultiplePrecursors_{label}_Count": summary.multiple_precursors,
            })
            if level == 2:
                if intensity.size:
                    q5, q95 = np.quantile(intensity, [0.05, 0.95])
                    result["ExtentPrecursorIntensity_95over5_MS2"] = float(
                        q95 / q5,
                    ) if q5 > 0 else None
                else:
                    result["ExtentPrecursorIntensity_95over5_MS2"] = None
        return result

    def _charge_metrics(self, summary: LevelSummary) -> dict[str, Any]:
        charges = np.asarray(summary.charges)
        known = charges[charges > 0]
        bins = [int(np.count_nonzero(charges == k)) for k in range(1, 6)]
        bins.extend([int(np.count_nonzero(charges >= 6)), int(np.count_nonzero(charges <= 0))])
        return {
            "MS2_PrecursorCharge_Fractions": {
                "charge_state": ["1", "2", "3", "4", "5", ">=6", "unknown"],
                "count": bins,
                "fraction": [n / summary.count if summary.count else None for n in bins],
            },
            "ChargeMean": _mean(known),
            "ChargeMedian": float(np.median(known)) if known.size else None,
            "ChargeMin": int(known.min()) if known.size else None,
            "ChargeMax": int(known.max()) if known.size else None,
            "ChargeRatio_3over2": bins[2] / bins[1] if bins[1] and bins[2] else None,
            "ChargeRatio_4over2": bins[3] / bins[1] if bins[1] and bins[3] else None,
        }
