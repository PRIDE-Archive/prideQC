"""Experimental one-pass TMT evidence from explicitly tagged profile MS2/MS3.

Uses OpenMS PeakPickerHiRes on ephemeral reporter-region arrays. Native centroid
TMT evidence remains a separate collector and is never modified by this module.
"""

from __future__ import annotations

from dataclasses import replace
from typing import Any

import numpy as np

from prideqc.models import Annotation, EvidenceKind, Spectrum
from prideqc.tmt import TMTReporterCollector


class TMTProfileReporterCollector:
    """Profile-only reporter evidence, deliberately separate from centroid QC.

    Peak picking is restricted to 125--137 m/z (reporter catalog region) with
    a 1 Da safety margin. It is experimental and does not validate the true
    labeling kit, channel occupancy, isotope corrections or SDRF assignments.
    """

    def __init__(
        self,
        *,
        oms: Any | None = None,
        ppm: float = 10.0,
        min_relative_intensity: float = 0.001,
        minimum_channel_hits: int = 3,
    ) -> None:
        if oms is None:
            import pyopenms as runtime_oms

            # pyOpenMS exports compiled extension classes dynamically. Its
            # typing information may not declare these runtime members.
            runtime_module: Any = runtime_oms
        else:
            runtime_module = oms
        self._oms: Any = runtime_module

        # Reuse runtime channel tables and matching algorithm. Each accepted
        # picked profile scan is delegated as an ephemeral peak list, not as an
        # observed vendor centroid scan. Output is explicitly relabeled below.
        self._matcher = TMTReporterCollector(
            oms=self._oms,
            ppm=ppm,
            min_relative_intensity=min_relative_intensity,
            minimum_channel_hits=minimum_channel_hits,
        )
        self._picker = self._oms.PeakPickerHiRes()
        masses = [
            channel.mz
            for level in self._matcher._levels.values()
            for method in level["methods"].values()
            for channel in method["channels"]
        ]
        self._range = (min(masses) - 1.0, max(masses) + 1.0)
        self._observable: dict[tuple[int, str, str], int] = {}
        self._stats: dict[int, dict[str, int]] = {
            level: {
                "profile_spectra_scanned": 0,
                "profile_spectra_evaluable": 0,
                "profile_spectra_without_reporter_region": 0,
                "profile_spectra_invalid_samples": 0,
                "profile_peak_picker_failures": 0,
                "profile_spectra_with_no_picked_peaks": 0,
                "non_profile_spectra_excluded": 0,
            }
            for level in (2, 3)
        }

    def consume_spectrum(self, spectrum: Spectrum) -> None:
        stats = self._stats.get(spectrum.ms_level)
        if stats is None:
            return
        if spectrum.representation != "profile":
            stats["non_profile_spectra_excluded"] += 1
            return
        stats["profile_spectra_scanned"] += 1
        mz = np.asarray(spectrum.mz, dtype=float)
        intensity = np.asarray(spectrum.intensity, dtype=float)
        if mz.ndim != 1 or intensity.shape != mz.shape or mz.size < 3:
            stats["profile_spectra_invalid_samples"] += 1
            return
        valid = np.isfinite(mz) & (mz > 0) & np.isfinite(intensity) & (intensity >= 0)
        if np.count_nonzero(valid) < 3:
            stats["profile_spectra_invalid_samples"] += 1
            return
        mz, intensity = mz[valid], intensity[valid]
        # Reject duplicate or unsorted profile axes; do not fabricate resolution.
        if np.any(np.diff(mz) <= 0):
            stats["profile_spectra_invalid_samples"] += 1
            return
        left, right = self._range
        region = (mz >= left) & (mz <= right)
        if np.count_nonzero(region) < 3:
            stats["profile_spectra_without_reporter_region"] += 1
            return
        mz_region, intensity_region = mz[region], intensity[region]
        sampled_channels: list[tuple[int, str, str]] = []
        sampled_masses: list[float] = []
        for method_name, method in self._matcher._levels[spectrum.ms_level]["methods"].items():
            for channel in method["channels"]:
                lo = int(np.searchsorted(mz_region, channel.mz - 0.02))
                hi = int(np.searchsorted(mz_region, channel.mz + 0.02, side="right"))
                if hi - lo >= 3:
                    sampled_channels.append((spectrum.ms_level, method_name, channel.name))
                    sampled_masses.append(channel.mz)
        if not sampled_channels:
            stats["profile_spectra_without_reporter_region"] += 1
            return
        source = self._oms.MSSpectrum()
        source.set_peaks((mz_region, intensity_region))
        picked = self._oms.MSSpectrum()
        try:
            self._picker.pick(source, picked)
            picked_mz, picked_intensity = picked.get_peaks()
        except (RuntimeError, ValueError, TypeError, AttributeError):
            stats["profile_peak_picker_failures"] += 1
            return
        picked_mz = np.asarray(picked_mz, dtype=float)
        picked_intensity = np.asarray(picked_intensity, dtype=float)
        if picked_mz.ndim != 1 or picked_mz.shape != picked_intensity.shape:
            stats["profile_peak_picker_failures"] += 1
            return
        # No picked reporter-region maxima is a measurable absence of signal,
        # NOT a peak-picker failure, if the source contained valid profile data.
        stats["profile_spectra_evaluable"] += 1
        # Per-channel scan-range accounting avoids calling an inaccessible
        # reporter channel "unused" when the acquisition excludes its m/z.
        for key in sampled_channels:
            self._observable[key] = self._observable.get(key, 0) + 1
        if picked_mz.size:
            # Match the centroid collector's full-scan base-peak threshold;
            # restricting peak picking to reporter m/z must not quietly lower
            # the relative-intensity cutoff to the local reporter-window base.
            cutoff = float(np.max(intensity)) * self._matcher.min_relative_intensity
            retained = (
                np.isfinite(picked_mz)
                & np.isfinite(picked_intensity)
                & (picked_intensity >= cutoff)
                & (picked_intensity > 0)
            )
            # Require actual local profile sampling for every matched reporter.
            # Otherwise a partially acquired reporter range could manufacture
            # apparent evidence for a channel with no usable m/z grid points.
            sampled = np.zeros(len(picked_mz), dtype=bool)
            for target in set(sampled_masses):
                sampled |= np.abs(picked_mz - target) <= target * self._matcher.ppm / 1e6
            retained &= sampled
            picked_mz, picked_intensity = picked_mz[retained], picked_intensity[retained]
        if not picked_mz.size:
            stats["profile_spectra_with_no_picked_peaks"] += 1
        self._matcher.consume_spectrum(
            replace(spectrum, mz=picked_mz, intensity=picked_intensity, representation="centroid")
        )

    def annotations(self) -> list[Annotation]:
        original = self._matcher.annotations()[0]
        values = original.value
        levels = {}
        evaluable_total = 0
        support = 0
        for ms_level, stats in self._stats.items():
            level = values["levels"][f"MS{ms_level}"]
            eligible = stats["profile_spectra_evaluable"]
            evaluable_total += eligible
            methods = {}
            for name, method in level["methods"].items():
                hits = method["reporter_positive_spectra"]
                support += hits  # Overlapping catalogs; not a unique-spectrum count.
                channels = []
                for channel in method["channels"]:
                    renamed = dict(channel)
                    fraction = renamed.pop("fraction_of_centroid_spectra")
                    renamed["fraction_of_evaluable_profile_spectra"] = fraction
                    observed = self._observable.get((ms_level, name, channel["channel"]), 0)
                    renamed["profile_spectra_with_local_mz_sampling"] = observed
                    renamed["fraction_of_locally_sampled_profile_spectra"] = (
                        renamed["spectra_with_hit"] / observed if observed else None
                    )
                    channels.append(renamed)
                methods[name] = {
                    "reporter_positive_spectra": hits,
                    "positive_fraction_of_evaluable_profile": (
                        hits / eligible if eligible else None
                    ),
                    "channels": channels,
                }
            levels[f"MS{ms_level}"] = {**stats, "methods": methods}
        evidence = {
            "schema_version": 1,
            "scientific_status": "EXPERIMENTAL_PROFILE_PEAK_PICKING_NOT_FOR_SDRF_WRITEBACK",
            "peak_picking_method": "OpenMS PeakPickerHiRes, ephemeral 1 Da padded reporter region",
            "profile_peak_picking_mz_range": [float(x) for x in self._range],
            "pyopenms_version": self._matcher.pyopenms_version,
            "channel_catalog_source": values["channel_catalog_source"],
            "parameters": values["parameters"],
            "unavailable_methods": values["unavailable_methods"],
            "levels": levels,
        }
        return [
            Annotation(
                "tmt_profile_reporter_evidence",
                evidence,
                EvidenceKind.INFERRED
                if support
                else (EvidenceKind.OBSERVED if evaluable_total else EvidenceKind.UNAVAILABLE),
                "one-pass experimental OpenMS PeakPickerHiRes profile TMT reporter screen v1",
                "Profile-derived channel evidence only; not native centroid measurements. "
                "No TMT kit, active sample-channel occupancy or SDRF labels inferred.",
                support=None,
                total=evaluable_total,
            )
        ]
