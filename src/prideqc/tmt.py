"""One-pass, identification-free TMT reporter-ion evidence (not sample labels).

Channel m/z values are taken from the running pyOpenMS quantitation methods.
No library of hard-coded channel masses, isotope correction, plex classification,
or SDRF label writeback is performed.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any

import numpy as np

from prideqc.models import Annotation, EvidenceKind, Spectrum

TMT_METHOD_CLASSES = {
    "TMT6": "TMTSixPlexQuantitationMethod",
    "TMT10": "TMTTenPlexQuantitationMethod",
    "TMT11": "TMTElevenPlexQuantitationMethod",
    "TMTpro16": "TMTSixteenPlexQuantitationMethod",
    "TMTpro18": "TMTEighteenPlexQuantitationMethod",
}


@dataclass(slots=True)
class _Channel:
    name: str
    mz: float
    spectra_with_hit: int = 0
    intensity_total: float = 0.0
    absolute_ppm_total: float = 0.0

    def report(self, centroid_spectra: int) -> dict[str, Any]:
        hits = self.spectra_with_hit
        return {
            "channel": self.name,
            "expected_mz": self.mz,
            "spectra_with_hit": hits,
            "fraction_of_centroid_spectra": hits / centroid_spectra if centroid_spectra else None,
            "uncorrected_total_matched_intensity": self.intensity_total,
            "mean_absolute_mass_error_ppm": self.absolute_ppm_total / hits if hits else None,
        }


class TMTReporterCollector:
    """Report channel-specific MS2/MS3 centroid peak evidence in one reader pass.

    Calls are optional, enabled via ``--tmt-reporter-evidence``. Data from profile
    or unknown spectrum types are *excluded*, never reported as TMT-negative.
    Matched peaks are not corrected for isotope impurity or interference and
    shared reporter masses cannot uniquely identify a labeling kit or occupancy.
    """

    def __init__(
        self,
        *,
        ppm: float = 10.0,
        min_relative_intensity: float = 0.001,
        minimum_channel_hits: int = 3,
        oms: Any | None = None,
    ) -> None:
        if not math.isfinite(ppm) or ppm <= 0:
            raise ValueError("TMT reporter tolerance must be finite and positive")
        if not math.isfinite(min_relative_intensity) or not 0 < min_relative_intensity <= 1:
            raise ValueError("TMT relative intensity must be finite and in (0, 1]")
        if minimum_channel_hits < 1:
            raise ValueError("minimum_channel_hits must be positive")
        if oms is None:
            import pyopenms as oms  # Imported only when the collector is enabled.

        self.ppm = ppm
        self.min_relative_intensity = min_relative_intensity
        self.minimum_channel_hits = minimum_channel_hits
        self.pyopenms_version = str(getattr(oms, "__version__", "unknown"))
        self.catalog_unavailable: dict[str, str] = {}
        catalog: dict[str, list[tuple[str, float]]] = {}
        for method, class_name in TMT_METHOD_CLASSES.items():
            cls = getattr(oms, class_name, None)
            if cls is None:
                self.catalog_unavailable[method] = f"missing pyOpenMS class {class_name}"
                continue
            try:
                channels = [
                    (str(channel.name), float(channel.center))
                    for channel in cls().getChannelInformation()
                ]
                if not channels or any(
                    not name or not math.isfinite(mz) or mz <= 0 for name, mz in channels
                ):
                    raise ValueError("missing or invalid channel name/mass")
                if len({name for name, _ in channels}) != len(channels):
                    raise ValueError("duplicate channel name")
                catalog[method] = channels
            except (AttributeError, TypeError, ValueError, RuntimeError) as exc:
                self.catalog_unavailable[method] = f"{type(exc).__name__}: {exc}"
        if not catalog:
            raise RuntimeError(
                "No valid TMT reporter channel tables in this pyOpenMS build: "
                f"{self.catalog_unavailable}"
            )

        self._levels: dict[int, dict[str, Any]] = {}
        for level in (2, 3):
            self._levels[level] = {
                "spectra_scanned": 0,
                "centroid_spectra": 0,
                "profile_skipped": 0,
                "unknown_skipped": 0,
                "methods": {
                    method: {
                        "reporter_positive_spectra": 0,
                        "channels": [_Channel(name, mz) for name, mz in channels],
                    }
                    for method, channels in catalog.items()
                },
            }

    def consume_spectrum(self, spectrum: Spectrum) -> None:
        level = self._levels.get(spectrum.ms_level)
        if level is None:
            return
        level["spectra_scanned"] += 1
        if spectrum.representation != "centroid":
            skip_name = (
                "profile_skipped" if spectrum.representation == "profile" else "unknown_skipped"
            )
            level[skip_name] += 1
            return
        level["centroid_spectra"] += 1
        mz = np.asarray(spectrum.mz, dtype=float)
        intensity = np.asarray(spectrum.intensity, dtype=float)
        if mz.shape != intensity.shape or mz.ndim != 1:
            raise ValueError("Spectrum m/z and intensity must be 1D arrays of equal length")
        valid = np.isfinite(mz) & (mz > 0) & np.isfinite(intensity) & (intensity > 0)
        if not np.any(valid):
            return
        mz, intensity = mz[valid], intensity[valid]
        keep = intensity >= float(np.max(intensity)) * self.min_relative_intensity
        mz, intensity = mz[keep], intensity[keep]
        if not len(mz):
            return
        order = np.argsort(mz)
        mz, intensity = mz[order], intensity[order]
        for method in level["methods"].values():
            matched = 0
            for channel in method["channels"]:
                margin = channel.mz * self.ppm / 1e6
                left = int(np.searchsorted(mz, channel.mz - margin, side="left"))
                right = int(np.searchsorted(mz, channel.mz + margin, side="right"))
                if left == right:
                    continue
                local_index = int(np.argmax(intensity[left:right])) + left
                observed = float(mz[local_index])
                channel.spectra_with_hit += 1
                channel.intensity_total += float(intensity[local_index])
                channel.absolute_ppm_total += abs(observed - channel.mz) / channel.mz * 1e6
                matched += 1
            if matched >= self.minimum_channel_hits:
                method["reporter_positive_spectra"] += 1

    def annotations(self) -> list[Annotation]:
        levels: dict[str, Any] = {}
        support = 0
        centroid_total = 0
        for ms_level, level in self._levels.items():
            centroid = level["centroid_spectra"]
            centroid_total += centroid
            methods = {}
            for name, method in level["methods"].items():
                hits = method["reporter_positive_spectra"]
                support += hits  # May overlap across catalogs; NOT an independent spectrum count.
                methods[name] = {
                    "reporter_positive_spectra": hits,
                    "positive_fraction_of_centroid": hits / centroid if centroid else None,
                    "channels": [channel.report(centroid) for channel in method["channels"]],
                }
            levels[f"MS{ms_level}"] = {
                "spectra_scanned": level["spectra_scanned"],
                "centroid_spectra": centroid,
                "profile_skipped": level["profile_skipped"],
                "unknown_skipped": level["unknown_skipped"],
                "methods": methods,
            }
        evidence = {
            "schema_version": 1,
            "scientific_status": "EXPLORATORY_NOT_FOR_AUTOMATIC_SDRF_WRITEBACK",
            "channel_catalog_source": (
                "runtime pyOpenMS quantitation method getChannelInformation()"
            ),
            "pyopenms_version": self.pyopenms_version,
            "parameters": {
                "ppm": self.ppm,
                "min_relative_intensity_fraction_of_spectrum_base_peak": (
                    self.min_relative_intensity
                ),
                "minimum_channel_hits_for_reporter_positive_spectrum": self.minimum_channel_hits,
            },
            "unavailable_methods": self.catalog_unavailable,
            "levels": levels,
        }
        return [
            Annotation(
                "tmt_reporter_evidence",
                evidence,
                EvidenceKind.INFERRED
                if support
                else (EvidenceKind.OBSERVED if centroid_total else EvidenceKind.UNAVAILABLE),
                "one-pass ID-free centroid MS2/MS3 TMT reporter peak screen v1",
                (
                    "Shared reporter masses make method scores non-exclusive. "
                    "Not isotope-corrected; "
                    "profile/unknown spectra are excluded. No TMT kit, channel occupancy, "
                    "sample mapping, or SDRF labels are inferred."
                ),
                support=None,  # Method counts overlap; do not sum them as unique spectra.
                total=centroid_total,
            )
        ]
