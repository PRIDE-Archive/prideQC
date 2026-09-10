"""Technical evidence and conservative, auditable SDRF suggestions."""

from __future__ import annotations

import math
from collections import Counter

import numpy as np

from prideqc.metrics import RunSummary
from prideqc.models import Annotation, CVTerm, EvidenceKind, RunMetadata, Spectrum

OBSERVED = EvidenceKind.OBSERVED
INFERRED = EvidenceKind.INFERRED
UNAVAILABLE = EvidenceKind.UNAVAILABLE


class TechnicalAnnotator:
    """Separate measured header/scan facts from acquisition-mode heuristics."""

    def annotate(self, metadata: RunMetadata, run: RunSummary) -> list[Annotation]:
        annotations = []
        for name, terms, column in (
            ("instrument", metadata.instruments, "comment[instrument]"),
            ("mass_analyzers", metadata.analyzers, None),
            ("ionization", metadata.ionization, None),
        ):
            annotations.append(Annotation(
                name, [{"accession": t.accession, "name": t.name} for t in terms] or None,
                OBSERVED if terms else UNAVAILABLE, "mzML instrumentConfiguration CV terms",
                "All configurations are reported; only an unambiguous instrument is proposed.",
                sdrf_column=column,
                sdrf_value=terms[0].sdrf_value() if len(terms) == 1 and column else None,
            ))
        annotations.append(Annotation(
            "instrument_details",
            metadata.instrument_details or None,
            OBSERVED if metadata.instrument_details else UNAVAILABLE,
            "OpenMS ExperimentalSettings instrument fields",
            "Plain instrument fields are preserved as observed metadata; no CV accession is inferred.",
        ))
        annotations.append(Annotation("instrument_serial_numbers", metadata.serial_numbers or None,
                                      OBSERVED if metadata.serial_numbers else UNAVAILABLE,
                                      "mzML instrument serial number"))
        for level, summary in sorted(run.levels.items()):
            if level < 2:
                continue
            methods = [
                {"accession": key[0], "name": key[1], "count": count}
                for key, count in sorted(summary.activation.items())
            ]
            only = next(iter(summary.activation)) if len(summary.activation) == 1 else None
            complete = summary.missing_activation == 0
            annotations.append(Annotation(
                f"dissociation_ms{level}", methods or None, OBSERVED if methods else UNAVAILABLE,
                "mzML precursor activation; scan counts per method",
                "Co-occurring methods are preserved; missing activation prevents automatic fill.",
                support=summary.count - summary.missing_activation, total=summary.count,
                sdrf_column="comment[dissociation method]" if level == 2 else None,
                sdrf_value=CVTerm(*only).sdrf_value() if only and complete and level == 2 else None,
            ))
            energies = [f"{value:g} {unit}" for value, unit in sorted(summary.collision_energy)]
            annotations.append(Annotation(
                f"collision_energy_ms{level}", energies or None,
                OBSERVED if energies else UNAVAILABLE, "mzML activation energy with preserved units",
                "Multiple recorded values do not establish stepped energy within each scan.",
                support=summary.count - summary.missing_collision_energy, total=summary.count,
                sdrf_column="comment[collision energy]" if level == 2 else None,
                sdrf_value=";".join(
                    energies,
                ) if energies and level == 2 and summary.missing_collision_energy == 0 else None,
            ))
        annotations.extend(self._acquisition(run))
        for name, column in (
            ("precursor_mass_tolerance", "comment[precursor mass tolerance]"),
            ("fragment_mass_tolerance", "comment[fragment mass tolerance]"),
        ):
            annotations.append(Annotation(
                name, None, UNAVAILABLE, "No calibrated mass-error estimator configured",
                "Search tolerances cannot be measured from instrument class or isolation width.",
                sdrf_column=column,
            ))
        return annotations

    def _acquisition(self, run: RunSummary) -> list[Annotation]:
        summary = run.levels.get(2)
        widths = np.asarray(summary.isolation_widths) if summary else np.array([])
        total = summary.count if summary else 0
        label, accession, support = None, None, 0
        # These thresholds retain TechSDRF's basic heuristic, but require coverage
        # and report ambiguity rather than treating isolation width as ground truth.
        if widths.size >= 100 and total and widths.size / total >= 0.9:
            narrow = int(np.count_nonzero(widths <= 3))
            wide = int(np.count_nonzero(widths >= 15))
            if narrow / widths.size >= 0.9:
                label, accession, support = "Data-dependent acquisition", "PRIDE:0000627", narrow
            elif wide / widths.size >= 0.9:
                label, accession, support = "Data-independent acquisition", "PRIDE:0000450", wide
        return [Annotation(
            "acquisition_method", label, INFERRED if label else UNAVAILABLE,
            "isolation-width heuristic v1: >=100 MS2, >=90% coverage and consensus",
            (
                "Narrow windows can also be PRM or narrow-window DIA. Mixed/intermediate "
                "windows and sparse runs abstain. Width unit is Th (m/z)."
            ),
            support=support, total=total,
            sdrf_column="comment[proteomics data acquisition method]",
            sdrf_value=CVTerm(accession, label).sdrf_value() if accession and label else None,
        )]


class DiagnosticIonCollector:
    """Optional lightweight signature screen; never asserts a PTM or label plex.

    Count centroid MS2/MS3 spectra with several reporter/diagnostic ions. No
    database search, mass-shift pairing, random null model, or stale dependency.
    """

    TARGETS = {
        "TMT_family": (126.127726, 127.131081, 128.134436, 129.137790, 130.141145, 131.138180),
        "iTRAQ_family": (114.1112, 115.1083, 116.1116, 117.1150),
        "glycan_oxonium": (204.0867, 366.1395),
    }

    def __init__(self, ppm: float = 20.0, min_relative_intensity: float = 0.01) -> None:
        if not math.isfinite(ppm) or ppm <= 0:
            raise ValueError("Diagnostic tolerance must be finite and positive.")
        if not 0 < min_relative_intensity <= 1:
            raise ValueError("Relative intensity threshold must be in (0, 1].")
        self.ppm = ppm
        self.min_relative_intensity = min_relative_intensity
        self.counts: Counter[str] = Counter()
        self.total = 0
        self.skipped_profile = 0

    def consume_spectrum(self, spectrum: Spectrum) -> None:
        if spectrum.ms_level not in (2, 3):
            return
        if spectrum.representation != "centroid":
            self.skipped_profile += 1
            return
        self.total += 1
        mz, intensity = spectrum.mz, spectrum.intensity
        valid = np.isfinite(mz) & (mz > 0) & np.isfinite(intensity) & (intensity > 0)
        mz, intensity = mz[valid], intensity[valid]
        if not mz.size:
            return
        mz = np.sort(mz[intensity >= intensity.max() * self.min_relative_intensity])
        for name, targets in self.TARGETS.items():
            masses = np.array(targets)
            delta = masses * self.ppm / 1e6
            matches = np.searchsorted(
                mz,
                masses + delta,
                side="right",
            ) > np.searchsorted(
                mz,
                masses - delta,
                side="left",
            )
            required = 2 if name == "glycan_oxonium" else 3
            self.counts[name] += int(np.count_nonzero(matches) >= required)

    def annotations(self) -> list[Annotation]:
        return [Annotation(
            f"signature_{name}", {"matching_spectra": self.counts[name], "eligible_spectra": self.total,
                                  "excluded_profile_or_unknown": self.skipped_profile},
            INFERRED if self.counts[name] else UNAVAILABLE,
            f"diagnostic ion screen v1; {self.ppm:g} ppm; relative intensity >= {self.min_relative_intensity:g}",
            (
                "Candidate evidence only; not proof of a modification, enrichment, labeling "
                "channel, or plex. Never written into SDRF automatically."
            ),
            support=self.counts[name], total=self.total,
        ) for name in self.TARGETS]
