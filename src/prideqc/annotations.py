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
        if metadata.acquisition_methods:
            # An explicit term may be recorded, but disagreement with a supported
            # scan-pattern inference blocks automatic classification/writeback.
            inferred_acquisition = next(
                (a for a in annotations if a.field == "acquisition_method"), None
            )
            annotations = [a for a in annotations if a.field not in {
                "acquisition_method", "acquisition_method_evidence"
            }]
            recorded_methods = list(dict.fromkeys(metadata.acquisition_methods))
            inferred_label = (
                str(inferred_acquisition.value)
                if inferred_acquisition is not None and inferred_acquisition.value is not None
                else None
            )
            conflicting = (
                len(recorded_methods) == 1
                and inferred_label is not None
                and inferred_label != recorded_methods[0].name
            )
            if conflicting:
                term = recorded_methods[0]
                annotations.append(Annotation(
                    "acquisition_method", None, UNAVAILABLE,
                    "conflict between recorded CV term and spectral inference",
                    "No automatic SDRF writeback when the RAW-header method conflicts with "
                    "the DDA/DIA isolation-pattern classifier.",
                    sdrf_column="comment[proteomics data acquisition method]",
                ))
                annotations.append(Annotation(
                    "acquisition_method_conflict",
                    {"recorded": {"accession": term.accession, "name": term.name},
                     "inferred": inferred_label},
                    INFERRED, "recorded acquisition versus DDA/DIA spectral heuristic",
                    "The recorded term and spectral inference disagree; human review required.",
                ))
            elif len(recorded_methods) == 1:
                term = recorded_methods[0]
                annotations.append(Annotation(
                    "acquisition_method", term.name, OBSERVED,
                    "explicit mzML fileContent PRIDE acquisition CV term",
                    "Recorded header metadata, not a scan-pattern inference.",
                    sdrf_column="comment[proteomics data acquisition method]",
                    sdrf_value=term.sdrf_value(),
                ))
                annotations.append(Annotation(
                    "acquisition_method_evidence",
                    {"cv_term": {"accession": term.accession, "name": term.name},
                     "classification": "OBSERVED_MZML_FILECONTENT",
                     "evidence_source": "mzML fileContent CV term"},
                    OBSERVED, "exact mzML fileContent CV term",
                    "The file records this method; no scan-pattern inference used.",
                ))
            else:
                annotations.append(Annotation(
                    "acquisition_method", None, UNAVAILABLE,
                    "conflicting mzML fileContent acquisition terms",
                    "Multiple acquisition terms prevent automatic classification or SDRF fill.",
                    sdrf_column="comment[proteomics data acquisition method]",
                ))
                annotations.append(Annotation(
                    "acquisition_recorded_terms",
                    [{"accession": t.accession, "name": t.name} for t in recorded_methods],
                    OBSERVED, "mzML fileContent CV terms",
                    "Conflict retained; no single acquisition method selected.",
                ))
        annotations.append(Annotation(
            "labeling_strategy", None, UNAVAILABLE,
            "No verified sample-channel label assignments in RAW spectrum metadata",
            "Reporter-ion matches are candidate evidence only; SILAC, label-free status, "
            "channel identities and multiplex design require submitter/SDRF evidence.",
            sdrf_column="comment[label]",
        ))
        annotations.append(Annotation(
            "chromatography_type", None, UNAVAILABLE,
            "Chromatography separation chemistry is not inferable from scan RT alone",
            "Run duration and scan counts are QC measurements, not evidence for "
            "reversed-phase, HILIC, SCX, LC flow regime or column type.",
        ))
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
        label, accession, support, support_total = None, None, 0, total
        method = (
            "acquisition heuristic v2: >=100 MS2, >=90% width coverage; fixed-cycle DIA "
            "requires >=20 cycles, >=8 MS2/cycle, >=90% cycle/set/order consensus, target-grid "
            "agreement, and median width >=5 Th; narrow DDA requires >=500 unique targets"
        )
        detail = (
            "High-specificity inference only. Wide-window consensus supports DIA. Stable repeated "
            "MS1-delimited cycles with non-narrow isolation also support DIA. Narrow isolation is "
            "called DDA only when many distinct precursor targets are observed; small fixed-target "
            "narrow runs abstain because PRM and narrow-window DIA are ambiguous. Width unit is Th (m/z)."
        )

        if summary is not None and widths.size >= 100 and total and widths.size / total >= 0.9:
            narrow = int(np.count_nonzero(widths <= 3))
            wide = int(np.count_nonzero(widths >= 15))
            narrow_fraction = narrow / widths.size
            wide_fraction = wide / widths.size
            median_width = float(np.median(widths))

            cycle = run.acquisition_cycles.classifier_features()
            cycle_count = int(cycle["count"] or 0)
            cycle_mode = int(cycle["mode"] or 0)
            target_coverage = cycle["target_coverage"]
            target_eligible_fraction = cycle["eligible_fraction"]
            cycle_modal_fraction, set_modal_fraction, order_modal_fraction = (
                run.acquisition_cycles.classifier_consensus()
            )

            rounded_targets = np.round(np.asarray(summary.isolation_precursor_mz), 1)
            unique_targets = int(np.unique(rounded_targets).size) if rounded_targets.size else 0
            target_grid_tolerance = max(1, int(round(cycle_mode * 0.1)))

            stable_cycle = (
                cycle_count >= 20
                and cycle_mode >= 8
                and 8 <= unique_targets <= 200
                and abs(unique_targets - cycle_mode) <= target_grid_tolerance
                and cycle_modal_fraction is not None and cycle_modal_fraction >= 0.9
                and target_coverage is not None and target_coverage >= 0.9
                and target_eligible_fraction is not None and target_eligible_fraction >= 0.9
                and set_modal_fraction is not None and set_modal_fraction >= 0.9
                and order_modal_fraction is not None and order_modal_fraction >= 0.9
            )

            if wide_fraction >= 0.9:
                label, accession, support = "Data-independent acquisition", "PRIDE:0000450", wide
            elif stable_cycle and median_width >= 5 and order_modal_fraction is not None:
                label, accession = "Data-independent acquisition", "PRIDE:0000450"
                modal_cycles = int(round(order_modal_fraction * cycle_count))
                support_total = total
                support = min(total, modal_cycles * cycle_mode)
            elif narrow_fraction >= 0.9 and unique_targets >= 500:
                label, accession, support = "Data-dependent acquisition", "PRIDE:0000627", narrow

        annotations = [Annotation(
            "acquisition_method", label, INFERRED if label else UNAVAILABLE,
            method,
            detail,
            support=support, total=support_total,
            sdrf_column="comment[proteomics data acquisition method]",
            sdrf_value=CVTerm(accession, label).sdrf_value() if accession and label else None,
        )]
        if accession and label:
            annotations.append(Annotation(
                "acquisition_method_evidence",
                {
                    "cv_term": {"accession": accession, "name": label},
                    "evidence_source": "MS2 precursor isolation widths and acquisition-cycle structure",
                    "supporting_scans": support,
                    "assessed_ms2_scans": support_total,
                    "classification": "INFERRED_NOT_OBSERVED",
                    "classifier": "acquisition heuristic v2",
                },
                INFERRED, method,
                "A data-derived run classification; not a recorded instrument acquisition setting. "
                "PRM, SRM/MRM and narrow DIA are not inferred from narrow-window scans alone.",
                support=support, total=support_total,
            ))
        return annotations


def labeling_candidates(evidence: list[Annotation]) -> Annotation | None:
    """Non-exclusive reporter-based candidate terms; never a verified sample label."""
    candidates: list[dict[str, object]] = []
    tmt = next((a for a in evidence if a.field in {
        "tmt_reporter_evidence", "tmt_profile_reporter_evidence"
    } and a.kind == INFERRED and isinstance(a.value, dict)), None)
    if tmt is not None:
        levels = tmt.value.get("levels", {})
        matched: dict[str, int] = {}
        if isinstance(levels, dict):
            for level in levels.values():
                if not isinstance(level, dict):
                    continue
                methods = level.get("methods", {})
                if not isinstance(methods, dict):
                    continue
                for name, info in methods.items():
                    if isinstance(info, dict):
                        count = info.get("reporter_positive_spectra", 0)
                        if isinstance(count, int) and count > 0:
                            matched[name] = matched.get(name, 0) + count
        if matched:
            candidates.append({
                "cv_term": {"accession": "PRIDE:0000314", "name": "TMT"},
                "source": tmt.field,
                "reporter_positive_spectra_per_nonexclusive_method": matched,
            })
    if not any(
        item.get("source") in {"tmt_reporter_evidence", "tmt_profile_reporter_evidence"}
        for item in candidates
    ):
        tmt_signature = next(
            (a for a in evidence if a.field == "signature_TMT_family"
             and a.kind == INFERRED and a.support is not None and a.support > 0),
            None,
        )
        if tmt_signature is not None:
            candidates.append({
                "cv_term": {"accession": "PRIDE:0000314", "name": "TMT"},
                "source": tmt_signature.field,
                "matching_spectra": tmt_signature.support,
                "eligible_spectra": tmt_signature.total,
            })
    itraq = next((a for a in evidence if a.field == "signature_iTRAQ_family"
                 and a.kind == INFERRED and a.support is not None and a.support > 0), None)
    if itraq is not None:
        candidates.append({
            "cv_term": {"accession": "PRIDE:0000313", "name": "iTRAQ"},
            "source": itraq.field,
            "matching_spectra": itraq.support,
            "eligible_spectra": itraq.total,
        })
    if not candidates:
        return None
    return Annotation(
        "labeling_strategy_candidates",
        {"scientific_status": "CANDIDATE_ONLY_NOT_VERIFIED_SAMPLE_LABEL",
         "candidate_terms": candidates,
         "not_inferred": ["SILAC", "label-free", "TMT plex", "sample-to-channel assignment"]},
        INFERRED,
        "Optional ID-free reporter-ion collectors; non-exclusive candidate vocabulary mappings",
        "Neither the chemical labelling strategy nor the biological sample/channel assignment "
        "is established by reporter signals alone. No SDRF writeback.",
    )


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
