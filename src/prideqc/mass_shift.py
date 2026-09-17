"""Fast identification-free recurrent precursor mass-shift scouting.

The collector in this module is deliberately not a peptide/PTM search engine.
It first requires fragment-spectrum relatedness, then clusters recurring absolute
neutral precursor-mass differences and finally annotates those clusters against
the modification knowledge shipped with the pinned OpenMS/pyOpenMS runtime.

Peak arrays remain ephemeral. Only compact top-peak fingerprints and a bounded
inverted index are retained. The resulting evidence is suitable for QC and
hypothesis generation; it does not identify peptide sequences or localize a
modification site.
"""

from __future__ import annotations

import math
from collections import Counter, defaultdict, deque
from collections.abc import Iterable
from dataclasses import dataclass, field
from statistics import median
from typing import TYPE_CHECKING, Any

import numpy as np

from prideqc.models import Annotation, EvidenceKind, Spectrum

if TYPE_CHECKING:
    from prideqc.mass_error import RepeatSpectrumMassErrorCollector

# CODATA/OpenMS-compatible fallbacks are used only when a caller does not supply
# the constants and pyOpenMS is unavailable (principally pure unit tests).
_FALLBACK_PROTON_MASS_U = 1.007276466621
_FALLBACK_C13C12_MASSDIFF_U = 1.00335483507

# Common positive-mode adduct substitutions. These are classification guards,
# not modification identifications.
_ADDUCT_SHIFTS_DA = {
    "sodium-for-proton": 21.981943,
    "potassium-for-proton": 37.955882,
}


@dataclass(frozen=True, slots=True)
class ModificationRecord:
    """One mass-compatible modification entry derived from OpenMS/UniMod."""

    accession: str
    name: str
    delta_mass_da: float
    origins: tuple[str, ...] = ()
    term_specificities: tuple[str, ...] = ()
    source_classification: str = ""


@dataclass(frozen=True, slots=True)
class _IndexedSpectrum:
    identifier: int
    rt: float
    precursor_mz: float
    charge: int
    neutral_mass: float
    mz: np.ndarray
    intensity: np.ndarray
    bins: tuple[int, ...]


@dataclass(frozen=True, slots=True)
class _MassShiftObservation:
    delta_mass_da: float
    left_id: int
    right_id: int
    mean_neutral_mass: float
    similarity: float


@dataclass(frozen=True, slots=True)
class _Cluster:
    center_da: float
    sigma_da: float
    minimum_da: float
    maximum_da: float
    pair_support: int
    unique_spectrum_support: int
    median_similarity: float
    mean_neutral_mass: float


def _text(value: Any) -> str:
    if isinstance(value, bytes):
        return value.decode(errors="replace").strip()
    return str(value).strip()


def _openms_constants(oms: Any | None = None) -> tuple[float, float]:
    """Return proton and C13-C12 masses from the active pyOpenMS build."""

    active_oms: Any
    if oms is None:
        try:
            import pyopenms
        except ImportError:
            return _FALLBACK_PROTON_MASS_U, _FALLBACK_C13C12_MASSDIFF_U
        active_oms = pyopenms
    else:
        active_oms = oms
    constants = getattr(active_oms, "Constants", None)
    proton = float(getattr(constants, "PROTON_MASS_U", _FALLBACK_PROTON_MASS_U))
    isotope = float(
        getattr(constants, "C13C12_MASSDIFF_U", _FALLBACK_C13C12_MASSDIFF_U)
    )
    return proton, isotope


def load_openms_modifications(oms: Any | None = None) -> tuple[ModificationRecord, ...]:
    """Freeze the UniMod-bearing entries exposed by the active OpenMS runtime.

    OpenMS can contain one record per residue/terminal specificity for the same
    UniMod concept. Those are collapsed to one record per accession/name/mass
    while preserving all observed origins and terminal specificities.
    """

    active_oms: Any
    if oms is None:
        import pyopenms

        active_oms = pyopenms
    else:
        active_oms = oms

    database = active_oms.ModificationsDB()
    grouped: dict[tuple[str, str, float], dict[str, Any]] = {}
    for index in range(int(database.getNumberOfModifications())):
        modification = database.getModification(index)
        accession = _text(modification.getUniModAccession())
        if not accession or not accession.casefold().startswith("unimod:"):
            continue
        delta = float(modification.getDiffMonoMass())
        if not math.isfinite(delta) or delta == 0.0:
            continue
        name = _text(modification.getFullName()) or _text(modification.getId())
        key = (accession, name, round(delta, 9))
        record = grouped.setdefault(
            key,
            {
                "accession": accession,
                "name": name,
                "delta": delta,
                "origins": set(),
                "terms": set(),
                "classification": "",
            },
        )
        try:
            origin = _text(modification.getOrigin())
            if origin and origin not in {"X", "."}:
                record["origins"].add(origin)
        except (AttributeError, RuntimeError, TypeError, ValueError):
            pass
        try:
            specificity = modification.getTermSpecificity()
            specificity_name = _text(modification.getTermSpecificityName(specificity))
            if specificity_name:
                record["terms"].add(specificity_name)
        except (AttributeError, RuntimeError, TypeError, ValueError):
            pass
        try:
            classification = modification.getSourceClassification()
            classification_name = _text(
                modification.getSourceClassificationName(classification)
            )
            if classification_name:
                record["classification"] = classification_name
        except (AttributeError, RuntimeError, TypeError, ValueError):
            pass

    records = [
        ModificationRecord(
            accession=value["accession"],
            name=value["name"],
            delta_mass_da=float(value["delta"]),
            origins=tuple(sorted(value["origins"])),
            term_specificities=tuple(sorted(value["terms"])),
            source_classification=str(value["classification"]),
        )
        for value in grouped.values()
    ]
    return tuple(
        sorted(records, key=lambda item: (abs(item.delta_mass_da), item.accession, item.name))
    )


def _top_peaks(
    spectrum: Spectrum,
    limit: int,
    *,
    precursor_exclusion_da: float,
) -> tuple[np.ndarray, np.ndarray]:
    mz = np.asarray(spectrum.mz, dtype=float)
    intensity = np.asarray(spectrum.intensity, dtype=float)
    valid = np.isfinite(mz) & (mz > 0) & np.isfinite(intensity) & (intensity > 0)
    first = spectrum.precursors[0] if spectrum.precursors else None
    if first is not None and math.isfinite(first.mz) and first.mz > 0:
        valid &= np.abs(mz - first.mz) > precursor_exclusion_da
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


def _matched_intensity_dot(
    left_mz: np.ndarray,
    left_intensity: np.ndarray,
    right_mz: np.ndarray,
    right_intensity: np.ndarray,
    *,
    offset_da: float,
    tolerance_da: float,
) -> tuple[int, float]:
    """Greedy one-to-one match after shifting the left spectrum by offset."""

    i = j = 0
    count = 0
    dot = 0.0
    while i < left_mz.size and j < right_mz.size:
        delta = float(right_mz[j] - (left_mz[i] + offset_da))
        if abs(delta) <= tolerance_da:
            # Resolve a local nearest-neighbour ambiguity before committing.
            if j + 1 < right_mz.size:
                next_delta = float(right_mz[j + 1] - (left_mz[i] + offset_da))
                if abs(next_delta) < abs(delta) and abs(next_delta) <= tolerance_da:
                    j += 1
                    continue
            if i + 1 < left_mz.size:
                next_delta = float(right_mz[j] - (left_mz[i + 1] + offset_da))
                if abs(next_delta) < abs(delta) and abs(next_delta) <= tolerance_da:
                    i += 1
                    continue
            count += 1
            dot += float(left_intensity[i] * right_intensity[j])
            i += 1
            j += 1
        elif delta < -tolerance_da:
            j += 1
        else:
            i += 1
    return count, dot


def _relatedness(
    left: _IndexedSpectrum,
    right: _IndexedSpectrum,
    signed_delta_da: float,
    *,
    fragment_match_da: float,
) -> tuple[int, int, float]:
    """Score unchanged plus modification-shifted fragment evidence.

    The shifted component tests both the full precursor delta (typical singly
    charged product ion) and half the delta (a common doubly charged product-ion
    case). This is only a relatedness screen; it does not localize which fragment
    carries a modification.
    """

    unchanged_count, unchanged_dot = _matched_intensity_dot(
        left.mz,
        left.intensity,
        right.mz,
        right.intensity,
        offset_da=0.0,
        tolerance_da=fragment_match_da,
    )
    shifted = [
        _matched_intensity_dot(
            left.mz,
            left.intensity,
            right.mz,
            right.intensity,
            offset_da=offset,
            tolerance_da=fragment_match_da,
        )
        for offset in (signed_delta_da, signed_delta_da / 2.0)
    ]
    shifted_count, shifted_dot = max(shifted, key=lambda item: (item[1], item[0]))
    # The two terms generally represent disjoint peak families for the mass
    # shifts of interest. Cap at one because the source vectors are normalized.
    similarity = min(1.0, unchanged_dot + shifted_dot)
    return unchanged_count, unchanged_count + shifted_count, similarity


def _cluster_observations(
    observations: Iterable[_MassShiftObservation],
    tolerance_da: float,
    *,
    minimum_pairs: int,
    minimum_unique_spectra: int,
) -> list[_Cluster]:
    ordered = sorted(observations, key=lambda item: item.delta_mass_da)
    if not ordered:
        return []
    groups: list[list[_MassShiftObservation]] = []
    current: list[_MassShiftObservation] = []
    center = 0.0
    for observation in ordered:
        if not current:
            current = [observation]
            center = observation.delta_mass_da
            continue
        if abs(observation.delta_mass_da - center) <= tolerance_da:
            current.append(observation)
            center += (observation.delta_mass_da - center) / len(current)
        else:
            groups.append(current)
            current = [observation]
            center = observation.delta_mass_da
    groups.append(current)

    clusters: list[_Cluster] = []
    for group in groups:
        spectrum_ids = {item.left_id for item in group} | {item.right_id for item in group}
        if len(group) < minimum_pairs or len(spectrum_ids) < minimum_unique_spectra:
            continue
        values = np.asarray([item.delta_mass_da for item in group], dtype=float)
        center_da = float(np.median(values))
        mad = float(np.median(np.abs(values - center_da)))
        sigma = 1.4826 * mad
        clusters.append(
            _Cluster(
                center_da=center_da,
                sigma_da=sigma,
                minimum_da=float(values.min()),
                maximum_da=float(values.max()),
                pair_support=len(group),
                unique_spectrum_support=len(spectrum_ids),
                median_similarity=float(median(item.similarity for item in group)),
                mean_neutral_mass=float(median(item.mean_neutral_mass for item in group)),
            )
        )
    return clusters


def _artifact_classification(
    mass_da: float,
    tolerance_da: float,
    isotope_mass_da: float,
) -> dict[str, Any] | None:
    for isotope_count in range(1, 5):
        theoretical = isotope_count * isotope_mass_da
        if abs(mass_da - theoretical) <= tolerance_da:
            return {
                "classification": "isotope-like",
                "name": f"C13 isotope-selection shift x{isotope_count}",
                "theoretical_delta_mass_da": theoretical,
                "residual_da": mass_da - theoretical,
            }
    for name, theoretical in _ADDUCT_SHIFTS_DA.items():
        if abs(mass_da - theoretical) <= tolerance_da:
            return {
                "classification": "adduct-like",
                "name": name,
                "theoretical_delta_mass_da": theoretical,
                "residual_da": mass_da - theoretical,
            }
    return None


def _matching_modifications(
    mass_da: float,
    tolerance_da: float,
    records: Iterable[ModificationRecord],
    *,
    maximum_candidates: int,
) -> list[dict[str, Any]]:
    matches: list[tuple[float, ModificationRecord]] = []
    for record in records:
        residual = mass_da - abs(record.delta_mass_da)
        if abs(residual) <= tolerance_da:
            matches.append((residual, record))
    matches.sort(key=lambda item: (abs(item[0]), item[1].accession, item[1].name))
    return [
        {
            "unimod_accession": record.accession,
            "name": record.name,
            "theoretical_delta_mass_da": record.delta_mass_da,
            "absolute_theoretical_delta_mass_da": abs(record.delta_mass_da),
            "residual_da": residual,
            "origins": list(record.origins),
            "term_specificities": list(record.term_specificities),
            "source_classification": record.source_classification or None,
        }
        for residual, record in matches[:maximum_candidates]
    ]


@dataclass(slots=True)
class MassShiftCollector:
    """Detect recurrent modification-like neutral precursor mass shifts.

    Candidate retrieval is an inverted index over coarse bins from the strongest
    centroid fragment peaks. Only spectra sharing several bins are scored. The
    scorer combines unchanged fragment matches with matches shifted by the
    precursor delta or half-delta. Accepted deltas are clustered after the run.

    When a frozen v19 mass-error collector is supplied, its precursor precision
    is used read-only to widen the delta clustering/UniMod window as needed. The
    mass-error estimator itself is never changed by this collector.
    """

    precision_source: RepeatSpectrumMassErrorCollector | None = None
    modifications: tuple[ModificationRecord, ...] | None = None
    oms: Any | None = None
    top_peaks: int = 60
    fragment_bin_da: float = 1.0
    fragment_match_da: float = 0.05
    precursor_exclusion_da: float = 2.0
    minimum_fragment_peaks: int = 12
    minimum_shared_bins: int = 5
    minimum_unchanged_matches: int = 4
    minimum_total_matches: int = 10
    minimum_similarity: float = 0.25
    minimum_shift_da: float = 0.5
    maximum_shift_da: float = 500.0
    candidate_rt_window_seconds: float = 900.0
    maximum_candidates_per_spectrum: int = 24
    maximum_pairs_per_spectrum: int = 2
    maximum_indexed_spectra: int = 20_000
    maximum_postings_per_bin: int = 64
    maximum_accepted_pairs: int = 100_000
    minimum_cluster_pairs: int = 8
    minimum_cluster_unique_spectra: int = 6
    base_cluster_tolerance_da: float = 0.01
    base_unimod_tolerance_da: float = 0.02
    precision_sigma_multiplier: float = 6.0
    maximum_modification_candidates: int = 12
    _spectra: dict[int, _IndexedSpectrum] = field(default_factory=dict)
    _active: deque[int] = field(default_factory=deque)
    _postings: defaultdict[int, deque[int]] = field(
        default_factory=lambda: defaultdict(deque)
    )
    _observations: list[_MassShiftObservation] = field(default_factory=list)
    _next_identifier: int = 0
    total_ms2: int = 0
    eligible_ms2: int = 0
    scored_pairs: int = 0
    accepted_pairs: int = 0
    index_evictions: int = 0
    excluded_profile_or_unknown: int = 0
    excluded_negative_polarity: int = 0
    excluded_missing_precursor: int = 0
    excluded_low_fragment_peaks: int = 0
    accepted_pair_cap_hits: int = 0
    _catalog_error: str | None = None
    proton_mass_u: float = field(init=False)
    isotope_mass_u: float = field(init=False)

    def __post_init__(self) -> None:
        for name, value in (
            ("fragment_bin_da", self.fragment_bin_da),
            ("fragment_match_da", self.fragment_match_da),
            ("precursor_exclusion_da", self.precursor_exclusion_da),
            ("minimum_shift_da", self.minimum_shift_da),
            ("maximum_shift_da", self.maximum_shift_da),
            ("candidate_rt_window_seconds", self.candidate_rt_window_seconds),
            ("base_cluster_tolerance_da", self.base_cluster_tolerance_da),
            ("base_unimod_tolerance_da", self.base_unimod_tolerance_da),
            ("precision_sigma_multiplier", self.precision_sigma_multiplier),
        ):
            if not math.isfinite(value) or value <= 0:
                raise ValueError(f"{name} must be finite and positive")
        if self.maximum_shift_da <= self.minimum_shift_da:
            raise ValueError("maximum_shift_da must be greater than minimum_shift_da")
        for name, value in (
            ("top_peaks", self.top_peaks),
            ("minimum_fragment_peaks", self.minimum_fragment_peaks),
            ("minimum_shared_bins", self.minimum_shared_bins),
            ("minimum_unchanged_matches", self.minimum_unchanged_matches),
            ("minimum_total_matches", self.minimum_total_matches),
            ("maximum_candidates_per_spectrum", self.maximum_candidates_per_spectrum),
            ("maximum_pairs_per_spectrum", self.maximum_pairs_per_spectrum),
            ("maximum_indexed_spectra", self.maximum_indexed_spectra),
            ("maximum_postings_per_bin", self.maximum_postings_per_bin),
            ("maximum_accepted_pairs", self.maximum_accepted_pairs),
            ("minimum_cluster_pairs", self.minimum_cluster_pairs),
            ("minimum_cluster_unique_spectra", self.minimum_cluster_unique_spectra),
            ("maximum_modification_candidates", self.maximum_modification_candidates),
        ):
            if value < 1:
                raise ValueError(f"{name} must be positive")
        if not 0 < self.minimum_similarity <= 1:
            raise ValueError("minimum_similarity must be in (0, 1]")
        self.proton_mass_u, self.isotope_mass_u = _openms_constants(self.oms)

    def _expire_if_needed(self) -> None:
        while len(self._active) >= self.maximum_indexed_spectra:
            identifier = self._active.popleft()
            if self._spectra.pop(identifier, None) is not None:
                self.index_evictions += 1

    def _bins(self, mz: np.ndarray) -> tuple[int, ...]:
        return tuple(sorted({int(math.floor(value / self.fragment_bin_da)) for value in mz}))

    def _candidate_ids(self, bins: tuple[int, ...]) -> list[tuple[int, int]]:
        counts: Counter[int] = Counter()
        for bin_id in bins:
            posting = self._postings[bin_id]
            while posting and posting[0] not in self._spectra:
                posting.popleft()
            for identifier in posting:
                if identifier in self._spectra:
                    counts[identifier] += 1
        return sorted(
            (
                (identifier, shared)
                for identifier, shared in counts.items()
                if shared >= self.minimum_shared_bins
            ),
            key=lambda item: (-item[1], item[0]),
        )[: self.maximum_candidates_per_spectrum]

    def _index(self, item: _IndexedSpectrum) -> None:
        self._expire_if_needed()
        self._spectra[item.identifier] = item
        self._active.append(item.identifier)
        for bin_id in item.bins:
            posting = self._postings[bin_id]
            posting.append(item.identifier)
            while len(posting) > self.maximum_postings_per_bin:
                posting.popleft()

    def consume_spectrum(self, spectrum: Spectrum) -> None:
        if spectrum.ms_level != 2:
            return
        self.total_ms2 += 1
        representation = spectrum.representation
        if representation == "unknown":
            representation = spectrum.estimated_representation or "unknown"
        if representation != "centroid":
            self.excluded_profile_or_unknown += 1
            return
        if spectrum.polarity == "negative":
            self.excluded_negative_polarity += 1
            return
        precursor = spectrum.precursors[0] if spectrum.precursors else None
        if (
            precursor is None
            or precursor.charge <= 0
            or not math.isfinite(precursor.mz)
            or precursor.mz <= 0
        ):
            self.excluded_missing_precursor += 1
            return
        mz, intensity = _top_peaks(
            spectrum,
            self.top_peaks,
            precursor_exclusion_da=self.precursor_exclusion_da,
        )
        if mz.size < self.minimum_fragment_peaks:
            self.excluded_low_fragment_peaks += 1
            return
        self.eligible_ms2 += 1
        identifier = self._next_identifier
        self._next_identifier += 1
        neutral_mass = precursor.mz * precursor.charge - precursor.charge * self.proton_mass_u
        current = _IndexedSpectrum(
            identifier=identifier,
            rt=float(spectrum.rt),
            precursor_mz=float(precursor.mz),
            charge=int(precursor.charge),
            neutral_mass=float(neutral_mass),
            mz=mz,
            intensity=intensity,
            bins=self._bins(mz),
        )

        accepted: list[tuple[float, _MassShiftObservation]] = []
        for candidate_id, _shared in self._candidate_ids(current.bins):
            candidate = self._spectra.get(candidate_id)
            if candidate is None:
                continue
            if abs(current.rt - candidate.rt) > self.candidate_rt_window_seconds:
                continue
            signed_delta = current.neutral_mass - candidate.neutral_mass
            absolute_delta = abs(signed_delta)
            if not self.minimum_shift_da <= absolute_delta <= self.maximum_shift_da:
                continue
            self.scored_pairs += 1
            unchanged, total, similarity = _relatedness(
                candidate,
                current,
                signed_delta,
                fragment_match_da=self.fragment_match_da,
            )
            if (
                unchanged < self.minimum_unchanged_matches
                or total < self.minimum_total_matches
                or similarity < self.minimum_similarity
            ):
                continue
            observation = _MassShiftObservation(
                delta_mass_da=float(absolute_delta),
                left_id=candidate.identifier,
                right_id=current.identifier,
                mean_neutral_mass=(candidate.neutral_mass + current.neutral_mass) / 2.0,
                similarity=similarity,
            )
            accepted.append((similarity, observation))

        for _score, observation in sorted(
            accepted,
            key=lambda item: (-item[0], item[1].left_id),
        )[: self.maximum_pairs_per_spectrum]:
            if len(self._observations) >= self.maximum_accepted_pairs:
                self.accepted_pair_cap_hits += 1
                break
            self._observations.append(observation)
            self.accepted_pairs += 1

        self._index(current)

    def _effective_tolerances(self) -> tuple[float, float, dict[str, Any]]:
        calibration: dict[str, Any] = {
            "source": "fixed-fallback",
            "cluster_tolerance_da": self.base_cluster_tolerance_da,
            "unimod_tolerance_da": self.base_unimod_tolerance_da,
        }
        cluster_tolerance = self.base_cluster_tolerance_da
        if self.precision_source is not None:
            precision = self.precision_source.precursor_precision_ppm()
            if (
                precision is not None
                and self.precision_source.precursor_paired_spectra
                >= self.precision_source.min_tolerance_pairs
                and self.precision_source.precursor_clusters_used
                >= self.precision_source.min_tolerance_clusters
                and self._observations
            ):
                sigma_ppm = float(precision["single_measurement_sigma"])
                representative_mass = float(
                    median(item.mean_neutral_mass for item in self._observations)
                )
                pair_sigma_da = (
                    math.sqrt(2.0) * sigma_ppm * representative_mass / 1_000_000.0
                )
                calibrated = self.precision_sigma_multiplier * pair_sigma_da
                cluster_tolerance = max(cluster_tolerance, calibrated)
                calibration = {
                    "source": "v19-repeat-precursor-precision-read-only",
                    "single_measurement_sigma_ppm": sigma_ppm,
                    "representative_neutral_mass_da": representative_mass,
                    "pair_sigma_da": pair_sigma_da,
                    "sigma_multiplier": self.precision_sigma_multiplier,
                    "cluster_tolerance_da": cluster_tolerance,
                    "unimod_tolerance_da": max(
                        self.base_unimod_tolerance_da, cluster_tolerance
                    ),
                }
        return (
            cluster_tolerance,
            max(self.base_unimod_tolerance_da, cluster_tolerance),
            calibration,
        )

    def _modification_records(self) -> tuple[ModificationRecord, ...]:
        if self.modifications is not None:
            return self.modifications
        try:
            self.modifications = load_openms_modifications(self.oms)
        except (ImportError, AttributeError, RuntimeError, TypeError, ValueError) as exc:
            self._catalog_error = f"{type(exc).__name__}: {exc}"
            self.modifications = ()
        return self.modifications

    def annotations(self) -> list[Annotation]:
        cluster_tolerance, unimod_tolerance, calibration = self._effective_tolerances()
        clusters = _cluster_observations(
            self._observations,
            cluster_tolerance,
            minimum_pairs=self.minimum_cluster_pairs,
            minimum_unique_spectra=self.minimum_cluster_unique_spectra,
        )
        records = self._modification_records()
        payload: list[dict[str, Any]] = []
        for cluster in clusters:
            artifact = _artifact_classification(
                cluster.center_da,
                unimod_tolerance,
                self.isotope_mass_u,
            )
            candidates = [] if artifact else _matching_modifications(
                cluster.center_da,
                unimod_tolerance,
                records,
                maximum_candidates=self.maximum_modification_candidates,
            )
            if artifact:
                classification = artifact["classification"]
            elif candidates:
                classifications = {
                    str(item.get("source_classification") or "").casefold()
                    for item in candidates
                }
                classification = (
                    "putative-ptm"
                    if any("post" in item for item in classifications)
                    else "putative-modification"
                )
            else:
                classification = "unknown"
            if cluster.pair_support >= 50 and cluster.unique_spectrum_support >= 20:
                confidence = "high-support"
            elif cluster.pair_support >= 15 and cluster.unique_spectrum_support >= 10:
                confidence = "moderate-support"
            else:
                confidence = "preliminary-support"
            payload.append(
                {
                    "delta_mass_da": cluster.center_da,
                    "cluster_sigma_da": cluster.sigma_da,
                    "cluster_min_da": cluster.minimum_da,
                    "cluster_max_da": cluster.maximum_da,
                    "pair_support": cluster.pair_support,
                    "unique_spectrum_support": cluster.unique_spectrum_support,
                    "median_spectral_similarity": cluster.median_similarity,
                    "classification": classification,
                    "confidence": confidence,
                    "match_tolerance_da": unimod_tolerance,
                    "artifact_candidate": artifact,
                    "unimod_candidates": candidates,
                    "orientation": (
                        "absolute neutral-mass difference; heavier/lighter modification "
                        "direction unresolved"
                    ),
                }
            )
        payload.sort(
            key=lambda item: (-int(item["pair_support"]), float(item["delta_mass_da"]))
        )
        method = (
            "identification-free recurrent mass-shift scout v22.0; centroid MS2; bounded "
            f"top-{self.top_peaks} fragment index; >= {self.minimum_shared_bins} shared coarse "
            "fragment bins; unchanged plus precursor-delta/half-delta fragment matching; "
            "absolute charge-aware neutral precursor deltas; OpenMS/UniMod annotation"
        )
        detail = (
            "Candidate chemistry only. No peptide database search or sequence/site localization is "
            "performed. Isotope/adduct-like shifts are classified before UniMod annotation; all "
            "mass-compatible UniMod candidates are preserved up to the configured cap, and "
            "unmatched recurrent shifts remain visible. The pair-based scout preferentially detects "
            "variable/co-occurring modified and reference-like forms; a modification present on every "
            "corresponding peptide without a related reference form may be invisible."
        )
        diagnostics = {
            "total_ms2": self.total_ms2,
            "eligible_centroid_ms2": self.eligible_ms2,
            "indexed_spectra_current": len(self._spectra),
            "index_evictions": self.index_evictions,
            "scored_candidate_pairs": self.scored_pairs,
            "accepted_related_pairs": self.accepted_pairs,
            "accepted_pair_cap_hits": self.accepted_pair_cap_hits,
            "recurrent_clusters": len(payload),
            "excluded_profile_or_unknown": self.excluded_profile_or_unknown,
            "excluded_negative_polarity": self.excluded_negative_polarity,
            "excluded_missing_precursor": self.excluded_missing_precursor,
            "excluded_low_fragment_peaks": self.excluded_low_fragment_peaks,
            "calibration": calibration,
            "openms_unimod_records": len(records),
            "openms_unimod_error": self._catalog_error,
            "proton_mass_u": self.proton_mass_u,
            "c13_c12_mass_difference_u": self.isotope_mass_u,
        }
        return [
            Annotation(
                "putative_modification_mass_shifts",
                payload or None,
                EvidenceKind.INFERRED if payload else EvidenceKind.UNAVAILABLE,
                method,
                detail,
                support=sum(int(item["pair_support"]) for item in payload),
                total=self.eligible_ms2,
            ),
            Annotation(
                "mass_shift_scout_diagnostics",
                diagnostics,
                EvidenceKind.INFERRED if self.eligible_ms2 else EvidenceKind.UNAVAILABLE,
                method,
                "Runtime/support accounting only; never an SDRF annotation.",
                support=self.accepted_pairs,
                total=self.eligible_ms2,
            ),
        ]
