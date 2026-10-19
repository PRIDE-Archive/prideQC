"""mzQC 1.0 serialization without a runtime mzQC/ontology dependency."""

from __future__ import annotations

import hashlib
import json
import re
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from prideqc import __version__
from prideqc.io import atomic_text, json_safe, write_json
from prideqc.models import AnalysisResult, Annotation, CVTerm, Metric

SECOND = CVTerm("UO:0000010", "second")
COUNT = CVTerm("UO:0000189", "count unit")
HZ = CVTerm("UO:0000106", "hertz")
FRACTION = CVTerm("UO:0000191", "fraction")
MZ = CVTerm("MS:1000040", "m/z")


@dataclass(frozen=True, slots=True)
class MetricDefinition:
    accession: str
    name: str
    description: str
    shape: str = "scalar"
    unit: CVTerm | None = None


# Only map terms with an unambiguous value contract. Other quantities have
# explicit local definitions, rather than borrowing a vaguely related accession.
STANDARD = {
    "DIAIsolationWindowCount": MetricDefinition(
        "MS:4000194", "DIA isolation window count",
        "Distinct DIA MS2 isolation windows by lower/upper m/z boundaries and recorded "
        "FAIMS or window ion-mobility parameters; requires complete window metadata.",
        unit=COUNT,
    ),
    "DIAIsolationWindowWidths": MetricDefinition(
        "MS:4000195", "DIA isolation window m/z widths",
        "Minimum and maximum widths of distinct DIA isolation windows, in m/z.",
        "tuple", MZ,
    ),
    "ChromatographyDuration": MetricDefinition(
        "MS:4000053",
        "chromatography duration",
        "Last minus first finite spectrum RT, over all MS levels.",
        unit=SECOND,
    ),
    "NumberOfSpectra_MS1": MetricDefinition(
        "MS:4000059",
        "number of MS1 spectra",
        "Count of all MS1 spectra, including empty scans.",
        unit=COUNT,
    ),
    "NumberOfSpectra_MS2": MetricDefinition(
        "MS:4000060",
        "number of MS2 spectra",
        "Count of all MS2 spectra, including empty scans.",
        unit=COUNT,
    ),
    "PeakDensity_MS1_Quantiles": MetricDefinition(
        "MS:4000061",
        "MS1 density quantiles",
        "25th, 50th, 75th percentiles of valid peaks per MS1 spectrum; empty scans included.",
        "tuple",
    ),
    "PeakDensity_MS2_Quantiles": MetricDefinition(
        "MS:4000062",
        "MS2 density quantiles",
        "25th, 50th, 75th percentiles of valid peaks per MS2 spectrum; empty scans included.",
        "tuple",
    ),
    "FastestFrequency_MS1": MetricDefinition(
        "MS:4000065",
        "fastest frequency for MS level 1 collection",
        "Maximum scan count in a closed 60-second window divided by 60, including short runs.",
        unit=HZ,
    ),
    "FastestFrequency_MS2": MetricDefinition(
        "MS:4000066",
        "fastest frequency for MS level 2 collection",
        "Maximum scan count in a closed 60-second window divided by 60, including short runs.",
        unit=HZ,
    ),
    "RT_MS1_Quantiles": MetricDefinition(
        "MS:4000184",
        "MS1 quantile RT fraction",
        "Four consecutive scan-RT quartile widths divided by the MS1 time span.",
        "tuple",
        FRACTION,
    ),
    "RT_MS2_Quantiles": MetricDefinition(
        "MS:4000185",
        "MS2 quantile RT fraction",
        "Four consecutive scan-RT quartile widths divided by the MS2 time span.",
        "tuple",
        FRACTION,
    ),
    "RT_TIC_Quantiles": MetricDefinition(
        "MS:4000183",
        "TIC quantile RT fraction",
        (
            "Four intervals at first scan-wise cumulative MS1 TIC quartile crossings, divided "
            "by the MS1 time span; this uses accumulation, not integration."
        ),
        "tuple",
        FRACTION,
    ),
    "ChargeMean": MetricDefinition(
        "MS:4000173",
        "mean MS2 precursor charge in all spectra",
        "Mean of recorded positive first-precursor MS2 charge states; unknown charges excluded.",
    ),
    "ChargeMedian": MetricDefinition(
        "MS:4000175",
        "median MS2 precursor charge in all spectra",
        "Median of recorded positive first-precursor MS2 charge states; unknown charges excluded.",
    ),
    "ChargeRatio_3over2": MetricDefinition(
        "MS:4000169",
        "ratio of 3+ over 2+ of all MS2 known precursor charges",
        "MS2 count at charge 3 divided by count at charge 2; unavailable if either charge is absent.",
    ),
    "ChargeRatio_4over2": MetricDefinition(
        "MS:4000171",
        "ratio of 4+ over 2+ of all MS2 known precursor charges",
        "MS2 count at charge 4 divided by count at charge 2; unavailable if either charge is absent.",
    ),
}

DEFINITIONS = {
    "NumberOfMSLevels": "Number of distinct observed MS levels.",
    "FAIMS_CV_Values": "Distinct finite FAIMS compensation voltages, sorted, in volts; not drift times.",
    "FAIMS_CV_Count": "Number of distinct finite FAIMS compensation voltages.",
    "FAIMS_CV_Range": "Minimum and maximum finite FAIMS compensation voltages, in volts.",
    "MS1_to_MS2_Ratio": "MS1 scan count divided by MS2 scan count; undefined when denominator is zero.",
    "AcquisitionCycle_Count": (
        "Count of MS1-delimited cycles containing at least one MS2 scan. MS2 scans before the "
        "first MS1 are ignored; a trailing cycle is included when it contains MS2 scans."
    ),
    "AcquisitionCycle_MS2Count_Quantiles": (
        "25th, 50th, 75th percentiles of MS2 scans per observed MS1-delimited cycle."
    ),
    "AcquisitionCycle_MS2Count_Mode": (
        "Most frequent MS2 scan count per observed MS1-delimited cycle; lower count wins ties."
    ),
    "AcquisitionCycle_TargetCoverageFraction": (
        "Fraction of MS2 scans inside observed MS1-delimited cycles with a finite positive first-"
        "precursor m/z paired with a positive isolation width."
    ),
    "AcquisitionCycle_TargetEligibleFraction": (
        "Fraction of observed MS1-delimited cycles eligible for target-set/order comparison."
    ),
    "AcquisitionCycle_UniqueTargetCount_Quantiles": (
        "25th, 50th, 75th percentiles of unique first-precursor isolation targets per eligible "
        "cycle after 0.1 Th quantization."
    ),
    "AvgCycleTime_MS1": "Mean strictly positive adjacent sorted MS1 RT differences, in seconds.",
}

LEVEL_DEFINITIONS = {
    "NumberOfSpectra_MS{n}": "Count of spectra at MS level {n}, including empty scans.",
    "PeakDensity_MS{n}_Quantiles": "25th, 50th, 75th percentiles of valid stored peak counts, including zero, MS{n}.",
    "FastestFrequency_MS{n}": "Maximum scan count in a closed 60-second window divided by 60, in Hz, MS{n}.",
    "RT_MS{n}_Quantiles": "Four scan-RT quartile interval fractions of that level's RT span, MS{n}.",
    "ObservedMzRange_MS{n}": (
        "Min/max valid peak m/z; this is an observed range, not the instrument scan window, MS{n}."
    ),

}


ANNOTATION_DESCRIPTIONS = {
    "ionization": (
        "ionization source CV evidence",
        "Recorded ion source PSI-MS CV terms from instrument configuration; "
        "not inferred from spectra or instrument model.",
    ),
    "mass_analyzers": (
        "mass analyzer CV evidence",
        "Recorded mass analyzer PSI-MS CV terms from instrument configuration; "
        "all recorded configurations are preserved without guessing.",
    ),
    "acquisition_method_conflict": (
        "acquisition method evidence conflict",
        "Recorded acquisition CV term disagrees with scan-pattern DDA/DIA inference. "
        "No automatic SDRF proposal can be made from this conflicting evidence.",
    ),
    "acquisition_method_evidence": (
        "acquisition method inference evidence",
        "Ontology-linked, conservative DDA/DIA inference from observed isolation windows "
        "and cycles. Support, denominator and classifier version retained. "
        "Not a verified instrument setting; no PRM/SRM/MRM inference from width alone.",
    ),
    "labeling_strategy_candidates": (
        "labelling strategy candidate evidence",
        "Exploratory TMT/iTRAQ reporter-ion evidence mapped to PRIDE family terms. "
        "Non-exclusive; does not establish SILAC, label-free status, plex or SDRF channels.",
    ),
    "estimated_precursor_mass_error_ppm": (
        "estimated precursor mass error precision",
        (
            "Robust repeat-observation estimate of single-measurement precursor mass precision "
            "in parts per million. This is measurement-precision evidence, not a recovered "
            "historical database-search tolerance."
        ),
    ),
    "estimated_precursor_mass_error_da": (
        "estimated precursor mass error precision in daltons",
        (
            "Robust repeat-observation estimate of single-measurement precursor mass precision "
            "in daltons. This is measurement-precision evidence, not a recovered historical "
            "database-search tolerance."
        ),
    ),
    "estimated_fragment_mass_error_ppm": (
        "estimated fragment mass error precision",
        (
            "Robust repeated-spectrum estimate of single-measurement fragment mass precision in "
            "parts per million, including the estimator's resolution-regime diagnostics."
        ),
    ),
    "estimated_fragment_mass_error_da": (
        "estimated fragment mass error precision in daltons",
        (
            "Robust repeated-spectrum estimate of single-measurement fragment mass precision in "
            "daltons, including the estimator's resolution-regime diagnostics."
        ),
    ),
    "suggested_precursor_search_tolerance_ppm": (
        "suggested precursor search tolerance",
        (
            "Precision-derived precursor search-tolerance starting point in parts per million. "
            "It is not a recovered historical search setting and is not guaranteed search-optimal."
        ),
    ),
    "suggested_fragment_search_tolerance_ppm": (
        "suggested fragment search tolerance",
        (
            "Precision-derived high-resolution fragment search-tolerance starting point in parts "
            "per million. It is emitted only when the estimator resolves a supported "
            "high-resolution regime."
        ),
    ),
    "suggested_fragment_search_tolerance_da": (
        "suggested fragment search tolerance in daltons",
        (
            "Precision-derived low-resolution fragment search-tolerance starting point in "
            "daltons. "
            "It is emitted only when the estimator resolves a supported uncensored "
            "low-resolution regime."
        ),
    ),
    "tmt_profile_reporter_evidence": (
        "TMT experimental profile reporter-ion evidence",
        (
            "Identification-free per-channel OpenMS PeakPickerHiRes profile MS2/MS3 "
            "reporter evidence, tracked separately from vendor centroid observations. "
            "Experimental only; no verified labeling kit, sample mapping or SDRF writeback."
        ),
    ),
    "tmt_reporter_evidence": (
        "TMT reporter-ion channel evidence",
        (
            "Identification-free, per-channel centroid MS2/MS3 reporter m/z matches and "
            "uncorrected signal totals using runtime pyOpenMS quantitation masses. "
            "Candidate evidence only; no labeling-kit or biological sample-channel assignment."
        ),
    ),
    "mass_error_estimator_diagnostics": (
        "mass error estimator diagnostics",
        (
            "Support and accounting diagnostics for the active repeat-observation mass-error "
            "backend, including paired-spectrum counts and resolution-regime evidence."
        ),
    ),
    "putative_modification_mass_shifts": (
        "putative modification mass shifts",
        (
            "Identification-free recurrent absolute neutral precursor-mass differences from "
            "fragment-related MS2 spectra after centroid normalization (including ephemeral "
            "OpenMS PeakPickerHiRes processing for profile MS2), with artifact-family "
            "classification and QC-filtered mass-compatible OpenMS/UniMod candidates. These are "
            "hypotheses, not localized PTM identifications."
        ),
    ),
    "mass_shift_scout_diagnostics": (
        "mass shift scout diagnostics",
        (
            "Runtime, support, profile-centroiding and reporting-accounting diagnostics for the "
            "identification-free recurrent mass-shift scout. Raw cluster counts are diagnostic; "
            "the bounded reported cluster set is the QC-facing result."
        ),
    ),
}


FILE_FORMATS = {
    ".raw": ("MS:1000563", "Thermo RAW format"),
    ".d": ("MS:1002817", "Bruker TDF format"),
    ".d.zip": ("MS:1002817", "Bruker TDF format"),
}


def file_format(path: Path) -> dict[str, str]:
    """Return the PSI-MS file-format term for an input or original source."""
    name = path.name.casefold()
    accession_name = next(
        (value for suffix, value in FILE_FORMATS.items() if name.endswith(suffix)),
        ("MS:1000584", "mzML format"),
    )
    return {"accession": accession_name[0], "name": accession_name[1]}


def _local_accession(key: str) -> str:
    return "QCPRIDE:" + hashlib.sha256(key.encode()).hexdigest()[:20].upper()


def annotation_definition(annotation: Annotation) -> MetricDefinition | None:
    """Map selected evidence annotations into structured mzQC metrics."""

    item = ANNOTATION_DESCRIPTIONS.get(annotation.field)
    if item is None:
        return None
    name, description = item
    return MetricDefinition(
        _local_accession(f"annotation:{annotation.field}"),
        name,
        description,
        "table",
    )


def definition(metric: Metric) -> MetricDefinition:
    if metric.key in STANDARD:
        return STANDARD[metric.key]
    description = DEFINITIONS.get(metric.key)
    if description is None:
        for pattern, text in LEVEL_DEFINITIONS.items():
            match = re.fullmatch(re.escape(pattern).replace(r"\{n\}", r"([1-9][0-9]*)"), metric.key)
            if match:
                description = text.format(n=match[1])
                break
    if description is None:
        raise ValueError(f"Unregistered QC metric: {metric.key}")
    # Shape is a property of the term, not of a particular observed value.
    # In particular, an unavailable tuple must not become a scalar CV term.
    if metric.key == "FAIMS_CV_Values" or any(
        part in metric.key for part in ("Range", "Quantiles", "LogRatios", "ScanWindow")
    ):
        shape = "tuple"
    else:
        shape = "scalar"
    return MetricDefinition(_local_accession(metric.key), metric.key, description, shape)


class MzQCWriter:
    """One runQuality per file, with provenance and a real companion local CV."""

    def build(self, result: AnalysisResult, local_cv: Path) -> dict[str, Any]:
        metrics = []
        for metric in result.metrics:
            term = definition(metric)
            value = json_safe(metric.value)
            if value is None:
                # Missing observations are explicit in summary JSON/TSV. Do not
                # turn missing values into falsely measured zero-valued metrics.
                continue
            item = {
                "accession": term.accession,
                "name": term.name,
                "description": f"{term.description} Internal key: {metric.key}.",
                "value": value,
            }
            if term.unit:
                item["unit"] = {"accession": term.unit.accession, "name": term.unit.name}
            metrics.append(item)
        for annotation in result.annotations:
            annotation_term = annotation_definition(annotation)
            value = json_safe(annotation.value)
            if annotation_term is None or value is None:
                continue
            metrics.append(
                {
                    "accession": annotation_term.accession,
                    "name": annotation_term.name,
                    "description": (
                        f"{annotation_term.description} Evidence field: {annotation.field}."
                    ),
                    "value": value,
                }
            )
        properties: list[dict[str, Any]] = [
            {"accession": term.accession, "name": term.name} for term in result.metadata.instruments
        ]
        if result.source_path:
            properties.append(
                {
                    "accession": "QCPRIDE:SOURCEFILE",
                    "name": "original vendor file",
                    "value": str(result.source_path),
                }
            )
        project_accession = (result.project_accession or "").strip().upper()
        if project_accession.startswith("PXD"):
            properties.append(
                {
                    "accession": "MS:1001919",
                    "name": "ProteomeXchange accession number",
                    "value": project_accession,
                }
            )
        provenance_path = result.source_path or result.input_path
        input_file: dict[str, Any] = {
            "name": provenance_path.name,
            "location": provenance_path.resolve().as_uri(),
            "fileFormat": file_format(provenance_path),
        }
        if properties:
            input_file["fileProperties"] = properties
        return {
            "mzQC": {
                "version": "1.0.0",
                "creationDate": datetime.now(UTC).isoformat(),
                "description": (
                    "prideqc ID-free analysis; missing values and annotation evidence are in the "
                    "companion summary."
                ),
                "runQualities": [
                    {
                        "metadata": {
                            "label": result.input_path.name,
                            "inputFiles": [input_file],
                            "analysisSoftware": [
                                {
                                    "accession": "QCPRIDE:SOFTWARE",
                                    "name": "prideqc",
                                    "version": __version__,
                                },
                                {
                                    "accession": "MS:1000531",
                                    "name": "software",
                                    "version": result.engine_version,
                                    "value": "pyOpenMS",
                                    "description": "pyOpenMS bindings used for mzML decoding",
                                    "uri": "https://openms.de",
                                },
                            ],
                        },
                        "qualityMetrics": metrics,
                    }
                ],
                "controlledVocabularies": [
                    {
                        "name": "PSI-MS",
                        "uri": "https://raw.githubusercontent.com/HUPO-PSI/psi-ms-CV/master/psi-ms.obo",
                    },
                    {"name": "Unit Ontology", "uri": "http://purl.obolibrary.org/obo/uo.obo"},
                    {
                        "name": "PRIDE",
                        "uri": (
                            "https://raw.githubusercontent.com/PRIDE-Archive/"
                            "pride-ontology/master/pride_cv.obo"
                        ),
                    },
                    {
                        "name": "prideqc local vocabulary",
                        "version": "1",
                        "uri": local_cv.resolve().as_uri(),
                    },
                ],
            }
        }

    def write(self, result: AnalysisResult, path: Path) -> None:
        cv_path = path.with_suffix(".obo")
        self.write_vocabulary(result.metrics, cv_path, result.annotations)
        write_json(path, self.build(result, cv_path))

    def write_vocabulary(
        self,
        metrics: list[Metric],
        path: Path,
        annotations: list[Annotation] | None = None,
    ) -> None:
        terms = [definition(metric) for metric in metrics]
        terms.extend(
            term
            for annotation in annotations or []
            if (term := annotation_definition(annotation)) is not None
            and annotation.value is not None
        )
        with atomic_text(path) as handle:
            handle.write("format-version: 1.2\ndata-version: 1\nontology: prideqc\n")
            for accession, name in (
                ("QCPRIDE:SOFTWARE", "prideqc"),
                ("QCPRIDE:SOURCEFILE", "original vendor file"),
            ):
                handle.write(f"\n[Term]\nid: {accession}\nname: {name}\n")
            for term in terms:
                if not term.accession.startswith("QCPRIDE:"):
                    continue
                parent = {"scalar": "MS:4000003", "tuple": "MS:4000004", "table": "MS:4000005"}[
                    term.shape
                ]
                description = json.dumps(term.description)
                handle.write(
                    f"\n[Term]\nid: {term.accession}\nname: {term.name}\ndef: {description} "
                    f"[QCPRIDE:SOFTWARE]\nis_a: {parent}\n",
                )
