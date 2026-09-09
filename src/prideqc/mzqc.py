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
from prideqc.models import AnalysisResult, CVTerm, Metric

SECOND = CVTerm("UO:0000010", "second")
COUNT = CVTerm("UO:0000189", "count unit")
HZ = CVTerm("UO:0000106", "hertz")
FRACTION = CVTerm("UO:0000191", "fraction")


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
}

DEFINITIONS = {
    "NumberOfMSLevels": "Number of distinct observed MS levels.",
    "NumberOfSpectra": "Total spectrum count over all levels.",
    "NumberOfSpectralPeaks": "Total valid stored peak data points over all levels; not chromatographic features.",
    "NumberOfChromatograms": "Number of stored chromatograms.",
    "NumberOfChromatogramDataPoints": "Total stored chromatogram array lengths, not chromatographic peaks.",
    "ChromatogramTypes": "Counts grouped by OpenMS ChromatogramType enum value; engine version is recorded.",
    "ChromatogramRTRange": "Minimum and maximum finite chromatogram time, in seconds.",
    "BasePeak_All_Max": (
        "Largest valid intensity across all spectra, including MS3 and higher; arbitrary "
        "intensity units."
    ),
    "FAIMS_CV_Values": "Distinct finite FAIMS compensation voltages, sorted, in volts; not drift times.",
    "FAIMS_CV_Count": "Number of distinct finite FAIMS compensation voltages.",
    "FAIMS_CV_Range": "Minimum and maximum finite FAIMS compensation voltages, in volts.",
    "MS1_to_MS2_Ratio": "MS1 scan count divided by MS2 scan count; undefined when denominator is zero.",
    "AvgCycleTime_MS1": "Mean strictly positive adjacent sorted MS1 RT differences, in seconds.",
    "MedianTIC_in_RT_MS1_IQR": (
        "Median TIC in middle two scan-index quartile groups (recycle-and-sort assignment), "
        "MS1. ID-free proxy; not an identified-peptide metric."
    ),
    "TIC_MS1_MedianInHalfRange": (
        "Median TIC in the shortest window containing ceil(n/2) finite-RT MS1 scans; first "
        "window wins ties. ID-free proxy."
    ),
    "ExtentPrecursorIntensity_95over5_MS2": (
        "95th/5th percentile ratio of positive recorded first-precursor MS2 intensities; "
        "ID-free proxy, no fragment-TIC fallback."
    ),
    "MS2_PrecursorCharge_Fractions": (
        "Counts and fractions of first-precursor charges, denominator all MS2 scans; missing "
        "or nonpositive charge is unknown."
    ),
    "ChargeMean": "Mean positive first-precursor MS2 charge; unknown charges excluded.",
    "ChargeMedian": "Median positive first-precursor MS2 charge; unknown charges excluded.",
    "ChargeMin": "Minimum positive first-precursor MS2 charge.",
    "ChargeMax": "Maximum positive first-precursor MS2 charge.",
    "ChargeRatio_3over2": (
        "Count of MS2 charge 3 divided by charge 2; unavailable unless both are present (rawQC "
        "convention)."
    ),
    "ChargeRatio_4over2": (
        "Count of MS2 charge 4 divided by charge 2; unavailable unless both are present (rawQC "
        "convention)."
    ),
}

LEVEL_DEFINITIONS = {
    "NumberOfSpectra_MS{n}": "Count of spectra at MS level {n}, including empty scans.",
    "RtRange_MS{n}": "Minimum and maximum finite RT in seconds, MS{n}.",
    "EmptyScans_MS{n}": "MS{n} scans with no valid positive intensity (zero TIC).",
    "PeakDensity_MS{n}_Quantiles": "25th, 50th, 75th percentiles of valid stored peak counts, including zero, MS{n}.",
    "BasePeak_MS{n}_Mean": (
        "Mean of per-spectrum maximum valid intensity, excluding scans with no valid peaks, "
        "MS{n}."
    ),
    "TIC_MS{n}_Median": "Median sum of valid peak intensities per scan, MS{n}.",
    "TIC_MS{n}_CV": "Sample standard deviation (ddof=1) divided by mean TIC, MS{n}.",
    "TIC_MS{n}_Area": (
        "Trapezoidal TIC integral against RT, stable time order, intensity times seconds, "
        "MS{n}."
    ),
    "TIC_MS{n}_Area_RTQuantiles": (
        "Four exact integrals of piecewise-linear TIC split at scan-RT quartiles, intensity "
        "times seconds, MS{n}."
    ),
    "TIC_MS{n}_QuartileLogRatios": "Natural logs of TIC Q2/Q1, Q3/Q2, Q4/Q3 (Q4=max), MS{n}.",
    "TIC_MS{n}_ChangeQuartileLogRatios": "Natural logs of absolute adjacent TIC change Q2/Q1, Q3/Q2, Q4/Q3, MS{n}.",
    "ScanRate_MS{n}": "Count of finite-RT scans divided by that level's RT span, in scans per minute, MS{n}.",
    "FastestFrequency_MS{n}": "Maximum scan count in a closed 60-second window divided by 60, in Hz, MS{n}.",
    "RT_MS{n}_Quantiles": "Four scan-RT quartile interval fractions of that level's RT span, MS{n}.",
    "RT_MS{n}_IQR": "75th minus 25th percentile of finite retention times, in seconds, MS{n}.",
    "RT_MS{n}_IQRRate": (
        "Count of finite-RT scans in the closed RT interquartile interval divided by its width "
        "in seconds, MS{n}; ID-free proxy."
    ),
    "ObservedMzRange_MS{n}": (
        "Min/max valid peak m/z; this is an observed range, not the instrument scan window, "
        "MS{n}."
    ),
    "ScanWindow_MS{n}": "Envelope of recorded scan-window bounds in Th (m/z), MS{n}.",
    "PeakTypes_MS{n}": "Counts of centroid, profile and unknown scan representations, MS{n}.",
    "EstimatedPeakTypes_MS{n}": (
        "Optional OpenMS PeakTypeEstimator counts for scans with more than ten peaks, MS{n}; "
        "empty if estimation was not requested or eligible."
    ),
    "Polarity_MS{n}": "Counts by recorded polarity, including unknown, MS{n}.",
    "InvalidPeakCount_MS{n}": (
        "Count of excluded peaks with non-finite/nonpositive m/z or non-finite/negative "
        "intensity, MS{n}."
    ),
    "TIC_MS{n}_SignalJump10x_Count": "Adjacent sorted-TIC ratios >=10, with strictly positive previous TIC, MS{n}.",
    "TIC_MS{n}_SignalFall10x_Count": (
        "Adjacent sorted-TIC ratios <=0.1, with strictly positive previous TIC; zero next TIC "
        "counts, MS{n}."
    ),
    "MzRange_MS{n}": "Range of positive finite first-precursor m/z, not fragment peak m/z, MS{n}.",
    "PrecursorMz_MS{n}_Median": "Median positive finite first-precursor m/z, MS{n}.",
    "PrecursorIntensity_MS{n}_Quantiles": (
        "Quartiles of recorded positive first-precursor intensity, MS{n}; no fragment-TIC "
        "substitution."
    ),
    "PrecursorIntensity_MS{n}_Mean": "Mean recorded positive first-precursor intensity, MS{n}.",
    "PrecursorIntensity_MS{n}_Sd": "Sample SD of recorded positive first-precursor intensity (ddof=1), MS{n}.",
    "PrecursorIntensity_MS{n}_MissingCount": (
        "Scans lacking a positive finite first-precursor intensity, including the OpenMS unset "
        "default zero, MS{n}."
    ),
    "IsolationWidth_MS{n}_Median": "Median positive first-precursor lower+upper isolation offsets in Th (m/z), MS{n}.",
    "MultiplePrecursors_MS{n}_Count": (
        "Scans with multiple precursors; first precursor supplies precursor-distribution "
        "metrics, MS{n}."
    ),
}


def _local_accession(key: str) -> str:
    return "QCPRIDE:" + hashlib.sha256(key.encode()).hexdigest()[:20].upper()


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
    if metric.key in {"ChromatogramTypes", "MS2_PrecursorCharge_Fractions"} or metric.key.startswith(
        ("PeakTypes_", "EstimatedPeakTypes_", "Polarity_"),
    ):
        shape = "table"
    elif metric.key == "FAIMS_CV_Values" or any(part in metric.key for part in ("Range", "Quantiles", "LogRatios", "ScanWindow")):
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
            item = {"accession": term.accession, "name": term.name,
                    "description": f"{term.description} Internal key: {metric.key}.", "value": value}
            if term.unit:
                item["unit"] = {"accession": term.unit.accession, "name": term.unit.name}
            metrics.append(item)
        properties: list[dict[str, Any]] = [
            {"accession": term.accession, "name": term.name}
            for term in result.metadata.instruments
        ]
        if result.source_path:
            properties.append({"accession": "QCPRIDE:SOURCEFILE", "name": "original vendor file",
                               "value": str(result.source_path)})
        input_file: dict[str, Any] = {
            "name": result.input_path.name, "location": result.input_path.resolve().as_uri(),
            "fileFormat": {"accession": "MS:1000584", "name": "mzML format"},
        }
        if properties:
            input_file["fileProperties"] = properties
        return {"mzQC": {
            "version": "1.0.0", "creationDate": datetime.now(UTC).isoformat(),
            "description": (
                "prideqc ID-free analysis; missing values and annotation evidence are in the "
                "companion summary."
            ),
            "runQualities": [{"metadata": {
                "label": result.input_path.name, "inputFiles": [input_file],
                "analysisSoftware": [
                    {"accession": "QCPRIDE:SOFTWARE", "name": "prideqc", "version": __version__},
                    {"accession": "MS:1000531", "name": "software", "version": result.engine_version,
                     "value": "pyOpenMS", "description": "pyOpenMS bindings used for mzML decoding", "uri": "https://openms.de"},
                ],
            }, "qualityMetrics": metrics}],
            "controlledVocabularies": [
                {"name": "PSI-MS", "uri": "https://raw.githubusercontent.com/HUPO-PSI/psi-ms-CV/master/psi-ms.obo"},
                {"name": "Unit Ontology", "uri": "http://purl.obolibrary.org/obo/uo.obo"},
                {"name": "prideqc local vocabulary", "version": "1", "uri": local_cv.resolve().as_uri()},
            ],
        }}

    def write(self, result: AnalysisResult, path: Path) -> None:
        cv_path = path.with_suffix(".obo")
        self.write_vocabulary(result.metrics, cv_path)
        write_json(path, self.build(result, cv_path))

    def write_vocabulary(self, metrics: list[Metric], path: Path) -> None:
        terms = [definition(metric) for metric in metrics]
        with atomic_text(path) as handle:
            handle.write("format-version: 1.2\ndata-version: 1\nontology: prideqc\n")
            for accession, name in (("QCPRIDE:SOFTWARE", "prideqc"), ("QCPRIDE:SOURCEFILE", "original vendor file")):
                handle.write(f"\n[Term]\nid: {accession}\nname: {name}\n")
            for term in terms:
                if not term.accession.startswith("QCPRIDE:"):
                    continue
                parent = {"scalar": "MS:4000003", "tuple": "MS:4000004", "table": "MS:4000005"}[term.shape]
                handle.write(
                    f"\n[Term]\nid: {term.accession}\nname: {term.name}\ndef: {json.dumps(
                        term.description,
                    )} [QCPRIDE:SOFTWARE]\nis_a: {parent}\n",
                )
