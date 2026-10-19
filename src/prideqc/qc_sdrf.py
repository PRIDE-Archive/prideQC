"""Standards-pending, read-only SDRF proposals for measured per-run QC.

The proposed headers are NOT recognised by the SDRF-Proteomics v1.1.0 term
registry. This module never mutates an SDRF. Once upstream standardisation is
complete, ``metric_annotations`` provides typed input for SDRFDocument.annotate.
"""

from __future__ import annotations

import csv
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Literal

from prideqc import __version__
from prideqc.io import atomic_text, write_json
from prideqc.models import AnalysisResult, Annotation, EvidenceKind
from prideqc.sdrf import MISSING, SDRFDocument, _equal, file_name


@dataclass(frozen=True, slots=True)
class QCSDRFMetric:
    """Proposed mapping only; not an approved SDRF column registration."""

    key: str
    accession: str
    label: str
    proposed_header: str
    unit_accession: str
    unit: str
    value_type: Literal["nonnegative_integer", "nonnegative_duration_seconds"]


QC_SDRF_PROPOSALS = (
    QCSDRFMetric(
        "NumberOfSpectra_MS1",
        "MS:4000059",
        "number of MS1 spectra",
        "comment[number of MS1 spectra]",
        "UO:0000189",
        "count",
        "nonnegative_integer",
    ),
    QCSDRFMetric(
        "NumberOfSpectra_MS2",
        "MS:4000060",
        "number of MS2 spectra",
        "comment[number of MS2 spectra]",
        "UO:0000189",
        "count",
        "nonnegative_integer",
    ),
    QCSDRFMetric(
        "ChromatographyDuration",
        "MS:4000067",
        "MS run duration",
        "comment[MS run duration]",
        "UO:0000010",
        "second",
        "nonnegative_duration_seconds",
    ),
)


@dataclass(frozen=True, slots=True)
class QCSDRFProposal:
    row: int | None
    data_file: str
    analyzed_file: str
    metric_key: str
    psi_ms_accession: str
    proposed_column: str
    value: str
    unit_accession: str
    unit: str
    status: str
    existing: str = ""
    reason: str = ""
    method: str = "prideqc.QCMetricCalculator.calculate"


def _format_value(spec: QCSDRFMetric, value: Any) -> str | None:
    """Strict measurement grammar, deliberately different from cohort grouping.

    Zero is valid for all three metrics; NaN/inf, booleans and negative values
    are never coerced to zero. The duration's ' s' suffix is *proposed* syntax.
    """
    if spec.value_type == "nonnegative_integer":
        if isinstance(value, bool) or not isinstance(value, int) or value < 0:
            return None
        return str(value)
    if (
        not isinstance(value, (int, float))
        or isinstance(value, bool)
        or not math.isfinite(value)
        or value < 0
    ):
        return None
    return f"{format(value, '.15g')} s"


def metric_annotations(result: AnalysisResult) -> list[Annotation]:
    """Build candidate annotations without modifying the result or its SDRF."""
    values: dict[str, list[Any]] = {}
    for metric in result.metrics:
        values.setdefault(metric.key, []).append(metric.value)
    output: list[Annotation] = []
    for spec in QC_SDRF_PROPOSALS:
        candidates = values.get(spec.key, [])
        if len(candidates) != 1:
            continue
        text = _format_value(spec, candidates[0])
        if text is None:
            continue
        output.append(
            Annotation(
                field=f"qc_{spec.key}",
                value={
                    "value": candidates[0],
                    "psi_ms_accession": spec.accession,
                    "unit_accession": spec.unit_accession,
                    "unit": spec.unit,
                    "registration_status": "proposed_not_registered",
                },
                kind=EvidenceKind.OBSERVED,
                method="prideqc.QCMetricCalculator.calculate",
                detail="Per-file observed spectrum measurement; SDRF header pending upstream acceptance",
                sdrf_column=spec.proposed_header,
                sdrf_value=text,
            )
        )
    return output


def collect_qc_sdrf_proposals(
    results: list[AnalysisResult],
    document: SDRFDocument | None = None,
    aliases: dict[str, str] | None = None,
) -> list[QCSDRFProposal]:
    """Inspect original SDRF cells using the standard exact file resolver.

    Every matching SDRF row receives an independent proposal; repeated sample
    rows are intentionally retained. No cell is ever changed here.
    """
    targets: list[tuple[int | None, str, list[str] | None, AnalysisResult | None]]
    if document is not None:
        lookup = document._result_lookup(results, aliases)
        targets = [
            (
                row_number,
                row[document.file_column],
                row,
                lookup.get(file_name(row[document.file_column])),
            )
            for row_number, row in enumerate(document.rows, start=2)
        ]
    else:
        if aliases:
            raise ValueError("File-map aliases require an original SDRF")
        targets = [(None, result.input_path.name, None, result) for result in results]

    proposals: list[QCSDRFProposal] = []
    for row_number, data_file, original, result in targets:
        values: dict[str, list[Any]] = {}
        if result is not None:
            for metric in result.metrics:
                values.setdefault(metric.key, []).append(metric.value)
        for spec in QC_SDRF_PROPOSALS:
            status, proposed, previous, reason = "", "", "", ""
            if result is None:
                status, reason = "unmatched", "No analyzed result matches this SDRF data file"
            else:
                candidates = values.get(spec.key, [])
                if not candidates:
                    status, reason = "missing_metric", "Source summary does not contain this metric"
                elif len(candidates) != 1:
                    status, reason = (
                        "ambiguous_metric",
                        "Source summary contains duplicate metric keys",
                    )
                elif candidates[0] is None:
                    status, reason = "unavailable_metric", "No measured value available"
                else:
                    proposed = _format_value(spec, candidates[0]) or ""
                    if not proposed:
                        status, reason = (
                            "invalid_metric",
                            "Expected nonnegative finite measurement of the declared type",
                        )
                    else:
                        status = "proposed_unregistered"
                        if original is not None and document is not None:
                            indices = document.indices(spec.proposed_header)
                            if len(indices) > 1:
                                status, reason = (
                                    "ambiguous_column",
                                    "Duplicate proposed-header columns",
                                )
                            elif indices:
                                previous = original[indices[0]]
                                if _equal(previous, proposed):
                                    status = "existing_match_unregistered"
                                elif previous.strip().casefold() in MISSING:
                                    status = "existing_missing_unregistered"
                                else:
                                    status, reason = (
                                        "existing_conflict_unregistered",
                                        "Existing value preserved",
                                    )
            proposals.append(
                QCSDRFProposal(
                    row=row_number,
                    data_file=data_file,
                    analyzed_file=result.input_path.name if result is not None else "",
                    metric_key=spec.key,
                    psi_ms_accession=spec.accession,
                    proposed_column=spec.proposed_header,
                    value=proposed,
                    unit_accession=spec.unit_accession,
                    unit=spec.unit,
                    status=status,
                    existing=previous,
                    reason=reason,
                )
            )
    return proposals


def write_qc_sdrf_proposals(
    output: Path,
    results: list[AnalysisResult],
    document: SDRFDocument | None = None,
    aliases: dict[str, str] | None = None,
) -> list[QCSDRFProposal]:
    """Write a TSV and provenance-rich JSON sidecar, never SDRF columns."""
    proposals = collect_qc_sdrf_proposals(results, document, aliases)
    fields = list(QCSDRFProposal.__dataclass_fields__)
    with atomic_text(output / "sdrf-qc-proposals.tsv") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(asdict(record) for record in proposals)
    write_json(
        output / "sdrf-qc-proposals.json",
        {
            "schema_version": 1,
            "prideqc_version": __version__,
            "status": "standards_pending_not_written_to_sdrf",
            "mapping": [asdict(spec) for spec in QC_SDRF_PROPOSALS],
            "proposals": [asdict(record) for record in proposals],
        },
    )
    return proposals
