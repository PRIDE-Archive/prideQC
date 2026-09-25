"""Lossless-column SDRF editing using the standard library's TSV support."""

from __future__ import annotations

import csv
import re
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path
from urllib.parse import unquote, urlsplit

from prideqc.io import atomic_text
from prideqc.models import AnalysisResult, Annotation, EvidenceKind

MISSING = frozenset({"", "not available", "not provided"})


def file_name(value: str) -> str:
    """Exact, case-sensitive basename; do not merge same-stem vendor files."""
    value = value.strip().replace("\\", "/")
    if "://" in value:
        value = unquote(urlsplit(value).path)
    return value.rsplit("/", 1)[-1]


def _equal(old: str, new: str) -> bool:
    # Compare complete CV accession sets before comparing normalized text.
    old_accessions = re.findall(r"(?:^|;)\s*AC=([^;]+)", old, flags=re.I)
    new_accessions = re.findall(r"(?:^|;)\s*AC=([^;]+)", new, flags=re.I)
    if old_accessions and new_accessions:
        return {x.strip().upper() for x in old_accessions} == {x.strip().upper() for x in new_accessions}
    return old.strip().casefold() == new.strip().casefold()


@dataclass(frozen=True, slots=True)
class SDRFChange:
    row: int
    data_file: str
    column: str
    previous: str
    proposed: str
    status: str
    evidence: str
    annotation_field: str = ""
    experiment_group: str = ""
    method: str = ""
    detail: str = ""
    support: int | None = None
    total: int | None = None


class SDRFDocument:
    """Preserve layout and apply evidence; sdrf-pipelines owns validation rules."""

    def __init__(self, columns: list[str], rows: list[list[str]], preamble: list[str] | None = None) -> None:
        self.columns = columns
        self.rows = rows
        self.preamble = preamble or []
        indices = self.indices("comment[data file]")
        if len(indices) != 1:
            raise ValueError("SDRF must have exactly one comment[data file] column.")
        self.file_column = indices[0]
        if any(len(row) != len(columns) for row in rows):
            raise ValueError(
                "SDRF contains a row with a different number of fields than its header.",
            )

    @classmethod
    def read(cls, path: Path) -> SDRFDocument:
        with path.open(encoding="utf-8-sig", newline="") as handle:
            preamble = []
            first = handle.readline()
            while first.startswith("#") or first in ("\n", "\r\n"):
                preamble.append(first)
                first = handle.readline()
            if not first:
                raise ValueError("SDRF is empty.")
            columns = next(csv.reader([first], delimiter="\t"))
            rows = list(csv.reader(handle, delimiter="\t"))
        return cls(columns, rows, preamble)

    def indices(self, column: str) -> list[int]:
        return [i for i, name in enumerate(
            self.columns,
        ) if name.strip().casefold() == column.casefold()]

    def data_files(self) -> list[str]:
        """Return file cells in row order, preserving case and repeated samples."""
        return [row[self.file_column] for row in self.rows]

    @staticmethod
    def _result_lookup(
        results: Iterable[AnalysisResult], aliases: dict[str, str] | None = None
    ) -> dict[str, AnalysisResult]:
        lookup: dict[str, AnalysisResult] = {}
        for result in results:
            names = [result.input_path.name, *result.metadata.source_files]
            if result.source_path:
                names.append(result.source_path.name)
            for name in names:
                key = file_name(name)
                if key in lookup and lookup[key] is not result:
                    raise ValueError(
                        f"Ambiguous data-file mapping: {key!r}; use unique input names.",
                    )
                lookup[key] = result
        for alias, target in (aliases or {}).items():
            if file_name(target) not in lookup:
                raise ValueError(f"File-map target was not analyzed: {target!r}")
            result = lookup[file_name(target)]
            key = file_name(alias)
            if key in lookup and lookup[key] is not result:
                raise ValueError(f"File-map alias conflicts with observed source: {key!r}")
            lookup[key] = result
        return lookup

    def resolve_results(
        self, results: Iterable[AnalysisResult], aliases: dict[str, str] | None = None
    ) -> tuple[list[AnalysisResult], list[str]]:
        """Return unique analyzed results referenced by the SDRF and missing file cells."""
        lookup = self._result_lookup(results, aliases)
        matched: list[AnalysisResult] = []
        missing: list[str] = []
        seen_results: set[int] = set()
        seen_missing: set[str] = set()
        for value in self.data_files():
            key = file_name(value)
            result = lookup.get(key)
            if result is None:
                if key not in seen_missing:
                    missing.append(value)
                    seen_missing.add(key)
                continue
            marker = id(result)
            if marker not in seen_results:
                matched.append(result)
                seen_results.add(marker)
        return matched, missing

    def _new_column_index(self) -> int:
        """Insert non-factor columns before the first factor-value column."""
        for index, name in enumerate(self.columns):
            if name.strip().casefold().startswith("factor value["):
                return index
        return len(self.columns)

    def _insert_column(self, column: str, default: str = "not available") -> int:
        index = self._new_column_index()
        self.columns.insert(index, column)
        for existing in self.rows:
            existing.insert(index, default)
        if index <= self.file_column:
            self.file_column += 1
        return index

    def set_annotation_tool(self, tool: str, version: str) -> str:
        """Set the single file-level SDRF annotation-tool provenance value."""
        tool_name = tool.strip()
        version_value = version.strip()
        if not tool_name or not version_value:
            raise ValueError("SDRF annotation tool name and version are required.")
        if not version_value.startswith("v"):
            version_value = f"v{version_value}"
        value = f"{tool_name} {version_value}"
        indices = self.indices("comment[sdrf annotation tool]")
        if len(indices) > 1:
            raise ValueError(
                "SDRF must not contain repeated comment[sdrf annotation tool] columns."
            )
        index = indices[0] if indices else self._insert_column("comment[sdrf annotation tool]")
        for row in self.rows:
            row[index] = value
        return value

    def annotate(
        self,
        results: Iterable[AnalysisResult],
        aliases: dict[str, str] | None = None,
        include_inferred: bool = False,
        overwrite: bool = False,
        *,
        include_inferred_fields: frozenset[str] = frozenset(),
        overwrite_fields: frozenset[str] = frozenset(),
        append_columns: frozenset[str] = frozenset(),
    ) -> list[SDRFChange]:
        """Apply evidence while preserving submitted values and repeated columns.

        ``include_inferred_fields`` and ``overwrite_fields`` provide a narrow
        opt-in path for confidence-gated cohort refinement without enabling or
        overwriting unrelated inferred annotations.  Columns listed in
        ``append_columns`` preserve existing assertions and use repeated SDRF
        columns for additional values.
        """
        results = list(results)
        lookup = self._result_lookup(results, aliases)
        append_column_names = {value.casefold() for value in append_columns}

        changes: list[SDRFChange] = []
        for row_number, row in enumerate(self.rows, start=2):
            name = row[self.file_column]
            matched_result = lookup.get(file_name(name))
            if matched_result is None:
                changes.append(
                    SDRFChange(
                        row_number, name, "", "", "", "unmatched",
                        "No analyzed file matched",
                    ),
                )
                continue
            for annotation in matched_result.annotations:
                column, proposed = annotation.sdrf_column, annotation.sdrf_value
                if not column or not proposed:
                    continue
                allowed_inferred = (
                    include_inferred or annotation.field in include_inferred_fields
                )
                if annotation.kind == EvidenceKind.INFERRED and not allowed_inferred:
                    existing_indices = self.indices(column)
                    matched_index = next(
                        (index for index in existing_indices if _equal(row[index], proposed)),
                        None,
                    )
                    if matched_index is not None:
                        changes.append(
                            self._change(
                                row_number,
                                name,
                                column,
                                row[matched_index],
                                proposed,
                                "match",
                                annotation,
                            )
                        )
                    else:
                        changes.append(
                            self._change(
                                row_number, name, column, "", proposed, "suggestion", annotation
                            )
                        )
                    continue
                if column.casefold() in append_column_names:
                    changes.append(
                        self._append_repeated_value(
                            row_number, row, name, column, proposed, annotation
                        )
                    )
                    continue
                indices = self.indices(column)
                if len(indices) > 1:
                    changes.append(self._change(
                        row_number, name, column, "", proposed, "ambiguous_column", annotation
                    ))
                    continue
                if not indices:
                    indices = [self._insert_column(column)]
                index = indices[0]
                previous = row[index]
                effective_overwrite = overwrite or annotation.field in overwrite_fields
                if _equal(previous, proposed):
                    status = "match"
                elif previous.strip().casefold() in MISSING:
                    row[index], status = proposed, "filled"
                elif effective_overwrite:
                    row[index], status = proposed, "replaced"
                else:
                    status = "conflict"
                changes.append(self._change(
                    row_number, name, column, previous, proposed, status, annotation
                ))
        return changes

    @staticmethod
    def _change(
        row_number: int,
        data_file: str,
        column: str,
        previous: str,
        proposed: str,
        status: str,
        annotation: Annotation,
    ) -> SDRFChange:
        field = str(getattr(annotation, "field", ""))
        value = getattr(annotation, "value", None)
        experiment_group = (
            str(value.get("experiment_group") or "") if isinstance(value, dict) else ""
        )
        kind = getattr(annotation, "kind", "")
        evidence = str(getattr(kind, "value", kind))
        return SDRFChange(
            row_number,
            data_file,
            column,
            previous,
            proposed,
            status,
            evidence,
            field,
            experiment_group,
            str(getattr(annotation, "method", "")),
            str(getattr(annotation, "detail", "")),
            getattr(annotation, "support", None),
            getattr(annotation, "total", None),
        )

    def _append_repeated_value(
        self,
        row_number: int,
        row: list[str],
        data_file: str,
        column: str,
        proposed: str,
        annotation: Annotation,
    ) -> SDRFChange:
        indices = self.indices(column)
        for index in indices:
            if _equal(row[index], proposed):
                return self._change(
                    row_number, data_file, column, row[index], proposed, "match", annotation
                )
        for index in indices:
            previous = row[index]
            if previous.strip().casefold() in MISSING:
                row[index] = proposed
                return self._change(
                    row_number, data_file, column, previous, proposed, "filled", annotation
                )
        index = self._insert_column(column)
        previous = row[index]
        row[index] = proposed
        return self._change(
            row_number, data_file, column, previous, proposed, "appended", annotation
        )

    def write(self, path: Path) -> None:
        with atomic_text(path) as handle:
            handle.writelines(self.preamble)
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(self.columns)
            writer.writerows(self.rows)
