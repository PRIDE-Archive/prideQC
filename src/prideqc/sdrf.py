"""Lossless-column SDRF editing using the standard library's TSV support."""

from __future__ import annotations

import csv
import re
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path
from urllib.parse import unquote, urlsplit

from prideqc.io import atomic_text
from prideqc.models import AnalysisResult, EvidenceKind

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

    def annotate(self, results: Iterable[AnalysisResult], aliases: dict[str, str] | None = None,
                 include_inferred: bool = False, overwrite: bool = False) -> list[SDRFChange]:
        """Fill missing cells by default; report disagreements without replacing.

        aliases maps an exact SDRF basename to the analyzed input's basename.
        Source-file names recorded by mzML are also legitimate exact aliases.
        Ambiguous mappings raise before any cell changes.
        """
        results = list(results)
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
        changes = []
        for row_number, row in enumerate(self.rows, start=2):
            name = row[self.file_column]
            matched_result = lookup.get(file_name(name))
            if matched_result is None:
                changes.append(
                    SDRFChange(
                        row_number,
                        name,
                        "",
                        "",
                        "",
                        "unmatched",
                        "No analyzed file matched",
                    ),
                )
                continue
            for annotation in matched_result.annotations:
                column, proposed = annotation.sdrf_column, annotation.sdrf_value
                if not column or not proposed:
                    continue
                if annotation.kind == EvidenceKind.INFERRED and not include_inferred:
                    changes.append(
                        SDRFChange(
                            row_number,
                            name,
                            column,
                            "",
                            proposed,
                            "suggestion",
                            annotation.kind.value,
                        ),
                    )
                    continue
                indices = self.indices(column)
                if len(indices) > 1:
                    changes.append(
                        SDRFChange(
                            row_number,
                            name,
                            column,
                            "",
                            proposed,
                            "ambiguous_column",
                            annotation.kind.value,
                        ),
                    )
                    continue
                if not indices:
                    self.columns.append(column)
                    for existing in self.rows:
                        existing.append("not available")
                    indices = [len(self.columns) - 1]
                index = indices[0]
                previous = row[index]
                if _equal(previous, proposed):
                    status = "match"
                elif previous.strip().casefold() in MISSING:
                    row[index], status = proposed, "filled"
                elif overwrite:
                    row[index], status = proposed, "replaced"
                else:
                    status = "conflict"
                changes.append(
                    SDRFChange(
                        row_number,
                        name,
                        column,
                        previous,
                        proposed,
                        status,
                        annotation.kind.value,
                    ),
                )
        return changes

    def write(self, path: Path) -> None:
        with atomic_text(path) as handle:
            handle.writelines(self.preamble)
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(self.columns)
            writer.writerows(self.rows)
