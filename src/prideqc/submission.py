"""Deterministic application of frozen SDRF adjudication decisions to a full SDRF.

This is intentionally separate from candidate generation and LLM adjudication.  It takes
one already-validated adjudication request/response pair, rebinds the request's target
runs to the canonical full SDRF by exact file semantics, verifies that the referenced
original cells have not drifted, applies only accepted supplied candidates, validates the
full output, and only then stamps prideQC annotation-tool provenance.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import shutil
import sys
from collections import Counter, defaultdict
from collections.abc import Mapping, Sequence
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

from prideqc import __version__
from prideqc.io import atomic_text, write_json
from prideqc.refinement_adjudication import (
    validate_llm_adjudication_request,
    validate_llm_refinement_decisions,
)
from prideqc.sdrf import MISSING, SDRFDocument, _equal, file_name
from prideqc.validation import SDRFPipelinesValidator, SDRFValidator

APPLICATION_SCHEMA_VERSION = "prideqc-sdrf-adjudication-application-v1"
ANNOTATION_TOOL_COLUMN = "comment[sdrf annotation tool]"
TECHNOLOGY_TYPE_COLUMN = "technology type"
ASSAY_NAME_COLUMN = "assay name"


@dataclass(frozen=True, slots=True)
class DecisionApplication:
    row: int
    data_file: str
    target_run: str
    decision_id: str
    decision: str
    target_field: str
    previous: str
    selected_value: str
    status: str
    experiment_group: str
    reason: str


@dataclass(frozen=True, slots=True)
class StructuralChange:
    kind: str
    column: str
    previous_position: int
    new_position: int
    detail: str


@dataclass(frozen=True, slots=True)
class ProvenanceChange:
    row: int
    data_file: str
    column: str
    previous: str
    current: str


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _load_object(path: Path, *, label: str) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"{label} must contain one JSON object")
    return value


def _load_file_map(path: Path | None) -> dict[str, str]:
    if path is None:
        return {}
    value = _load_object(path, label="--file-map")
    mapping: dict[str, str] = {}
    for source, target in value.items():
        if not isinstance(source, str) or not isinstance(target, str):
            raise ValueError("--file-map keys and values must be strings")
        source_name = file_name(source)
        target_name = file_name(target)
        if not source_name or not target_name:
            raise ValueError("--file-map keys and values must have non-empty basenames")
        if source_name in mapping and mapping[source_name] != target_name:
            raise ValueError(f"Conflicting --file-map source: {source_name!r}")
        mapping[source_name] = target_name
    return mapping


def _mapped_run(value: str, aliases: Mapping[str, str]) -> str:
    name = file_name(value)
    return aliases.get(name, name)


def _request_decisions(request: Mapping[str, Any]) -> dict[str, Mapping[str, Any]]:
    raw = request.get("decisions")
    if not isinstance(raw, list):
        return {}
    return {
        str(item["decision_id"]): item
        for item in raw
        if isinstance(item, Mapping) and item.get("decision_id")
    }


def _response_decisions(response: Mapping[str, Any]) -> dict[str, Mapping[str, Any]]:
    raw = response.get("decisions")
    if not isinstance(raw, list):
        return {}
    return {
        str(item["decision_id"]): item
        for item in raw
        if isinstance(item, Mapping) and item.get("decision_id")
    }


def _target_rows_by_run(
    document: SDRFDocument,
    target_runs: Sequence[str],
    aliases: Mapping[str, str],
) -> dict[str, list[int]]:
    targets = [file_name(str(value)) for value in target_runs]
    if not targets or any(not value for value in targets):
        raise ValueError("Adjudication decision must contain non-empty target runs")
    if len(set(targets)) != len(targets):
        raise ValueError("Adjudication decision target runs must be unique")

    rows: dict[str, list[int]] = {target: [] for target in targets}
    target_set = set(targets)
    for row_index, row in enumerate(document.rows):
        run = _mapped_run(row[document.file_column], aliases)
        if run in target_set:
            rows[run].append(row_index)

    missing = [target for target in targets if not rows[target]]
    if missing:
        raise ValueError(
            "Full SDRF does not contain adjudication target run(s): " + ", ".join(missing)
        )
    return rows


def _row_values(document: SDRFDocument, row_index: int, column: str) -> tuple[str, ...]:
    return tuple(document.rows[row_index][index] for index in document.indices(column))


def _verify_original_context(
    document: SDRFDocument,
    decision: Mapping[str, Any],
    rows_by_run: Mapping[str, Sequence[int]],
    aliases: Mapping[str, str],
) -> None:
    decision_id = str(decision.get("decision_id") or "")
    target_field = str(decision.get("target_field") or "")
    original = decision.get("original")
    if not isinstance(original, Mapping):
        raise ValueError(f"Decision {decision_id} is missing original SDRF context")

    expected_column_count = original.get("column_count")
    if not isinstance(expected_column_count, int):
        raise ValueError(f"Decision {decision_id} has invalid original column count")
    observed_indices = document.indices(target_field)
    if len(observed_indices) != expected_column_count:
        raise ValueError(
            f"Decision {decision_id} full-SDRF column drift for {target_field!r}: "
            f"expected {expected_column_count}, observed {len(observed_indices)}"
        )

    raw_rows = original.get("rows")
    if not isinstance(raw_rows, list):
        raise ValueError(f"Decision {decision_id} has invalid original row context")

    expected: defaultdict[str, list[tuple[str, ...]]] = defaultdict(list)
    for item in raw_rows:
        if not isinstance(item, Mapping):
            raise ValueError(f"Decision {decision_id} original row context must be objects")
        data_file = item.get("data_file")
        values = item.get("values")
        if not isinstance(data_file, str) or not isinstance(values, list):
            raise ValueError(f"Decision {decision_id} has malformed original row context")
        if any(not isinstance(value, str) for value in values):
            raise ValueError(f"Decision {decision_id} original values must be strings")
        expected[_mapped_run(data_file, aliases)].append(tuple(values))

    for run, row_indices in rows_by_run.items():
        expected_values = Counter(expected.get(run, []))
        observed_values = Counter(
            _row_values(document, row_index, target_field) for row_index in row_indices
        )
        if expected_values != observed_values:
            raise ValueError(
                f"Decision {decision_id} original SDRF values drifted for run {run!r}: "
                f"expected {sorted(expected_values.elements())!r}, "
                f"observed {sorted(observed_values.elements())!r}"
            )


def _insert_column(document: SDRFDocument, column: str, default: str = "not available") -> int:
    # Keep the same placement semantics as SDRFDocument annotation: non-factor metadata
    # goes immediately before the first factor-value column.
    insert_at = len(document.columns)
    for candidate_index, name in enumerate(document.columns):
        if name.strip().casefold().startswith("factor value["):
            insert_at = candidate_index
            break
    document.columns.insert(insert_at, column)
    for row in document.rows:
        row.insert(insert_at, default)
    if insert_at <= document.file_column:
        document.file_column += 1
    return insert_at


def _display_values(values: Sequence[str]) -> str:
    return " | ".join(values)


def _apply_single_value(
    document: SDRFDocument,
    *,
    row_index: int,
    target_field: str,
    selected: str,
) -> tuple[str, str]:
    indices = document.indices(target_field)
    if len(indices) > 1:
        raise ValueError(f"Target field {target_field!r} is ambiguous in the full SDRF")
    if not indices:
        indices = [_insert_column(document, target_field)]
    index = indices[0]
    previous = document.rows[row_index][index]
    if _equal(previous, selected):
        return previous, "accepted_match"
    document.rows[row_index][index] = selected
    status = "accepted_filled" if previous.strip().casefold() in MISSING else "accepted_replaced"
    return previous, status


def _apply_repeated_value(
    document: SDRFDocument,
    *,
    row_index: int,
    target_field: str,
    selected: str,
) -> tuple[str, str]:
    row = document.rows[row_index]
    indices = document.indices(target_field)
    for index in indices:
        if _equal(row[index], selected):
            return row[index], "accepted_match"
    for index in indices:
        previous = row[index]
        if previous.strip().casefold() in MISSING:
            row[index] = selected
            return previous, "accepted_filled"
    index = _insert_column(document, target_field)
    previous = document.rows[row_index][index]
    document.rows[row_index][index] = selected
    return previous, "accepted_appended"


def _move_unique_column_after(
    document: SDRFDocument,
    *,
    column: str,
    predecessor: str,
) -> StructuralChange | None:
    indices = document.indices(column)
    predecessor_indices = document.indices(predecessor)
    if len(indices) != 1 or len(predecessor_indices) != 1:
        raise ValueError(
            f"Column-order normalization requires exactly one {column!r} and one "
            f"{predecessor!r} column"
        )
    source = indices[0]
    predecessor_index = predecessor_indices[0]
    if source == predecessor_index + 1:
        return None

    column_name = document.columns.pop(source)
    values = [row.pop(source) for row in document.rows]
    if source < predecessor_index:
        predecessor_index -= 1
    destination = predecessor_index + 1
    document.columns.insert(destination, column_name)
    for row, value in zip(document.rows, values, strict=True):
        row.insert(destination, value)

    file_indices = document.indices("comment[data file]")
    if len(file_indices) != 1:
        raise ValueError("SDRF must retain exactly one comment[data file] column")
    document.file_column = file_indices[0]
    return StructuralChange(
        "move_column",
        column_name,
        source + 1,
        destination + 1,
        f"Moved {column_name!r} immediately after {document.columns[predecessor_index]!r}",
    )


def _write_application_table(path: Path, records: Sequence[DecisionApplication]) -> None:
    fields = [field.name for field in DecisionApplication.__dataclass_fields__.values()]
    with atomic_text(path) as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(asdict(record) for record in records)


def _write_structural_table(path: Path, records: Sequence[StructuralChange]) -> None:
    fields = [field.name for field in StructuralChange.__dataclass_fields__.values()]
    with atomic_text(path) as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(asdict(record) for record in records)


def _write_provenance_table(path: Path, records: Sequence[ProvenanceChange]) -> None:
    fields = [field.name for field in ProvenanceChange.__dataclass_fields__.values()]
    with atomic_text(path) as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(asdict(record) for record in records)


def _project_accession_guard(path: Path, request: Mapping[str, Any]) -> str:
    accession = str(request.get("project_accession") or "").strip().upper()
    if not re.fullmatch(r"PXD\d+", accession):
        raise ValueError("Adjudication request must contain one PXD project accession")
    observed = {match.upper() for match in re.findall(r"PXD\d+", str(path), re.I)}
    if observed and observed != {accession}:
        raise ValueError(
            f"Full SDRF path accession does not match adjudication request {accession}: "
            + ", ".join(sorted(observed))
        )
    return accession


def apply_adjudicated_sdrf(
    *,
    sdrf: Path | str,
    request_path: Path | str,
    decisions_path: Path | str,
    output_directory: Path | str,
    file_map: Path | str | None = None,
    normalize_technology_type_order: bool = False,
    overwrite: bool = False,
    validator: SDRFValidator | None = None,
) -> dict[str, Any]:
    """Apply frozen adjudication decisions to a canonical full SDRF.

    Benchmark row numbers are never used.  Each request decision is rebound to the full
    SDRF through ``target_runs`` and optional exact file-map aliases.  Before mutation,
    the target field values are compared against the request's persisted original SDRF
    context; any drift aborts the application.
    """
    source = Path(sdrf).resolve(strict=True)
    request_file = Path(request_path).resolve(strict=True)
    decisions_file = Path(decisions_path).resolve(strict=True)
    output = Path(output_directory).resolve()
    file_map_path = Path(file_map).resolve(strict=True) if file_map is not None else None

    if output.exists() and not output.is_dir():
        raise FileExistsError(f"Output path is not a directory: {output}")
    if output.exists() and any(output.iterdir()):
        if not overwrite:
            raise FileExistsError(
                f"Output directory is not empty: {output}. "
                "Choose a new directory or pass --overwrite."
            )
        for child in output.iterdir():
            if child.is_dir() and not child.is_symlink():
                shutil.rmtree(child)
            else:
                child.unlink()
    output.mkdir(parents=True, exist_ok=True)

    request = _load_object(request_file, label="--request")
    decisions = _load_object(decisions_file, label="--decisions")
    validate_llm_adjudication_request(request)
    validate_llm_refinement_decisions(decisions, request)
    if request.get("input_mode") != "sdrf-backed":
        raise ValueError("Full-SDRF application requires an sdrf-backed adjudication request")
    accession = _project_accession_guard(source, request)
    aliases = _load_file_map(file_map_path)

    validator_impl = validator or SDRFPipelinesValidator()
    input_validation = validator_impl.validate(source)
    shutil.copyfile(source, output / "original.sdrf.tsv")
    document = SDRFDocument.read(source)

    request_by_id = _request_decisions(request)
    response_by_id = _response_decisions(decisions)

    # Verify the entire frozen application plan before mutating any cell.
    row_plan: dict[str, dict[str, list[int]]] = {}
    for decision_id in sorted(request_by_id, key=str.casefold):
        item = request_by_id[decision_id]
        raw_targets = item.get("target_runs")
        if not isinstance(raw_targets, list) or any(not isinstance(v, str) for v in raw_targets):
            raise ValueError(f"Decision {decision_id} has invalid target runs")
        rows_by_run = _target_rows_by_run(document, raw_targets, aliases)
        _verify_original_context(document, item, rows_by_run, aliases)
        row_plan[decision_id] = rows_by_run

    structural_changes: list[StructuralChange] = []
    if normalize_technology_type_order:
        change = _move_unique_column_after(
            document,
            column=TECHNOLOGY_TYPE_COLUMN,
            predecessor=ASSAY_NAME_COLUMN,
        )
        if change is not None:
            structural_changes.append(change)

    applications: list[DecisionApplication] = []
    changed_cell_count = 0
    decision_counts = Counter(
        str(item.get("decision") or "") for item in response_by_id.values()
    )

    for decision_id in sorted(request_by_id, key=str.casefold):
        request_decision = request_by_id[decision_id]
        response_decision = response_by_id[decision_id]
        final_decision = str(response_decision.get("decision") or "")
        selected = response_decision.get("selected_value")
        reason = str(response_decision.get("reason") or "")
        target_field = str(request_decision.get("target_field") or "")
        experiment_group = str(request_decision.get("experiment_group") or "")
        write_semantics = str(request_decision.get("write_semantics") or "")
        rows_by_run = row_plan[decision_id]

        for target_run in sorted(rows_by_run, key=str.casefold):
            for row_index in rows_by_run[target_run]:
                row_number = row_index + 2
                data_file = document.rows[row_index][document.file_column]
                before_values = _row_values(document, row_index, target_field)

                if final_decision == "accept":
                    assert isinstance(selected, str)
                    if write_semantics == "fill_or_replace_canonical_value":
                        previous, status = _apply_single_value(
                            document,
                            row_index=row_index,
                            target_field=target_field,
                            selected=selected,
                        )
                    elif (
                        write_semantics
                        == "append_one_allowed_canonical_value_if_accepted_and_missing"
                    ):
                        previous, status = _apply_repeated_value(
                            document,
                            row_index=row_index,
                            target_field=target_field,
                            selected=selected,
                        )
                    else:
                        raise ValueError(
                            f"Decision {decision_id} has unsupported write semantics: "
                            f"{write_semantics!r}"
                        )
                    if status in {"accepted_filled", "accepted_replaced", "accepted_appended"}:
                        changed_cell_count += 1
                    selected_text = selected
                else:
                    previous = _display_values(before_values)
                    status = f"{final_decision}_preserved"
                    selected_text = ""

                applications.append(
                    DecisionApplication(
                        row_number,
                        data_file,
                        target_run,
                        decision_id,
                        final_decision,
                        target_field,
                        previous,
                        selected_text,
                        status,
                        experiment_group,
                        reason,
                    )
                )

    _write_application_table(output / "decision-application.tsv", applications)
    _write_structural_table(output / "structural-changes.tsv", structural_changes)

    changed_before_provenance = changed_cell_count > 0 or bool(structural_changes)
    candidate_path = output / "candidate.sdrf.tsv"
    document.write(candidate_path)
    candidate_validation = validator_impl.validate(candidate_path)

    annotation_tool: str | None = None
    provenance_changes: list[ProvenanceChange] = []
    final_validation = None
    submission_path: Path | None = None
    if candidate_validation.valid:
        if changed_before_provenance:
            existing_tool_indices = document.indices(ANNOTATION_TOOL_COLUMN)
            previous_tool_values = [
                (
                    row_index + 2,
                    row[document.file_column],
                    row[existing_tool_indices[0]] if existing_tool_indices else "",
                )
                for row_index, row in enumerate(document.rows)
            ]
            annotation_tool = document.set_annotation_tool("prideQC", __version__)
            for row_number, data_file, previous in previous_tool_values:
                if not _equal(previous, annotation_tool):
                    provenance_changes.append(
                        ProvenanceChange(
                            row_number,
                            data_file,
                            ANNOTATION_TOOL_COLUMN,
                            previous,
                            annotation_tool,
                        )
                    )
        final_path = output / f"{accession}.sdrf.tsv"
        document.write(final_path)
        final_validation = validator_impl.validate(final_path)
        if final_validation.valid:
            submission_path = final_path

    _write_provenance_table(output / "provenance-changes.tsv", provenance_changes)

    submission_ready = bool(
        changed_before_provenance
        and candidate_validation.valid
        and final_validation is not None
        and final_validation.valid
    )
    manifest = {
        "schema_version": APPLICATION_SCHEMA_VERSION,
        "prideqc_version": __version__,
        "project_accession": accession,
        "source_sdrf": {
            "path": str(source),
            "sha256": _sha256(source),
            "validation": input_validation.to_dict(),
        },
        "frozen_adjudication": {
            "request": str(request_file),
            "request_sha256": _sha256(request_file),
            "request_id": request.get("request_id"),
            "source_packet": request.get("source_packet"),
            "benchmark_sdrf": request.get("sdrf"),
            "decisions": str(decisions_file),
            "decisions_sha256": _sha256(decisions_file),
            "decision_counts": {
                "accept": decision_counts.get("accept", 0),
                "reject": decision_counts.get("reject", 0),
                "abstain": decision_counts.get("abstain", 0),
            },
        },
        "file_map": str(file_map_path) if file_map_path is not None else None,
        "normalize_technology_type_order": normalize_technology_type_order,
        "decision_application": "decision-application.tsv",
        "structural_changes": [asdict(item) for item in structural_changes],
        "provenance_changes": "provenance-changes.tsv",
        "provenance_change_count": len(provenance_changes),
        "changed_cell_count": changed_cell_count,
        "changed_before_provenance": changed_before_provenance,
        "candidate_sdrf": candidate_path.name,
        "candidate_validation": candidate_validation.to_dict(),
        "annotation_tool": annotation_tool,
        "submission_sdrf": submission_path.name if submission_path is not None else None,
        "submission_sha256": _sha256(submission_path) if submission_path is not None else None,
        "final_validation": final_validation.to_dict() if final_validation is not None else None,
        "submission_ready": submission_ready,
        "success": bool(candidate_validation.valid and final_validation and final_validation.valid),
    }
    write_json(output / "manifest.json", manifest)
    return manifest


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(
        prog="python -m prideqc.submission",
        description="Apply frozen prideQC adjudication decisions to a canonical full SDRF",
    )
    result.add_argument("--sdrf", required=True, type=Path, help="Canonical full original SDRF")
    result.add_argument(
        "--request",
        required=True,
        type=Path,
        help="Frozen adjudication request JSON",
    )
    result.add_argument("--decisions", required=True, type=Path, help="Frozen final decision JSON")
    result.add_argument("-o", "--output-dir", required=True, type=Path)
    result.add_argument(
        "--file-map",
        type=Path,
        help="Optional JSON object mapping exact SDRF filenames to adjudication target-run names",
    )
    result.add_argument(
        "--normalize-technology-type-order",
        action="store_true",
        help="Explicitly move the unique technology type column immediately after assay name",
    )
    result.add_argument("--sdrf-template", default="ms-proteomics")
    result.add_argument("--validate-ontology", action="store_true")
    result.add_argument("--overwrite", action="store_true")
    return result


def main(argv: Sequence[str] | None = None) -> int:
    arguments = parser().parse_args(argv)
    try:
        manifest = apply_adjudicated_sdrf(
            sdrf=arguments.sdrf,
            request_path=arguments.request,
            decisions_path=arguments.decisions,
            output_directory=arguments.output_dir,
            file_map=arguments.file_map,
            normalize_technology_type_order=arguments.normalize_technology_type_order,
            overwrite=arguments.overwrite,
            validator=SDRFPipelinesValidator(
                arguments.sdrf_template,
                ontology=arguments.validate_ontology,
            ),
        )
    except (OSError, ValueError, RuntimeError) as exc:
        print(f"prideqc submission: {exc}", file=sys.stderr)
        return 2

    print(
        f"Applied frozen adjudication for {manifest['project_accession']}: "
        f"changed_cells={manifest['changed_cell_count']} "
        f"submission_ready={str(manifest['submission_ready']).lower()}"
    )
    print(f"Results: {arguments.output_dir}")
    return 0 if manifest["success"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
