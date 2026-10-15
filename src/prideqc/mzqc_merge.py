"""Lossless, fail-closed consolidation of per-run mzQC 1.0 documents.

This is packaging, not scientific aggregation: no setQualities or recalculated
quality metrics are emitted. Per-run JSON subtrees are preserved verbatim.
"""

from __future__ import annotations

import hashlib
import json
from datetime import UTC, datetime
from pathlib import Path
from typing import Any
from urllib.parse import urlparse

from prideqc.io import atomic_text

LOCAL_VOCABULARY = "prideqc local vocabulary"


def _digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _terms(text: str, source: Path) -> dict[str, str]:
    """Collect OBO term stanzas by accession, detecting definitions that disagree."""
    _require("format-version:" in text, f"Invalid companion OBO header: {source}")
    definitions: dict[str, str] = {}
    for segment in text.replace("\r\n", "\n").split("\n[Term]\n")[1:]:
        block = "[Term]\n" + segment.strip() + "\n"
        accession = next(
            (line[4:].strip() for line in block.splitlines() if line.startswith("id: ")),
            "",
        )
        _require(bool(accession), f"OBO stanza without id in {source}")
        if accession in definitions:
            _require(definitions[accession] == block, f"Conflicting OBO term {accession}: {source}")
        definitions[accession] = block
    _require(bool(definitions), f"No OBO terms in {source}")
    return definitions


def _check_run(run: Any, source: Path) -> tuple[str, list[str]]:
    _require(isinstance(run, dict), f"Invalid runQuality object: {source}")
    metadata = run.get("metadata")
    _require(isinstance(metadata, dict), f"Missing run metadata: {source}")
    label = metadata.get("label")
    _require(isinstance(label, str) and bool(label.strip()), f"Missing run label: {source}")
    files = metadata.get("inputFiles")
    _require(isinstance(files, list) and bool(files), f"Missing inputFiles: {source}")
    software = metadata.get("analysisSoftware")
    _require(isinstance(software, list) and bool(software), f"Missing analysisSoftware: {source}")
    metrics = run.get("qualityMetrics")
    _require(isinstance(metrics, list) and bool(metrics), f"Missing qualityMetrics: {source}")
    names: list[str] = []
    for f in files:
        _require(isinstance(f, dict), f"Invalid inputFile: {source}")
        name = f.get("name")
        _require(isinstance(name, str) and bool(name.strip()), f"Missing inputFile.name: {source}")
        _require(isinstance(f.get("location"), str), f"Missing inputFile.location: {source}")
        _require(isinstance(f.get("fileFormat"), dict), f"Missing inputFile.fileFormat: {source}")
        names.append(name)
    for metric in metrics:
        _require(isinstance(metric, dict), f"Invalid qualityMetric in {source}")
        _require(
            isinstance(metric.get("accession"), str)
            and bool(metric["accession"])
            and isinstance(metric.get("name"), str)
            and bool(metric["name"])
            and "value" in metric,
            f"Invalid qualityMetric properties in {source}",
        )
    return label, names


def merge_mzqc(
    input_dir: Path,
    output: Path,
    *,
    accession: str,
    expected_runs: int | None = None,
    cv_uri: str | None = None,
) -> dict[str, Any]:
    """Merge exactly one runQuality per per-run file; refuse ambiguous inputs.

    The output .obo and .manifest.json live beside the merged .mzQC. Passing
    cv_uri is useful inside containers: it should resolve to the *host-side*
    final .obo path or a URL where that sidecar will actually be published.
    """
    _require(input_dir.is_dir(), f"Missing per-run input directory: {input_dir}")
    _require(
        not output.resolve().is_relative_to(input_dir.resolve()),
        "Merged output must be outside the per-run input directory",
    )
    _require(accession.startswith("PXD") and accession[3:].isdigit(), "Invalid PXD accession")
    _require(expected_runs is None or expected_runs > 0, "expected_runs must be positive")
    obo_output = output.with_suffix(".obo")
    manifest_output = output.with_suffix(".manifest.json")
    cv_uri = cv_uri or obo_output.absolute().as_uri()
    uri = urlparse(cv_uri)
    _require(bool(uri.scheme), "The merged CV URI must be an absolute URI")
    files = sorted(input_dir.rglob("*.mzQC"), key=lambda p: p.relative_to(input_dir).as_posix())
    _require(bool(files), f"No per-run mzQC found in {input_dir}")
    _require(
        expected_runs is None or len(files) == expected_runs,
        (
            f"Incomplete accession {accession}: expected {expected_runs}, "
            f"found {len(files)} mzQC files"
        ),
    )
    runs: list[dict[str, Any]] = []
    vocabularies: dict[str, dict[str, str]] = {}
    terms: dict[str, str] = {}
    labels: set[str] = set()
    input_names: set[str] = set()
    member_manifest: list[dict[str, str]] = []

    for path in files:
        contents = path.read_bytes()
        try:
            data = json.loads(contents)
        except (UnicodeError, json.JSONDecodeError) as exc:
            raise ValueError(f"Invalid JSON: {path}: {exc}") from exc
        _require(isinstance(data, dict) and set(data) == {"mzQC"}, f"Invalid mzQC root: {path}")
        body = data["mzQC"]
        _require(
            isinstance(body, dict) and body.get("version") == "1.0.0",
            f"Unsupported mzQC version: {path}",
        )
        _require("setQualities" not in body, f"Unexpected setQualities in per-run file: {path}")
        file_runs = body.get("runQualities")
        _require(
            isinstance(file_runs, list) and len(file_runs) == 1,
            f"Expected exactly one runQuality: {path}",
        )
        date_text = body.get("creationDate")
        _require(isinstance(date_text, str), f"Missing creationDate: {path}")
        try:
            date = datetime.fromisoformat(date_text.replace("Z", "+00:00"))
        except ValueError as exc:
            raise ValueError(f"Invalid creationDate in {path}: {date_text}") from exc
        _require(date.tzinfo is not None, f"mzQC creationDate lacks timezone: {path}")
        run = file_runs[0]
        label, names = _check_run(run, path)
        _require(label not in labels, f"Duplicate mzQC run label {label!r} at {path}")
        labels.add(label)
        for name in names:
            _require(name not in input_names, f"Duplicate inputFile.name {name!r} at {path}")
            input_names.add(name)
            for f in run["metadata"]["inputFiles"]:
                if f["name"] != name:
                    continue
                for prop in f.get("fileProperties", []):
                    if isinstance(prop, dict) and prop.get("accession") == "MS:1001919":
                        _require(
                            prop.get("value") == accession,
                            f"ProteomeXchange accession mismatch in {path}: {prop.get('value')}",
                        )
        cvs = body.get("controlledVocabularies")
        _require(isinstance(cvs, list) and bool(cvs), f"Missing controlledVocabularies: {path}")
        has_local = False
        for cv in cvs:
            _require(isinstance(cv, dict), f"Invalid controlledVocabulary: {path}")
            name, source_uri = cv.get("name"), cv.get("uri")
            _require(isinstance(name, str) and isinstance(source_uri, str), f"Invalid CV in {path}")
            if name == LOCAL_VOCABULARY:
                _require(
                    cv.get("version") == "1",
                    f"Unsupported prideqc local vocabulary version in {path}",
                )
                has_local = True
                continue
            _require(
                not source_uri.startswith("file:"),
                f"Unknown file-based controlled vocabulary {name!r} in {path}",
            )
            previous_cv = vocabularies.setdefault(name, cv)
            _require(previous_cv == cv, f"Conflicting controlled vocabulary {name!r} at {path}")
        _require(has_local, f"Missing {LOCAL_VOCABULARY} declaration: {path}")
        obo_source = path.with_suffix(".obo")
        _require(obo_source.is_file(), f"Missing companion OBO: {obo_source}")
        obo_bytes = obo_source.read_bytes()
        try:
            source_terms = _terms(obo_bytes.decode("utf-8"), obo_source)
        except UnicodeError as exc:
            raise ValueError(f"Invalid OBO encoding: {obo_source}") from exc
        for key, definition in source_terms.items():
            if key in terms:
                _require(
                    terms[key] == definition,
                    f"Conflicting QCPRIDE OBO definition {key} in {obo_source}",
                )
            terms[key] = definition
        for metric in run["qualityMetrics"]:
            accession_id = metric["accession"]
            if accession_id.startswith("QCPRIDE:"):
                _require(
                    accession_id in source_terms,
                    f"Missing QCPRIDE term {accession_id} in {obo_source}",
                )
        runs.append(run)
        member_manifest.append(
            {
                "mzqc": path.relative_to(input_dir).as_posix(),
                "mzqc_sha256": _digest(contents),
                "obo": obo_source.relative_to(input_dir).as_posix(),
                "obo_sha256": _digest(obo_bytes),
            }
        )

    _require(
        "QCPRIDE:SOFTWARE" in terms and "QCPRIDE:SOURCEFILE" in terms,
        "Missing required QCPRIDE metadata definitions",
    )
    vocabularies[LOCAL_VOCABULARY] = {
        "name": LOCAL_VOCABULARY,
        "version": "1",
        "uri": cv_uri,
    }
    manifest_base = {
        "accession": accession,
        "status": "complete" if expected_runs is not None else "count_unverified",
        "expected_runs": expected_runs,
        "merged_runs": len(runs),
        "cv_uri": cv_uri,
        "members": member_manifest,
    }
    if manifest_output.is_file() and output.is_file() and obo_output.is_file():
        previous_manifest: object
        try:
            previous_manifest = json.loads(manifest_output.read_text(encoding="utf-8"))
        except (ValueError, OSError):
            previous_manifest = None
        if isinstance(previous_manifest, dict) and all(
            previous_manifest.get(key) == value for key, value in manifest_base.items()
        ):
            if previous_manifest.get("mzqc_sha256") == _digest(
                output.read_bytes()
            ) and previous_manifest.get("obo_sha256") == _digest(obo_output.read_bytes()):
                return previous_manifest

    created_utc = datetime.now(UTC).isoformat()
    merged = {
        "mzQC": {
            "version": "1.0.0",
            "creationDate": created_utc,
            "description": (
                f"Lossless prideQC study consolidation for {accession}; "
                "per-run quality metrics unchanged; no setQualities computed. "
                "See companion manifest for file checksums."
            ),
            "runQualities": runs,
            "controlledVocabularies": [vocabularies[name] for name in sorted(vocabularies)],
        }
    }
    obo_text = "format-version: 1.2\ndata-version: 1\nontology: prideqc\n" + "".join(
        "\n" + terms[key] for key in sorted(terms)
    )
    mzqc_text = json.dumps(merged, indent=2, ensure_ascii=False, allow_nan=False) + "\n"
    mzqc_sha = _digest(mzqc_text.encode("utf-8"))
    obo_sha = _digest(obo_text.encode("utf-8"))
    manifest = {
        **manifest_base,
        "merged_utc": created_utc,
        "mzqc_sha256": mzqc_sha,
        "obo_sha256": obo_sha,
    }
    # Publish manifest last: it acts as an integrity/completeness marker in the
    # event of a crash between publishing the sidecar and the merged document.
    with atomic_text(obo_output) as handle:
        handle.write(obo_text)
    with atomic_text(output) as handle:
        handle.write(mzqc_text)
    with atomic_text(manifest_output) as handle:
        json.dump(manifest, handle, indent=2, ensure_ascii=False)
        handle.write("\n")
    return manifest
