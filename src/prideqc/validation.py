"""Delegate SDRF rules and ontology checks to sdrf-pipelines."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any, Protocol


def dependency_version(distribution: str) -> str:
    """Include installed dependency versions in provenance, including editable installs."""
    try:
        return version(distribution)
    except PackageNotFoundError:
        return "unknown"


@dataclass(frozen=True, slots=True)
class ValidationReport:
    path: str
    template: str
    ontology: bool
    issues: tuple[str, ...]
    engine_version: str

    @property
    def valid(self) -> bool:
        """Conservatively require no findings from the upstream validator."""
        return not self.issues

    def to_dict(self) -> dict[str, Any]:
        return {"engine": "sdrf-pipelines", **asdict(self), "valid": self.valid}


class SDRFValidator(Protocol):
    """Validation boundary shared by the CLI and workflow."""

    def validate(self, path: Path) -> ValidationReport: ...


class SDRFPipelinesValidator:
    """Use upstream templates without copying SDRF or controlled-vocabulary rules.

    Ontology validation is opt-in because it requires the ontology extra and may
    access external services. Import and validation failures are never success.
    """

    def __init__(self, template: str = "ms-proteomics", *, ontology: bool = False) -> None:
        if not template.strip():
            raise ValueError("An SDRF template is required.")
        self.template = template
        self.ontology = ontology

    def validate(self, path: Path) -> ValidationReport:
        try:
            from sdrf_pipelines.sdrf.sdrf import read_sdrf
        except ImportError as exc:
            raise RuntimeError("SDRF validation requires sdrf-pipelines; run `uv sync`.") from exc
        try:
            document = read_sdrf(str(path))
            findings = document.validate_sdrf(
                template=self.template, skip_ontology=not self.ontology
            )
            issues = tuple(str(getattr(item, "message", item)) for item in findings)
        except Exception as exc:
            hint = " Install ontology support with `uv sync --extra ontology`." if self.ontology else ""
            raise RuntimeError(f"sdrf-pipelines could not validate {path}: {exc}.{hint}") from exc
        return ValidationReport(
            str(path.resolve()), self.template, self.ontology, issues,
            dependency_version("sdrf-pipelines"),
        )
