"""Small domain objects; readers and analyzers do not depend on each other."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from enum import StrEnum
from pathlib import Path
from typing import Any, Protocol

import numpy as np
from numpy.typing import NDArray

FloatArray = NDArray[np.float64]


@dataclass(frozen=True, slots=True)
class CVTerm:
    accession: str
    name: str

    def sdrf_value(self) -> str:
        return f"NT={self.name};AC={self.accession}"


class EvidenceKind(StrEnum):
    OBSERVED = "observed"
    INFERRED = "inferred"
    UNAVAILABLE = "unavailable"


@dataclass(frozen=True, slots=True)
class Annotation:
    """A measurement or suggestion, with its method and support kept explicit."""

    field: str
    value: Any
    kind: EvidenceKind
    method: str
    detail: str = ""
    support: int | None = None
    total: int | None = None
    sdrf_column: str | None = None
    sdrf_value: str | None = None


@dataclass(frozen=True, slots=True)
class Precursor:
    mz: float
    charge: int = 0
    intensity: float = 0.0
    isolation_width: float | None = None
    activation: tuple[CVTerm, ...] = ()
    collision_energy: tuple[float, str] | None = None


@dataclass(frozen=True, slots=True)
class Spectrum:
    """Ephemeral arrays: consumers must not retain the spectrum or its peaks."""

    ms_level: int
    rt: float
    mz: FloatArray
    intensity: FloatArray
    precursors: tuple[Precursor, ...] = ()
    polarity: str = "unknown"
    representation: str = "unknown"
    scan_windows: tuple[tuple[float, float], ...] = ()
    faims_cv: float | None = None
    native_id: str = ""
    estimated_representation: str | None = None


@dataclass(slots=True)
class RunMetadata:
    instruments: list[CVTerm] = field(default_factory=list)
    analyzers: list[CVTerm] = field(default_factory=list)
    ionization: list[CVTerm] = field(default_factory=list)
    serial_numbers: list[str] = field(default_factory=list)
    source_files: list[str] = field(default_factory=list)
    started_at: str | None = None
    # Native vendor readers expose some instrument fields as plain strings
    # rather than PSI-MS CV terms. Preserve those values without fabricating
    # accessions from model-name heuristics.
    instrument_details: dict[str, str] = field(default_factory=dict)


@dataclass(frozen=True, slots=True)
class Metric:
    key: str
    value: Any


@dataclass(slots=True)
class AnalysisResult:
    input_path: Path
    metadata: RunMetadata
    metrics: list[Metric]
    annotations: list[Annotation]
    warnings: list[str]
    engine_version: str
    elapsed_seconds: float
    source_path: Path | None = None

    def to_dict(self) -> dict[str, Any]:
        result = asdict(self)
        result["input_path"] = str(self.input_path)
        result["source_path"] = str(self.source_path) if self.source_path else None
        return result


class SpectrumSink(Protocol):
    def consume_spectrum(self, spectrum: Spectrum) -> None: ...

    def consume_chromatogram(self, rt: FloatArray, kind: int) -> None: ...


class SpectrumReader(Protocol):
    """A vendor or open-format adapter streams into the same analysis sink."""

    engine_version: str

    def read(self, path: Path, sink: SpectrumSink) -> RunMetadata: ...


class EvidenceCollector(Protocol):
    """Optional estimators consume the same peak arrays without a second read."""

    def consume_spectrum(self, spectrum: Spectrum) -> None: ...

    def annotations(self) -> list[Annotation]: ...
