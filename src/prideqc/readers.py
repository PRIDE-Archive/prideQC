"""pyOpenMS streaming adapter and a bounded mzML header reader."""

from __future__ import annotations

import gzip
import math
import shutil
import tempfile
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any

import numpy as np

from prideqc.models import CVTerm, Precursor, RunMetadata, Spectrum, SpectrumSink

# OpenMS enum *names*, never unexplained integer constants across API releases.
ACTIVATION = {
    "CID": CVTerm("MS:1000133", "collision-induced dissociation"),
    "HCID": CVTerm("MS:1000422", "beam-type collision-induced dissociation"),
    "HCD": CVTerm("MS:1000422", "beam-type collision-induced dissociation"),
    "ETD": CVTerm("MS:1000598", "electron transfer dissociation"),
    "ECD": CVTerm("MS:1000250", "electron capture dissociation"),
    "IRMPD": CVTerm("MS:1000262", "infrared multiphoton dissociation"),
    "ETHCD": CVTerm("MS:1002631", "electron-transfer/higher-energy collision dissociation"),
    "ETCID": CVTerm("MS:1002632", "electron-transfer/collision-induced dissociation"),
}


def _tag(element: ET.Element) -> str:
    return element.tag.rsplit("}", 1)[-1]


def read_header(path: Path) -> RunMetadata:
    """Preserve instrument CV accessions instead of guessing from model substrings.

    Stop before the first spectrum/chromatogram array. This bounded header read
    does not decode or traverse the spectrum payload a second time.
    """
    metadata = RunMetadata()
    groups: dict[str, list[dict[str, str]]] = {}

    def terms(node: ET.Element) -> list[dict[str, str]]:
        result = []
        for child in node:
            if _tag(child) == "cvParam":
                result.append(child.attrib)
            elif _tag(child) == "referenceableParamGroupRef":
                result.extend(groups.get(child.get("ref", ""), []))
        return result

    with path.open("rb") as handle:
        for event, element in ET.iterparse(handle, events=("start", "end")):
            tag = _tag(element)
            if event == "start":
                if tag in {"spectrumList", "chromatogramList"}:
                    break
                if tag == "run":
                    metadata.started_at = element.get("startTimeStamp")
                continue
            if tag == "referenceableParamGroup":
                groups[element.get("id", "")] = terms(element)
            elif tag == "sourceFile":
                name = element.get("name")
                if name:
                    metadata.source_files.append(name)
            elif tag == "instrumentConfiguration":
                for term in terms(element):
                    acc, name = term.get("accession", ""), term.get("name", "")
                    if acc == "MS:1000529":
                        metadata.serial_numbers.append(term.get("value", ""))
                    elif acc.startswith("MS:") and not term.get("value"):
                        metadata.instruments.append(CVTerm(acc, name))
                for child in element.iter():
                    target = {"analyzer": metadata.analyzers, "source": metadata.ionization}.get(
                        _tag(child)
                    )
                    if target is not None:
                        for term in terms(child):
                            if not term.get(
                                "value",
                            ) and term.get("accession", "").startswith(
                                "MS:",
                            ):
                                target.append(CVTerm(term["accession"], term.get("name", "")))
                element.clear()
    for key in ("instruments", "analyzers", "ionization", "serial_numbers", "source_files"):
        setattr(metadata, key, list(dict.fromkeys(getattr(metadata, key))))
    return metadata


def _enum_int(value: Any) -> int:
    """Support the integer enums in 3.5 and Python Enum bindings in 3.6+."""
    return int(value.value if hasattr(value, "value") else value)


def _enum_names(enum: Any) -> dict[int, str]:
    result = {}
    for name in dir(enum):
        if name.startswith("_") or name.startswith("SIZE_OF"):
            continue
        try:
            result[_enum_int(getattr(enum, name))] = name
        except (AttributeError, ValueError, TypeError):
            continue
    return result


def _meta_number(obj: Any, keys: tuple[str, ...]) -> float | None:
    for key in keys:
        if obj.metaValueExists(key):
            try:
                value = float(obj.getMetaValue(key))
            except (TypeError, ValueError):
                continue
            if math.isfinite(value):
                return value
    return None


class _Consumer:
    """The camelCase methods are pyOpenMS's required callback interface."""

    def __init__(self, sink: SpectrumSink, oms: Any, estimate_peak_type: bool = False) -> None:
        self.sink = sink
        self.activations = _enum_names(oms.Precursor.ActivationMethod)
        self.polarities = _enum_names(oms.IonSource.Polarity)
        self.types = _enum_names(oms.SpectrumSettings.SpectrumType)
        self.faims_unit = getattr(oms.DriftTimeUnit, "FAIMS_COMPENSATION_VOLTAGE", None)
        self.estimator = oms.PeakTypeEstimator() if estimate_peak_type else None

    def setExperimentalSettings(self, settings: Any) -> None:
        pass

    def setExpectedSize(self, spectra: int, chromatograms: int) -> None:
        pass

    def consumeChromatogram(self, chromatogram: Any) -> None:
        rt, _ = chromatogram.get_peaks()
        self.sink.consume_chromatogram(
            np.asarray(rt, dtype=float),
            _enum_int(chromatogram.getChromatogramType()),
        )

    def consumeSpectrum(self, spectrum: Any) -> None:
        mz, intensity = spectrum.get_peaks()
        precursors = []
        for precursor in spectrum.getPrecursors():
            methods = tuple(
                ACTIVATION[name]
                for method in sorted(precursor.getActivationMethods(), key=_enum_int)
                if (name := self.activations.get(_enum_int(method), "").upper()) in ACTIVATION
            )
            normalized = _meta_number(
                precursor,
                ("MS:1000138", "percent collision energy", "normalized collision energy"),
            )
            energy = _meta_number(precursor, ("MS:1000045", "collision energy"))
            # OpenMS's activation-energy field corresponds to collision energy in eV.
            # Zero is also its unset default, so do not invent an observation for it.
            if energy is None and precursor.getActivationEnergy() > 0:
                energy = float(precursor.getActivationEnergy())
            collision = (normalized, "NCE") if normalized is not None else (
                (energy, "eV") if energy is not None else None
            )
            width = float(
                precursor.getIsolationWindowLowerOffset() + precursor.getIsolationWindowUpperOffset(),
            )
            precursors.append(Precursor(
                mz=float(precursor.getMZ()), charge=int(precursor.getCharge()),
                intensity=float(precursor.getIntensity()),
                isolation_width=width if math.isfinite(width) and width > 0 else None,
                activation=methods, collision_energy=collision,
            ))
        settings = spectrum.getInstrumentSettings()
        polarity = self.polarities.get(_enum_int(settings.getPolarity()), "unknown").lower()
        if polarity not in {"positive", "negative"}:
            polarity = "unknown"
        representation = self.types.get(_enum_int(spectrum.getType()), "unknown").lower()
        if representation not in {"centroid", "profile"}:
            representation = "unknown"
        estimated = None
        if self.estimator is not None and mz.size > 10:
            estimated = self.types.get(
                _enum_int(self.estimator.estimateType(spectrum)),
                "unknown",
            ).lower()
        faims = _meta_number(spectrum, ("MS:1001581", "FAIMS compensation voltage"))
        if faims is None and self.faims_unit is not None:
            if _enum_int(spectrum.getDriftTimeUnit()) == _enum_int(self.faims_unit):
                faims = float(spectrum.getDriftTime())
        self.sink.consume_spectrum(Spectrum(
            ms_level=int(spectrum.getMSLevel()), rt=float(spectrum.getRT()),
            mz=np.asarray(mz, dtype=float), intensity=np.asarray(intensity, dtype=float),
            precursors=tuple(precursors), polarity=polarity, representation=representation,
            scan_windows=tuple((float(w.begin), float(w.end)) for w in settings.getScanWindows()),
            faims_cv=faims, native_id=str(spectrum.getNativeID()),
            estimated_representation=estimated,
        ))


class PyOpenMSReader:
    """Read indexed and non-indexed mzML without holding MSExperiment in memory."""

    def __init__(self, *, estimate_peak_type: bool = False) -> None:
        try:
            import pyopenms
        except ImportError as exc:
            raise RuntimeError("pyOpenMS is required for mzML input; run `uv sync` first.") from exc
        self._oms = pyopenms
        self.engine_version = str(pyopenms.__version__)
        self.estimate_peak_type = estimate_peak_type

    def read(self, path: Path, sink: SpectrumSink) -> RunMetadata:
        path = path.resolve(strict=True)
        name = path.name.lower()
        if not name.endswith((".mzml", ".mzml.gz")):
            raise ValueError(f"Unsupported input {path.name!r}; convert vendor data to mzML first.")
        if name.endswith(".gz"):
            # Use disk, not RAM, for decompression; parser support is then identical.
            with tempfile.TemporaryDirectory(prefix="prideqc-gzip-") as folder:
                expanded = Path(folder) / path.stem
                with gzip.open(path, "rb") as source, expanded.open("wb") as destination:
                    shutil.copyfileobj(source, destination, length=1024 * 1024)
                return self._read_mzml(expanded, sink)
        return self._read_mzml(path, sink)

    def _read_mzml(self, path: Path, sink: SpectrumSink) -> RunMetadata:
        metadata = read_header(path)
        self._oms.MzMLFile().transform(
            str(path),
            _Consumer(sink, self._oms, self.estimate_peak_type),
        )
        return metadata
