"""pyOpenMS adapters for mzML and optional native vendor readers."""

from __future__ import annotations

import gzip
import math
import shutil
import tempfile
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any
from zipfile import ZipFile

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


class VendorReaderUnavailable(RuntimeError):
    """Raised when the installed pyOpenMS build has no native vendor reader."""


def vendor_format(path: Path) -> str | None:
    """Return the supported vendor format represented by *path*, if any."""
    name = path.name.casefold()
    if name.endswith(".raw"):
        return "Thermo RAW"
    if name.endswith(".d.zip") or name.endswith(".d"):
        return "Bruker TDF (.d)"
    return None


def _loader_for(oms: Any, format_name: str) -> Any | None:
    """Return the concrete native reader class exposed by pyOpenMS."""
    class_name = "ThermoRawFile" if format_name == "Thermo RAW" else "BrukerTimsFile"
    if not hasattr(oms, class_name):
        return None
    return getattr(oms, class_name, None)


def direct_vendor_support(oms: Any, path: Path) -> bool:
    """Feature-detect whether *oms* can directly read this vendor path."""
    format_name = vendor_format(path)
    return format_name is not None and _loader_for(oms, format_name) is not None


def _vendor_error(path: Path, version: str) -> VendorReaderUnavailable:
    format_name = vendor_format(path) or "vendor"
    return VendorReaderUnavailable(
        f"Direct reading of {format_name} input {path.name!r} is unavailable in pyOpenMS "
        f"{version}. Direct reading requires a pyOpenMS build exposing "
        "ThermoRawFile or BrukerTimsFile. Current development builds are available from "
        "https://pypi.openms.de/simple/pyopenms/. Alternatively upgrade pyOpenMS or "
        "provide --converter thermorawfileparser/msconvert."
    )


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


def _text_value(obj: Any, names: tuple[str, ...]) -> str | None:
    """Read a string-valued OpenMS property when the binding exposes it."""
    for name in names:
        getter = getattr(obj, name, None)
        if not callable(getter):
            continue
        try:
            value = getter()
        except (AttributeError, TypeError, RuntimeError):
            continue
        if isinstance(value, bytes):
            value = value.decode(errors="replace")
        text = str(value).strip()
        if text:
            return text
    return None


def _metadata_from_experiment(experiment: Any, source: Path) -> RunMetadata:
    """Copy metadata exposed by MSExperiment without guessing CV accessions."""
    metadata = RunMetadata(source_files=[source.name])
    settings_getter = getattr(experiment, "getExperimentalSettings", None)
    settings = settings_getter() if callable(settings_getter) else experiment

    source_getter = getattr(settings, "getSourceFiles", None)
    if callable(source_getter):
        try:
            for source_file in source_getter():
                name = _text_value(source_file, ("getNameOfFile", "getName", "getPathToFile"))
                if name:
                    metadata.source_files.append(name)
        except (AttributeError, TypeError, RuntimeError):
            pass

    date_getter = getattr(settings, "getDateTime", None)
    if callable(date_getter):
        try:
            value = date_getter()
            if value:
                metadata.started_at = str(value)
        except (AttributeError, TypeError, RuntimeError):
            pass

    instrument = None
    instrument_getter = getattr(settings, "getInstrument", None)
    if callable(instrument_getter):
        try:
            instrument = instrument_getter()
        except (AttributeError, TypeError, RuntimeError):
            instrument = None
        if instrument is not None:
            for key, names in (
                ("name", ("getName",)),
                ("vendor", ("getVendor",)),
                ("model", ("getModel",)),
                ("software", ("getSoftware",)),
            ):
                value = _text_value(instrument, names)
                if value:
                    metadata.instrument_details[key] = value
            serial = _text_value(instrument, ("getSerialNumber",))
            if serial is None:
                exists = getattr(instrument, "metaValueExists", None)
                getter = getattr(instrument, "getMetaValue", None)
                if callable(exists) and callable(getter):
                    for key in ("MS:1000529", "serial number"):
                        try:
                            if exists(key):
                                serial = str(getter(key))
                                break
                        except (AttributeError, TypeError, RuntimeError):
                            continue
            if serial:
                metadata.serial_numbers.append(serial)

    for key, getter_names in (
        ("analyzers", ("getMassAnalyzers",)),
        ("ionization", ("getIonSources",)),
    ):
        getter = getattr(instrument if instrument is not None else settings, getter_names[0], None)
        if not callable(getter):
            continue
        try:
            values = getter()
        except (AttributeError, TypeError, RuntimeError):
            continue
        target = getattr(metadata, key)
        for value in values:
            accession = _text_value(value, ("getAccession", "accession"))
            name = _text_value(value, ("getName", "name"))
            if accession and name and accession.startswith(("MS:", "UO:", "PRIDE:")):
                target.append(CVTerm(accession, name))

    for key in ("instruments", "analyzers", "ionization", "serial_numbers", "source_files"):
        setattr(metadata, key, list(dict.fromkeys(getattr(metadata, key))))
    return metadata


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
    """Read mzML and, when available, Thermo/Bruker vendor files."""

    def __init__(self, *, estimate_peak_type: bool = False, oms: Any | None = None) -> None:
        if oms is None:
            try:
                import pyopenms
            except ImportError as exc:
                raise RuntimeError("pyOpenMS is required for mzML input; run `uv sync` first.") from exc
            oms = pyopenms
        # pyOpenMS exposes a dynamic binding surface without complete typing
        # stubs (notably MzMLFile.transform and __version__). Keep the native
        # boundary typed as Any while the rest of the reader remains checked.
        self._oms: Any = oms
        self.engine_version = str(getattr(oms, "__version__", "unknown"))
        self.estimate_peak_type = estimate_peak_type

    def supports_direct(self, path: Path) -> bool:
        """Return whether a native reader is exposed for this vendor path."""
        return direct_vendor_support(self._oms, path)

    def unavailable_error(self, path: Path) -> VendorReaderUnavailable:
        """Build the actionable fallback error for an unsupported vendor path."""
        return _vendor_error(path, self.engine_version)

    def read(self, path: Path, sink: SpectrumSink) -> RunMetadata:
        path = path.resolve(strict=True)
        name = path.name.casefold()
        if name.endswith((".mzml", ".mzml.gz")):
            if name.endswith(".gz"):
                # Use disk, not RAM, for decompression; parser support is then identical.
                with tempfile.TemporaryDirectory(prefix="prideqc-gzip-") as folder:
                    expanded = Path(folder) / path.stem
                    with gzip.open(path, "rb") as source, expanded.open("wb") as destination:
                        shutil.copyfileobj(source, destination, length=1024 * 1024)
                    return self._read_mzml(expanded, sink)
            return self._read_mzml(path, sink)

        if vendor_format(path) is None:
            raise ValueError(f"Unsupported input {path.name!r}; convert vendor data to mzML first.")
        if not self.supports_direct(path):
            raise self.unavailable_error(path)
        return self._read_vendor(path, sink)

    def _read_mzml(self, path: Path, sink: SpectrumSink) -> RunMetadata:
        metadata = read_header(path)
        self._oms.MzMLFile().transform(
            str(path),
            _Consumer(sink, self._oms, self.estimate_peak_type),
        )
        return metadata

    def _load_vendor(self, path: Path, format_name: str) -> Any:
        """Use the public pyOpenMS vendor APIs without signature guessing."""
        if format_name == "Thermo RAW":
            experiment = self._oms.MSExperiment()
            self._oms.ThermoRawFile().load(str(path), experiment)
            return experiment
        return self._oms.BrukerTimsFile().load(str(path))

    def _consume_experiment(self, experiment: Any, sink: SpectrumSink) -> None:
        consumer = _Consumer(sink, self._oms, self.estimate_peak_type)
        spectra_getter = getattr(experiment, "getSpectra", None) or getattr(experiment, "get_spectra", None)
        chromatograms_getter = getattr(experiment, "getChromatograms", None) or getattr(
            experiment, "get_chromatograms", None
        )
        if callable(spectra_getter):
            spectra = spectra_getter()
        elif callable(getattr(experiment, "getNrSpectra", None)) and callable(
            getattr(experiment, "getSpectrum", None),
        ):
            spectra = (experiment.getSpectrum(index) for index in range(experiment.getNrSpectra()))
        else:
            spectra = experiment
        for spectrum in spectra:
            consumer.consumeSpectrum(spectrum)
        if callable(chromatograms_getter):
            for chromatogram in chromatograms_getter():
                consumer.consumeChromatogram(chromatogram)

    @staticmethod
    def _extract_d_archive(path: Path) -> tuple[Any, Path]:
        temporary = tempfile.TemporaryDirectory(prefix="prideqc-bruker-")
        root = Path(temporary.name)
        try:
            with ZipFile(path) as archive:
                for member in archive.infolist():
                    target = (root / member.filename).resolve()
                    if not target.is_relative_to(root):
                        raise ValueError("Bruker .d.zip contains an unsafe archive path.")
                archive.extractall(root)
            directories = [candidate for candidate in root.rglob("*")
                           if candidate.is_dir() and candidate.name.casefold().endswith(".d")]
            if len(directories) != 1:
                raise ValueError("Bruker .d.zip must contain exactly one .d directory.")
            return temporary, directories[0]
        except Exception:
            temporary.cleanup()
            raise

    def _read_vendor(self, path: Path, sink: SpectrumSink) -> RunMetadata:
        format_name = vendor_format(path)
        assert format_name is not None
        if _loader_for(self._oms, format_name) is None:
            raise self.unavailable_error(path)
        temporary = None
        try:
            if format_name == "Bruker TDF (.d)" and path.name.casefold().endswith(".d.zip"):
                temporary, extracted = self._extract_d_archive(path)
                experiment = self._load_vendor(extracted, format_name)
            else:
                experiment = self._load_vendor(path, format_name)
            # Consume before removing a temporary archive extraction: some
            # reader implementations expose lazy/on-disc experiment objects.
            self._consume_experiment(experiment, sink)
        finally:
            if temporary is not None:
                temporary.cleanup()
        # Vendor readers do not expose mzML's bounded XML header. Preserve all
        # metadata that the loaded experiment provides, without inventing CVs.
        return _metadata_from_experiment(experiment, path)
