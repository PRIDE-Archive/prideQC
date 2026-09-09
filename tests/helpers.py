"""Deterministic fixtures independent of the native reader."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from prideqc.models import CVTerm, Precursor, RunMetadata, Spectrum, SpectrumSink

HCD = CVTerm("MS:1000422", "beam-type collision-induced dissociation")
FUSION = CVTerm("MS:1002416", "Orbitrap Fusion")


def spectrum(rt: float, intensities: list[float], level: int = 1,
             charge: int | None = None, precursor_intensity: float = 100.0,
             width: float = 2.0, mz: list[float] | None = None,
             representation: str = "centroid") -> Spectrum:
    peaks = np.array(mz if mz is not None else [100.0 + i for i in range(len(intensities))])
    precursors = () if charge is None else (Precursor(
        500,
        charge,
        precursor_intensity,
        width,
        (HCD,),
        (30, "eV"),
    ),)
    return Spectrum(level, rt, peaks.astype(float), np.array(intensities, dtype=float),
                    precursors, "positive", representation)


class MemoryReader:
    engine_version = "test-reader"

    def __init__(self, spectra: list[Spectrum], metadata: RunMetadata | None = None) -> None:
        self.spectra = spectra
        self.metadata = metadata or RunMetadata(instruments=[FUSION])
        self.read_count = 0

    def read(self, path: Path, sink: SpectrumSink) -> RunMetadata:
        self.read_count += 1
        for item in self.spectra:
            sink.consume_spectrum(item)
        return self.metadata
