#!/usr/bin/env python3
"""Report pyOpenMS capabilities used by the v22.1 mass-shift scout.

Run inside the exact prideQC image/SIF. The probe is read-only apart from an
in-memory PeakPickerHiRes smoke test and prints JSON for experiment provenance.
"""

from __future__ import annotations

import json
import time
from typing import Any

import numpy as np


def _text(value: Any) -> str:
    if isinstance(value, bytes):
        return value.decode(errors="replace")
    return str(value)


def _peak_picker_smoke(oms: Any) -> dict[str, Any]:
    available = hasattr(oms, "PeakPickerHiRes") and hasattr(oms, "MSSpectrum")
    if not available:
        return {
            "available": False,
            "smoke_passed": False,
            "picked_peaks": 0,
            "error": "PeakPickerHiRes or MSSpectrum unavailable",
        }
    try:
        source = oms.MSSpectrum()
        source.set_peaks(
            (
                np.asarray([99.98, 99.99, 100.00, 100.01, 100.02], dtype=float),
                np.asarray([1.0, 8.0, 20.0, 8.0, 1.0], dtype=float),
            )
        )
        picked = oms.MSSpectrum()
        oms.PeakPickerHiRes().pick(source, picked)
        mz, _intensity = picked.get_peaks()
        return {
            "available": True,
            "smoke_passed": bool(len(mz)),
            "picked_peaks": int(len(mz)),
            "error": None,
        }
    except (AttributeError, RuntimeError, TypeError, ValueError) as exc:
        return {
            "available": True,
            "smoke_passed": False,
            "picked_peaks": 0,
            "error": f"{type(exc).__name__}: {exc}",
        }


def main() -> int:
    import pyopenms as oms
    from pyopenms.Constants import C13C12_MASSDIFF_U, PROTON_MASS_U

    started = time.perf_counter()
    database = oms.ModificationsDB()
    count = int(database.getNumberOfModifications())
    unimod = 0
    examples: list[dict[str, Any]] = []
    for index in range(count):
        modification = database.getModification(index)
        accession = _text(modification.getUniModAccession()).strip()
        if not accession.casefold().startswith("unimod:"):
            continue
        unimod += 1
        name = _text(modification.getFullName()).strip() or _text(
            modification.getId()
        ).strip()
        if len(examples) < 5:
            examples.append(
                {
                    "accession": accession,
                    "name": name,
                    "delta_mono_mass_da": float(modification.getDiffMonoMass()),
                }
            )

    output = {
        "pyopenms_version": getattr(oms, "__version__", None),
        "openms_proton_mass_u": float(PROTON_MASS_U),
        "openms_c13_c12_massdiff_u": float(C13C12_MASSDIFF_U),
        "modifications_db_entries": count,
        "unimod_bearing_entries": unimod,
        "modification_api": {
            "getNumberOfModifications": hasattr(database, "getNumberOfModifications"),
            "getModification": hasattr(database, "getModification"),
            "searchModificationsByDiffMonoMass": hasattr(
                database, "searchModificationsByDiffMonoMass"
            ),
            "getBestModificationByDiffMonoMass": hasattr(
                database, "getBestModificationByDiffMonoMass"
            ),
        },
        "peak_picker_hires": _peak_picker_smoke(oms),
        "spectrum_alignment_available": hasattr(oms, "SpectrumAlignment"),
        "spectrum_alignment_score_available": hasattr(oms, "SpectrumAlignmentScore"),
        "examples": examples,
        "elapsed_seconds": time.perf_counter() - started,
    }
    print(json.dumps(output, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
