#!/usr/bin/env python3
"""Report the pyOpenMS capabilities used by the v22 mass-shift scout.

Run this inside the exact prideQC image/SIF. It is intentionally read-only and
prints JSON so cluster logs can be archived as provenance for the implementation
benchmark.
"""

from __future__ import annotations

import json
import time
from typing import Any


def _text(value: Any) -> str:
    if isinstance(value, bytes):
        return value.decode(errors="replace")
    return str(value)


def main() -> int:
    import pyopenms as oms

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
        name = _text(modification.getFullName()).strip() or _text(modification.getId()).strip()
        if len(examples) < 5:
            examples.append(
                {
                    "accession": accession,
                    "name": name,
                    "delta_mono_mass_da": float(modification.getDiffMonoMass()),
                }
            )

    spectrum_alignment = getattr(oms, "SpectrumAlignment", None)
    spectrum_alignment_score = getattr(oms, "SpectrumAlignmentScore", None)
    output = {
        "pyopenms_version": getattr(oms, "__version__", None),
        "openms_proton_mass_u": float(oms.Constants.PROTON_MASS_U),
        "openms_c13_c12_massdiff_u": float(oms.Constants.C13C12_MASSDIFF_U),
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
        "spectrum_alignment_available": spectrum_alignment is not None,
        "spectrum_alignment_score_available": spectrum_alignment_score is not None,
        "examples": examples,
        "elapsed_seconds": time.perf_counter() - started,
    }
    print(json.dumps(output, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
