#!/usr/bin/env python3
"""Generate the public metrics reference directly from the versioned registry."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REGISTRY = ROOT / "docs" / "metrics.tsv"
DEST = ROOT / "docs" / "qc-metrics-reference.md"

ANNOTATIONS = [
    ("acquisition_method", "DDA/DIA/PRM/SRM run classification", "Explicit header CV or DDA/DIA inference; targeted scans abstain"),
    ("acquisition_method_evidence", "PRIDE CV term, support and classifier", "Structured mzQC; observed or inferred"),
    ("acquisition_method_conflict", "Contradictory RAW-header and spectral evidence", "No automatic classification/writeback"),
    ("labeling_strategy", "Verified sample labelling strategy", "Unavailable from RAW alone; no SDRF writeback"),
    ("labeling_strategy_candidates", "TMT/iTRAQ reporter evidence", "Optional, exploratory, non-exclusive; structured mzQC"),
    ("ionization", "Ion source PSI-MS CV terms", "Observed when present in instrument configuration"),
    ("chromatography_type", "LC separation type", "Unavailable unless external method metadata is supplied"),
    ("instrument", "Instrument PSI-MS CV", "Observed where supplied"),
    ("mass_analyzers", "Analyzer PSI-MS CV", "Observed where supplied"),
    ("dissociation_ms2", "MS2 dissociation methods", "Observed when scan metadata are complete"),
    ("collision_energy_ms2", "Recorded collision energies", "Observed where present"),
    ("tmt_reporter_evidence", "TMT centroid MS2/MS3 channel evidence", "Optional experimental collector"),
    ("tmt_profile_reporter_evidence", "TMT profile MS2/MS3 evidence", "Optional experimental collector"),
    ("signature_TMT_family", "Diagnostic TMT reporter-ion screen", "Optional, candidate-only"),
    ("signature_iTRAQ_family", "Diagnostic iTRAQ reporter-ion screen", "Optional, candidate-only"),
    ("estimated_precursor_mass_error_ppm", "Precursor mass precision", "Optional estimator"),
    ("estimated_fragment_mass_error_ppm", "Fragment mass precision", "Optional estimator"),
    ("putative_modification_mass_shifts", "Neutral mass-shift hypotheses", "Optional; not peptide IDs"),
]


def category(key: str) -> str:
    if "Charge" in key or "Precursor" in key or "Isolation" in key:
        return "Precursors and isolation"
    if key.startswith("AcquisitionCycle"):
        return "Acquisition cycles"
    if key.startswith("FAIMS"):
        return "Ion mobility / FAIMS"
    if "Chromatogram" in key or "Chromatography" in key:
        return "Chromatography"
    if "TIC" in key or "BasePeak" in key or "Signal" in key:
        return "Ion intensity and TIC"
    if "RT_" in key or "RtRange" in key or "Frequency" in key or "ScanRate" in key:
        return "Scan timing"
    if "Spectra" in key or "MSLevels" in key or "Peak" in key or "Polarity" in key:
        return "Spectra and peak summaries"
    return "Other run-level QC"


def render() -> str:
    with REGISTRY.open(encoding="utf-8", newline="") as handle:
        entries = list(csv.DictReader(handle, delimiter="\t"))
    keys = [row["key"] for row in entries]
    if len(set(keys)) != len(keys):
        raise ValueError("Duplicate keys in metrics.tsv")
    for row in entries:
        if not all(row.get(k) for k in ("key", "accession", "shape", "description")):
            raise ValueError(f"Incomplete registry entry: {row}")
    lines = [
        "# QC metrics and experimental annotations", "",
        "This page is generated from the tracked `docs/metrics.tsv` registry. The",
        "registered metrics below describe **QC measurements**; method/label",
        "annotations are separate evidence fields and are not necessarily quantitative.",
        "A registered metric can be absent from a particular mzQC when measurement",
        "is unavailable; optional collectors require their respective CLI flags.", "",
        f"**Registered core QC metric definitions: {len(entries)}.**", "",
        "## Experimental and technical annotation fields", "",
        "| Field | Information | Evidence policy |", "| --- | --- | --- |",
    ]
    for field, meaning, evidence in ANNOTATIONS:
        lines.append(f"| `{field}` | {meaning} | {evidence} |")
    lines.extend(["", "The PRIDE ontology groups acquisition types under `PRIDE:0000659`",
                  "and sample labels under `PRIDE:0000514`.", "",
                  "Acquisition classifications use `PRIDE:0000450` (DIA) and",
                  "`PRIDE:0000627` (DDA) only with sufficient spectral evidence.",
                  "PRM (`PRIDE:0000629`) and SRM/MRM (`PRIDE:0000630`) are not",
                  "automatically called from narrow-isolation scans. Exact PRIDE acquisition",
                  "terms recorded in mzML fileContent may establish PRM or SRM/MRM,",
                  "unless they conflict with spectral evidence. TMT (`PRIDE:0000314`)",
                  "and iTRAQ (`PRIDE:0000313`) may be **candidates**, never sample-label",
                  "assignments. SILAC (`PRIDE:0000315`) cannot be inferred from",
                  "ordinary precursor charge/isotope information alone. No RAW-only",
                  "absence of signals is interpreted as label-free.", "",
                  "## Registered quantitative QC metrics", ""])
    groups: dict[str, list[dict[str, str]]] = {}
    for entry in entries:
        groups.setdefault(category(entry["key"]), []).append(entry)
    for name, rows in sorted(groups.items()):
        lines.extend([f"### {name} ({len(rows)})", "", "| Key | CV accession | Shape | Definition |",
                      "| --- | --- | --- | --- |"])
        for row in rows:
            desc = row["description"].replace("|", "\\|").replace("\n", " ")
            lines.append(f"| `{row['key']}` | `{row['accession']}` | {row['shape']} | {desc} |")
        lines.append("")
    lines.extend(["## Provenance and caveats", "",
                  "Metrics are calculated per MS run and merged without recalculation.",
                  "`NumberOfSpectra_MS1`, `NumberOfSpectra_MS2`,",
                  "`MS2_PrecursorCharge_Fractions`, `ChromatographyDuration`,",
                  "and precursor isolation-width statistics are already part of the registry.",
                  "Ion-source CV metadata are observed only if the reader supplies it.",
                  "RT/duration do not reveal LC separation chemistry. Label/channel",
                  "annotations require verified submitter/SDRF evidence for writeback.", "",
                  "See also [QC workflow](qc-workflow.md) and [measurement evidence](measurement-evidence.md).", ""])
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="Fail if the generated doc is stale")
    args = parser.parse_args()
    content = render()
    if args.check:
        if not DEST.exists() or DEST.read_text(encoding="utf-8") != content:
            print(f"METRICS_REFERENCE_STALE={DEST}")
            return 1
        print("METRICS_REFERENCE=PASS")
    else:
        DEST.write_text(content, encoding="utf-8")
        print(f"METRICS_REFERENCE_WRITTEN={DEST}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
