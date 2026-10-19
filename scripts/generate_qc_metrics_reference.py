#!/usr/bin/env python3
"""Generate the public metrics reference directly from the versioned registry."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REGISTRY = ROOT / "docs" / "metrics.tsv"
DEST = ROOT / "docs" / "qc-metrics-reference.md"
PRESENTATION = ROOT / "docs" / "metric-presentation.tsv"

TIERS = ("core", "diagnostic", "extended")
APPLICABILITY = ("all", "ms1", "ms2", "ms1_ms2", "cycles", "faims", "chromatograms")
CORE_COUNT = 8

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


def load_catalog(
    registry: Path = REGISTRY, presentation: Path = PRESENTATION,
) -> list[tuple[dict[str, str], dict[str, str]]]:
    """Validate an exact, one-to-one display policy over the scientific registry.

    The display metadata does not define CV terms, compute values or modify mzQC.
    """
    with registry.open(encoding="utf-8", newline="") as handle:
        entries = list(csv.DictReader(handle, delimiter="\t"))
    with presentation.open(encoding="utf-8", newline="") as handle:
        policies = list(csv.DictReader(handle, delimiter="\t"))

    keys = [row["key"] for row in entries]
    display_keys = [row["key"] for row in policies]
    if len(set(keys)) != len(keys) or len(set(display_keys)) != len(display_keys):
        raise ValueError("Duplicate metric keys in registry or presentation policy")
    if keys != display_keys:
        raise ValueError("Metric presentation must match the registry exactly, in order")
    for metric, policy in zip(entries, policies, strict=True):
        if not all(metric.get(k) for k in ("key", "accession", "shape", "description")):
            raise ValueError(f"Incomplete metric definition: {metric}")
        if policy.get("tier") not in TIERS:
            raise ValueError(f"Unrecognized presentation tier for {metric['key']}")
        if policy.get("applicability") not in APPLICABILITY:
            raise ValueError(f"Unrecognized applicability for {metric['key']}")
        if policy["tier"] == "core" and not all(policy.get(k) for k in ("label", "caveat")):
            raise ValueError(f"Core metric lacks display label or caveat: {metric['key']}")
    if sum(p["tier"] == "core" for p in policies) != CORE_COUNT:
        raise ValueError(f"Expected exactly {CORE_COUNT} curated core metrics")
    return list(zip(entries, policies, strict=True))


def _cell(value: str) -> str:
    return value.replace("|", "\\|").replace("\n", " ")


def render() -> str:
    catalog = load_catalog()
    counts = {tier: sum(p["tier"] == tier for _, p in catalog) for tier in TIERS}
    lines = [
        "# QC metrics and experimental annotations", "",
        "This page is generated from the tracked `docs/metrics.tsv` measurement",
        "registry and `docs/metric-presentation.tsv` display policy. The latter",
        "selects which existing mzQC metrics are prominently displayed; it does",
        "**not** change acquisition, calculation, CV identity, serialization or merging.",
        "Experimental method and label annotations are recorded separately.", "",
        f"**Registered QC metric definitions: {len(catalog)}.**",
        f"Presentation tiers: **{counts['core']} core**, **{counts['diagnostic']} diagnostic**,",
        f"and **{counts['extended']} extended**. A registration is not a promise",
        "that a value exists for every run. Metrics absent or undefined in the source",
        "mzQC must remain absent/undefined; do not fill missing values with zero.", "",
        "## Recommended PRIDE website overview", "",
        f"These {CORE_COUNT} **core display candidates** cover run structure and acquisition.",
        "Show a value or visualization only when the run has the corresponding",
        "stored mzQC metric. For SRM/MRM or other MS1-free acquisitions, MS1 and",
        "metrics may be inapplicable; this is not QC failure.", "",
        "The overview uses readable names; the full registry keys, CV accessions,",
        "value shapes and definitions remain available in the sections below.", "",
        "| Core metric | Applies to | Interpretation / limitation |",
        "| --- | --- | --- |",
    ]
    for _metric, policy in catalog:
        if policy["tier"] == "core":
            lines.append(
                f"| {_cell(policy['label'])} | {policy['applicability']} | "
                f"{_cell(policy['caveat'])} |"
            )
    lines.extend([
        "", "Counts include empty scans where stated by the registry. Charge",
        "statistics and DIA windows are reported only when supported. Quantile",
        "and tuple metrics remain multivalued; the overview must not silently collapse",
        "them into an undocumented scalar or recompute them from mzQC.", "",
        "Acquisition method (DDA/DIA/PRM/SRM), labelling candidates, ionization,",
        "instrument and LC annotations belong in the experimental-context section,",
        "not in a numeric QC score. No universal composite QC score is defined.", "",
        "The number of distinct MS levels is no longer a QC metric: `MS:1000511`",
        "is a per-spectrum level indicator, not a run-level count. Observed FAIMS",
        "voltages are recorded individually as `MS:1001581` file metadata (volts).",
        "The ramp endpoint terms `MS:1003450` and `MS:1003451` are not emitted",
        "for stepped CV acquisitions. MS1-delimited cycle summaries remain internal",
        "evidence for the conservative DDA/DIA classifier, not published QC metrics.", "",
        "## Experimental and technical annotation fields", "",
        "| Field | Information | Evidence policy |", "| --- | --- | --- |",
    ])
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
                  "## Registered quantitative QC metrics", "",
                  "Tier and applicability are editorial display metadata, not additional",
                  "scientific measurements or controlled vocabulary terms.",
                  "Browse the compact three-column indexes by category, then expand",
                  "a category's definitions to see full CV accessions, value shapes",
                  "and the original metric descriptions.", ""])
    groups: dict[str, list[tuple[dict[str, str], dict[str, str]]]] = {}
    for metric, policy in catalog:
        groups.setdefault(category(metric["key"]), []).append((metric, policy))
    for name, rows in sorted(groups.items()):
        lines.extend([
            f"### {name} ({len(rows)})", "",
            "| Registry metric | Tier | Applies to |",
            "| --- | --- | --- |",
        ])
        for metric, policy in rows:
            lines.append(
                f"| `{metric['key']}` | {policy['tier']} | "
                f"{policy['applicability']} |"
            )
        lines.extend([
            "", ":::{dropdown} Full definitions and CV accessions — " + name,
            "",
        ])
        for metric, policy in rows:
            lines.extend([
                f"- **`{metric['key']}`** — {policy['tier']}; "
                f"{metric['shape']}; applies to {policy['applicability']}",
                f"  - CV: `{metric['accession']}`",
                f"  - Definition: {metric['description']}",
                "",
            ])
        lines.extend([":::", ""])
    lines.extend(["## Provenance and caveats", "",
                  "The selected set of calculated QC metrics remains in mzQC per run",
                  "and in losslessly consolidated accession mzQC. A website or",
                  "MongoDB projection may select **stored values only**, never generate",
                  "report-only scientific measurements. Grouped spectrum counts, MS2/MS3",
                  "representation, exploratory reporter-ion evidence and denominators",
                  "must retain their original definitions. A missing measurement is",
                  "not a zero value, low score, or proof that an experimental method",
                  "was absent.", "",
                  "`NumberOfSpectra_MS1`, `NumberOfSpectra_MS2`,",
                  "`DIAIsolationWindowWidths`, `ChromatographyDuration`,",
                  "and the supported DIA window statistics remain in the scientific registry.",
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
