# QC workflow

The central purpose of prideQC is to compute **quality-control metrics and technical evidence directly from raw mass-spectrometry data**.

```{figure} _static/figures/figure2_qc_metrics.svg
:alt: Synthetic LC-MS run with TIC, scan timing, base peaks, precursor charges, isolation windows, and repeated-target metrics.
:width: 100%
:class: doc-figure

**Synthetic run, real metric calculations.** A synthetic 90-minute DDA run is reduced by prideQC's `RunSummary` and `QCMetricCalculator` into TIC area/CV, base-peak intensity, acquisition timing, precursor-charge and isolation-window summaries, and repeated-target statistics.
```

## Inputs

prideQC can operate on:

- local mzML files;
- selected files downloaded from a PRIDE accession;
- supported vendor RAW inputs through available native readers or configured conversion tools.

An SDRF is optional for run-level QC. It becomes important when the goal is to map the resulting evidence back onto study metadata.

## Evidence families

### General QC metrics

prideQC summarizes run- and spectrum-level properties such as scan counts, peak statistics, retention-time behavior, chromatogram information, precursor characteristics, isolation-window properties, polarity, peak representation, and other technical descriptors.

### Acquisition evidence

The analysis can retain evidence about acquisition characteristics and diagnostic ions without forcing that evidence into a categorical metadata decision.

### Measurement precision

Repeated observations are used to estimate precursor and fragment mass-error precision. Supported estimates can be converted into conservative search-tolerance suggestions. See [Measurement precision & mass shifts](measurement-evidence.md).

### Recurrent neutral-mass shifts

Related spectra can reveal repeated neutral-mass-shift families. prideQC can report mass-compatible OpenMS / UniMod candidates while retaining ambiguity instead of treating mass compatibility as proof of modification identity.

## Structured outputs

prideQC writes evidence as machine-readable artifacts so that later stages do not need to scrape console logs.

For example:

```text
results/
├── metrics.tsv                    # study-wide flattened QC table
├── annotations.tsv                # study-wide evidence / annotation table
├── mass-shifts.tsv                # aggregate mass-shift evidence, when enabled
├── runA.mzML.mzQC                 # mzQC-compatible per-file report
├── runA.mzML.summary.json         # complete structured per-file summary
├── runA.mzML.mass-shifts.tsv      # per-file mass-shift evidence
└── manifest.json                  # run provenance and success/failure accounting
```

## From structured QC to visual reporting

The mzQC-compatible output can be aggregated into an interactive report using pmultiqc / MultiQC. See [Reporting](reporting.md).

## From QC to SDRF refinement

The run-level evidence remains separate from metadata write-back. A later cohort stage asks whether the same evidence is sufficiently supported across the relevant runs before proposing SDRF changes.

This separation is deliberate:

- **measurement** happens at file/run level;
- **study interpretation** happens at cohort level;
- **metadata changes** happen only after policy and validation gates.
