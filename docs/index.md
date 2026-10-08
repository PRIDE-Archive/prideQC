# prideQC documentation

<div class="hero-title">Quality control from raw mass-spectrometry data, with evidence-based SDRF refinement</div>
<div class="hero-subtitle">
prideQC computes QC metrics directly from mass-spectrometry data, produces structured and visual reports,
and then uses supported study-level evidence to improve technical SDRF metadata while preserving provenance.
</div>

```{figure} _static/figures/figure1_overview.svg
:alt: A raw mass-spectrometry input enters per-file QC, which creates structured evidence used independently by reporting and gated SDRF refinement.
:width: 100%
:class: doc-figure

prideQC analyzes raw data into per-file QC evidence. The resulting mzQC files can be summarized with pmultiqc, while structured run summaries feed cohort-level SDRF refinement. Uncertain proposals can abstain.
```

::::{grid} 1 2 2 2
:gutter: 3
:class-container: note-grid

:::{grid-item-card} Getting started
:link: getting-started
:link-type: doc
Install prideQC, analyze raw data, inspect outputs, and refine an SDRF.
:::

:::{grid-item-card} QC workflow
:link: qc-workflow
:link-type: doc
See what prideQC measures from spectra and how the evidence is represented.
:::

:::{grid-item-card} Measurement precision & mass shifts
:link: measurement-evidence
:link-type: doc
Understand mass-error estimation and recurrent neutral-mass-shift scouting with synthetic examples.
:::

:::{grid-item-card} Reporting
:link: reporting
:link-type: doc
Generate interactive study-level QC reports with pmultiqc / MultiQC.
:::

:::{grid-item-card} SDRF refinement
:link: sdrf-refinement
:link-type: doc
Learn how run-level QC evidence is synthesized into conservative, auditable metadata proposals.
:::

:::{grid-item-card} CLI reference
:link: cli-reference
:link-type: doc
Reference the main prideQC commands and options.
:::
::::

## What prideQC is designed to do

prideQC connects two tasks that are usually handled separately:

1. **Compute reproducible QC evidence from raw mass-spectrometry data.**
2. **Use sufficiently supported evidence to refine technical SDRF metadata.**

The first task is always primary. Metadata refinement is downstream of the measured QC evidence and can abstain when the evidence is incomplete, ambiguous, or inconsistent.

## Core capabilities

::::{grid} 1 1 2 2
:gutter: 3

:::{grid-item-card} Run-level QC
- scan and chromatogram metrics
- acquisition characteristics
- peak-type and diagnostic-ion evidence
- structured per-file summaries
:::

:::{grid-item-card} Measurement precision
- precursor mass-error precision
- fragment mass-error precision
- conservative tolerance suggestions
- explicit support and diagnostics
:::

:::{grid-item-card} Mass-shift evidence
- recurrent related-spectrum neutral-mass shifts
- OpenMS / UniMod mass-compatible candidates
- ambiguity and support retained explicitly
:::

:::{grid-item-card} Study-level refinement
- cohort synthesis across files
- confidence-gated SDRF proposals
- validation, provenance, and audit artifacts
- optional adjudication for unresolved cases
:::
::::

## Documentation map

```{toctree}
:maxdepth: 2
:hidden:

getting-started
qc-workflow
measurement-evidence
reporting
sdrf-refinement
cli-reference
```

## External resources

- prideQC source: <https://github.com/PRIDE-Archive/prideQC>
- PRIDE Archive: <https://www.ebi.ac.uk/pride/>
- OpenMS: <https://www.openms.de/>
- pyOpenMS: <https://pyopenms.readthedocs.io/>
- pmultiqc: <https://github.com/bigbio/pmultiqc>
- MultiQC: <https://seqera.io/multiqc/>
- SDRF-Pipelines: <https://github.com/bigbio/sdrf-pipelines>
