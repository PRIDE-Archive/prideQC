# prideQC Project Handoff — 2026-09-15 — v20 historical tolerance comparison

## Project goal

prideQC processes PRIDE Archive MS data at scale to collect QC metrics and technical metadata without modifying scientific data.

The benchmark uses SDRF-annotated datasets as ground truth. Historical SDRF precursor/fragment search tolerances are preserved annotations; they are not inferred from RAW data.

## Current production architecture

- Docker is the canonical build.
- Codon uses Singularity.
- One physical PRIDE data file per Slurm array task.
- `prideqc --workers 1` inside each task.
- RAW data is staged to node-local `/tmp`.
- Compact QC/provenance outputs are persisted to NFS.
- No source checkout or compilation is required on compute nodes.

## Ground-truth manifest

Current unique accessions: 7.

Physical manifest v6:
- 663 physical analysis tasks.
- PXD000612 contains 231 physical files after authoritative PRIDE filename alias reconciliation.
- 42 PXD000612 alias records are retained as provenance and deduplicated from physical analysis.

GT controls: 10 files across the 7 accessions.

## v19 mass-error/tolerance status

v19 is the validated baseline for measurement-derived fragment tolerance candidates.

Its major behavior:
- high-resolution runs use ppm candidates from the 0.2-Da matching distribution;
- low-resolution runs adaptively select 0.5 or 1.0-Da measurement windows;
- profile spectra are never modified; profile fragment centers are ephemeral evidence;
- ambiguous/intermediate regimes abstain;
- DIA PXD018830 controls abstain rather than inventing a tolerance;
- historical search tolerances are never written or reconstructed from RAW-derived estimates.

## 55-file holdout

Run root:
`results/mass-error-v19-holdout-55-nochecksum`

The complete 55-file holdout is usable:
- 55/55 scientific analyses completed.
- 5 `success / ok`.
- 50 `warning / refined_sdrf_validation`.
- 0 scientific-analysis failures.

Task 654 (`PXD018830/TDM_M1811_009.raw`) initially failed after an FTP timeout, then succeeded with the alternative download path and the full 3,299,258,779-byte RAW.

## 55-file accession summary

Current v19 summary:
- PXD000612: 10/10 inferred, ppm, high-resolution, 0.2-Da window, median candidate ~10.21 ppm.
- PXD000759: 10/10 inferred, Da, low-resolution, 0.5/1.0-Da windows, median candidate ~0.666 Da.
- PXD001819: 10/10 inferred, Da, low-resolution, 0.5-Da window, median candidate ~0.396 Da.
- PXD003772: current holdout summary contains 1 file; inferred ppm/high-resolution/0.2-Da, ~10.47 ppm. Accession-wide generalization is not yet strong from one file.
- PXD008934: 10/10 inferred, ppm, high-resolution, 0.2-Da window, median candidate ~14.30 ppm.
- PXD017618: 4/4 inferred, ppm, high-resolution, 0.2-Da window, median candidate ~7.83 ppm; 2 high-confidence and 2 moderate-confidence results.
- PXD018830: 10/10 unavailable, correctly abstaining for DIA.

The strongest conclusion is categorical consistency: 45/45 non-DIA files inferred; 10/10 DIA files abstained; unit/regime behavior matches the curated reference state for all 55 files.

## Current v20 iteration: historical SDRF comparison

New reporting utility:
`scripts/benchmarks/compare_historical_tolerances.py`

New documentation:
`docs/mass-error-historical-comparison.md`

The utility compares the v19 candidate to the historical SDRF fragment tolerance without changing v19.

### Comparison rules

1. Read the historical value from each completed result's original `input.sdrf.tsv`.
2. Fall back to the curated GT TSV only when necessary.
3. Calculate ratio/difference/percent difference only when candidate and historical values use the same unit.
4. Mark ppm-vs-Da cases as `unit-incompatible` rather than inventing a conversion based on an arbitrary m/z.
5. Preserve candidate-unavailable and historical-unavailable explicitly.
6. Produce both a per-file report and a 7-accession aggregate report.

This is intentionally a validation/reporting step, not a model-fitting or parameter-tuning step.

## Why this matters

The historical SDRF tolerance is a search-parameter annotation. The v19 candidate represents a RAW-derived measurement-precision-based starting value. A candidate that is numerically close to historical settings is encouraging, but numerical equality is not the target because the quantities are not identical.

The comparison should answer:
- Is v19 generally in the same scale as historical settings?
- Is there systematic under- or over-estimation?
- Does behavior differ by resolution regime or accession?
- Do unit choices remain scientifically coherent?
- Are any cases censored or unsupported?

## Next steps after comparison

1. Run the comparison on all 55 holdout files using the validated SIF Python environment.
2. Inspect the 7-row accession summary and all same-unit file-level ratios.
3. Quantify any systematic bias separately by ppm and Da regimes.
4. Do not tune the `6 x sigma` multiplier until the comparison is understood and a separate validation protocol is defined.
5. If a future tuning study is needed, split datasets into development/training and held-out validation so the 55-file report is not reused for parameter selection.

## Repository quality

The most recent size-guard fix narrows `fileSizeBytes` away from `None` before `int()` and has a regression test for missing API size metadata.

On Codon, use the Singularity image's Python for repository benchmark helper scripts. Do not use the login-node `/usr/bin/python` directly.

Local repository quality command remains:
`scripts/dev/quality.sh --fix`

## Required handoff/commit discipline

For repository changes, provide complete modified repository-relative files in a `.tar.gz`, not a patch, together with a tar extraction command.

After an overlay, always provide:
- extraction command;
- `git add` command;
- `git commit` command;
- relevant Codon experiment commands;
- updated handoff document.
