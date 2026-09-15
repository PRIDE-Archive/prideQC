# prideQC Project Handoff — v21 — 2026-09-15

## Project goal

Build a PRIDE-scale quality-control pipeline that can process mass-spectrometry data and collect QC metrics/technical metadata without modifying the underlying scientific data.

The benchmark strategy uses SDRF-annotated PRIDE datasets as ground truth and runs production-style per-file processing on Codon using Singularity + Slurm.

## Current baseline

### Runtime / deployment

- Codon uses **Singularity**, not Docker commands.
- Production execution is one PRIDE data file per Slurm array task.
- `prideqc --workers 1` is used for predictable node memory.
- Large vendor RAW files are staged under node-local `/tmp` and only compact result artifacts are persisted to NFS.
- The validated image is pulled from GHCR through the Codon Singularity/ORAS helper and receives a local `.sif.sha256` automatically.
- Current validated v19 runtime image:
  `prideqc_sha-85090462ce5ed9806764aaa5c71718bb6894d7ba.sif`
- pyOpenMS:
  `3.6.0.dev20260910`
- Native Thermo RAW and Bruker tims readers are available in the validated image.

### Ground truth / manifests

Ground-truth controls currently cover 10 files across 7 accessions:

- PXD008934: GT001
- PXD017618: GT002, GT003
- PXD001819: GT004
- PXD000759: GT005
- PXD000612: GT006
- PXD003772: GT007, GT008, GT009
- PXD018830: GT010

The reconciled v6 physical manifest contains 663 physical analysis tasks plus 42 alias records. PXD000612 contains the known 42 `pY_...` alternate filename references; they are treated as aliases rather than duplicate physical analyses.

## Acquisition inference status

The acquisition-method inference is v11 and is considered validated on the four-control set:

- narrow DDA controls infer Data-dependent acquisition;
- the wide-window fixed-cycle PXD018830 DIA control infers Data-independent acquisition;
- ambiguous/small fixed-target narrow runs abstain.

The four-control validation report passed all expected state checks.

## v19 mass-error baseline

v19 is the current frozen mass-error/tolerance baseline.

It provides measurement-precision evidence from repeated observations and produces a candidate fragment search tolerance only when support is sufficient.

Key policy decisions:

- high-resolution data use ppm tolerance candidates;
- low-resolution data use adaptive 0.5 or 1.0 Da measurement windows;
- profile spectra are never centroided in-place or modified; profile peak centers are ephemeral analytical estimates;
- ambiguous acquisition modes, especially DIA cases where a search tolerance cannot be justified from the available evidence, abstain;
- historical SDRF search parameters are not reconstructed from RAW data and are not written into SDRF automatically.

## 55-file holdout

Run root:
`results/mass-error-v19-holdout-55-nochecksum`

Completed scientific analyses:
- 55/55 usable analyses
- 5 `success / ok`
- 50 `warning / refined_sdrf_validation`
- 0 scientific-analysis failures after retrying the one transient PXD018830 download failure with the alternative `globus` protocol

Accession-level inference:
- PXD000612: 10/10 inferred, ppm, high-resolution, 0.2 Da window; median candidate ~10.21 ppm.
- PXD000759: 10/10 inferred, Da, low-resolution; adaptive 0.5/1.0 Da windows; median candidate ~0.666 Da.
- PXD001819: 10/10 inferred, Da, low-resolution, 0.5 Da window; median candidate ~0.396 Da.
- PXD003772: 1/1 inferred, ppm, high-resolution, 0.2 Da; ~10.47 ppm. This accession has weak holdout coverage.
- PXD008934: 10/10 inferred, ppm, high-resolution, 0.2 Da; median candidate ~14.30 ppm.
- PXD017618: 4/4 inferred, ppm, high-resolution, 0.2 Da; median candidate ~7.83 ppm; two high-confidence and two moderate-confidence cases.
- PXD018830: 10/10 unavailable; DIA cases correctly abstain.

The strongest holdout conclusion is categorical consistency: 45/45 non-DIA files infer a candidate, 10/10 DIA files abstain, and the expected state/unit regime is maintained across the benchmark.

## Historical SDRF comparison result

v20 introduced the reporting utility:
`scripts/benchmarks/compare_historical_tolerances.py`

The 55-file comparison produced:

- 24 same-unit comparisons
- 10 unit-incompatible comparisons
- 11 historical-unavailable comparisons
- 10 candidate-unavailable comparisons

Direct same-unit accession-level observations:

- PXD000759: historical 0.6 Da; median v19 candidate ~0.666 Da; ratio ~1.11x.
- PXD008934: historical 20 ppm; median v19 candidate ~14.30 ppm; ratio ~0.715x.
- PXD017618: historical median 7.5 ppm across the four comparable files; median v19 candidate ~7.83 ppm; ratio ~1.025x.

The PXD000612 cases are intentionally unit-incompatible because historical SDRF says 0.05 Da while v19 reports ppm. No arbitrary m/z conversion is used in the primary comparison.

PXD018830 has no historical fragment search tolerance and v19 correctly abstains.

## v21 purpose: stabilization, not endless iteration

v21 is intended to be the **last comparison/reporting iteration unless a concrete validation failure is found**.

It extends the v20 comparison utility in three bounded ways:

1. Historical provenance is explicit: per-run SDRF first, curated GT fallback, otherwise unavailable.
2. Unit-incompatible cases may receive a secondary mass-context equivalence using the observed MS2 m/z range midpoint. This is supporting evidence only and never replaces the same-unit comparison.
3. The accession summary reports historical-unavailable counts and mass-context ratio summaries in addition to the v20 statistics.

There is **no change to the v19 mass-error estimator or its `6 x sigma` multiplier** in v21.

## Scientific interpretation / decision rule

Historical SDRF tolerance is a recorded search parameter. The v19 candidate is derived from measured repeat-observation precision. They are related but are not assumed to be identical.

The existing holdout evidence is sufficient to say that the v19 candidate is generally in the same scale as historical settings for directly comparable files, while not being a trivial copy of the historical values.

The next scientific decision should therefore be a single explicit one:

> Is the measurement-derived candidate, together with its uncertainty/inlier diagnostics and resolution-aware matching window, sufficiently defensible as a general search-tolerance recommendation?

Do not keep adding heuristic multipliers to chase historical SDRF values. If further tuning is ever required, use a separate development set and keep the 55-file holdout untouched as validation.

## Repository quality

Before pushing repository changes, run:

`scripts/dev/quality.sh --fix`

and the focused benchmark test:

`uv run --no-sync python -m unittest tests.test_historical_tolerance_comparison -v`

The local quality gate should leave Ruff and mypy clean.

## v21 files

- `scripts/benchmarks/compare_historical_tolerances.py`
- `tests/test_historical_tolerance_comparison.py`
- `docs/mass-error-historical-comparison.md`
- `prideqc-project-handoff-20260915-v21.md`

## Recommended final workflow

1. Extract the v21 overlay into the current repository.
2. Run the quality gate and focused tests.
3. Commit and push.
4. Build/publish the Docker image and Singularity image through the existing CI pipeline.
5. Run the v21 comparison helper against the frozen 55-file holdout using the validated SIF Python environment.
6. Review the resulting accession summary and file-level comparison once.
7. Treat v19/v21 as the frozen baseline for the next broader benchmark rather than another cycle of heuristic tuning.

## Reproducibility artifacts

Existing v19 holdout reports:

- `benchmarks/manifests/mass-error-v19-holdout-55-files.tsv`
- `benchmarks/manifests/mass-error-v19-holdout-55-accessions.tsv`
- `benchmarks/manifests/mass-error-v19-holdout-55-historical-comparison-files.tsv`
- `benchmarks/manifests/mass-error-v19-holdout-55-historical-comparison-accessions.tsv`

The v20 reports should remain unchanged as the record of the same-unit comparison; v21 adds contextual mass-aware reporting without altering those baseline results.
