# Historical SDRF vs v19 fragment-tolerance comparison (v21)

v21 freezes the v19 mass-error estimator and turns the historical-tolerance comparison into a stable validation report. It is not a new estimator and it does not tune the `6 x sigma` multiplier.

Use `scripts/benchmarks/compare_historical_tolerances.py` to compare the RAW-derived candidate with the historical SDRF annotation.

## Comparison policy

- Prefer the fragment-tolerance value recorded in the result's original `input.sdrf.tsv`.
- Fall back to the curated ground-truth TSV only when a per-file SDRF value cannot be recovered.
- Label the historical provenance explicitly (`run-sdrf`, `ground-truth`, or unavailable).
- Calculate the primary ratio, absolute difference, and percent difference only when both values are present and use the same unit (`ppm` or `Da`).
- Never use an arbitrary m/z to convert the primary comparison between Da and ppm.
- For unit-incompatible cases, optionally report a **mass-context equivalence** using the observed MS2 m/z range stored in `metrics.tsv`. The conversion is contextual evidence only; it does not change the primary comparison category.
- Preserve candidate-unavailable, historical-unavailable, and unit-incompatible cases explicitly.

The mass-context equivalence uses the midpoint of `ObservedMzRange_MS2` and reports the corresponding candidate/historical ratio in one common unit. This is intentionally labelled as contextual because a tolerance expressed in Da maps to different ppm values across the observed mass range.

## v21 stabilization goal

The 55-file v19 holdout is the validation set. v21 should be treated as the reporting/finalization layer unless the holdout exposes a concrete scientific defect. Do not retune v19 from this same holdout.

The historical comparison currently shows 24 same-unit comparisons, 10 unit-incompatible comparisons, 11 historical-unavailable comparisons, and 10 candidate-unavailable DIA comparisons across the 55-file holdout.

## Codon usage

Run the utility inside the validated prideQC Singularity image so the helper uses the project's Python environment:

```bash
singularity exec \
  --bind /nfs/research/juan/DIA/singj/prideQC:/prideqc \
  "$SIF" \
  python /prideqc/scripts/benchmarks/compare_historical_tolerances.py \
  --gt /prideqc/benchmarks/data/prideqc_ground_truth.tsv \
  --run-root /prideqc/results/mass-error-v19-holdout-55-nochecksum \
  --files-output /prideqc/benchmarks/manifests/mass-error-v19-holdout-55-historical-comparison-files.tsv \
  --accessions-output /prideqc/benchmarks/manifests/mass-error-v19-holdout-55-historical-comparison-accessions.tsv
```

Inspect the accession report first, then the same-unit file-level cases. Treat mass-context ratios as supporting evidence rather than as replacements for same-unit comparisons.

The comparison is a validation report only. Do not feed its results back into v19 thresholds or multipliers without a separate development/validation protocol.


## v21 stabilization fix

The primary comparison remains unchanged: same-unit values are compared directly, while Da/ppm mismatches remain `unit-incompatible`.

For unit-incompatible pairs, v21 now reports mass-context equivalents at the observed MS2 m/z minimum, observed-range midpoint, and observed MS2 m/z maximum. The midpoint is explicitly a range midpoint; it is not a statistical m/z median because the compact metrics contract does not retain the full m/z distribution. These contextual equivalents are diagnostic only and do not change the primary comparison category or the estimator.
