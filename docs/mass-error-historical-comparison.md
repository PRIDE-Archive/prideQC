# Historical SDRF vs v19 fragment-tolerance comparison

The v19 fragment-tolerance candidate is derived from RAW measurement precision. The historical SDRF fragment tolerance is a recorded search parameter. They are related but are not assumed to be the same quantity.

Use `scripts/benchmarks/compare_historical_tolerances.py` to compare them without tuning v19.

## Comparison policy

- Parse the historical fragment tolerance from each result's original `input.sdrf.tsv`.
- Fall back to `benchmarks/data/prideqc_ground_truth.tsv` when the per-run SDRF cannot supply the value.
- Compare candidate and historical values only when both are available and use the same unit (`ppm` or `Da`).
- For same-unit values report candidate/historical ratio, absolute difference, and percent difference.
- Do not convert between ppm and Da at an arbitrary m/z.
- Preserve unavailable and unit-incompatible cases explicitly.

## Codon usage

Run the script inside the validated prideQC Singularity image so the helper uses the project's Python 3.12 environment:

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

The comparison is a validation report only. Do not feed its results back into v19 thresholds or multipliers without a separate, explicitly designed training/validation procedure.
