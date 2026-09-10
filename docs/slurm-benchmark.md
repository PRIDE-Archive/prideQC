# prideQC portable Slurm benchmark launcher

These launchers are intended to be portable Slurm infrastructure. Cluster-specific
paths and resource choices are supplied via environment variables at submission time.
Do not hard-code Codon paths into prideQC itself.

The file-level topology is the recommended benchmark mode: one PRIDE archive file per
Slurm array element and `prideqc --workers 1`. This avoids multiplying OpenMS memory
usage inside a single allocation.

The manifest preserves two filenames when present:

- `sdrf_data_file`: `comment[data file]`, used to associate QC evidence back to SDRF rows.
- `archive_file`: basename of `comment[file uri]`, used for the actual PRIDE download.

When they differ, the launcher supplies prideQC `--file-map` automatically. This handles
annotated datasets such as PXD001819, whose SDRF references an mzML analysis filename but
whose archive URI points to the corresponding RAW file.

Codon is a deployment profile, not a runtime requirement. Example site values remain
external to the scripts:

```bash
export PERSIST_ROOT=/nfs/research/juan/DIA/singj/prideQC
export PARTITION=research
export SCRATCH_ROOT=/tmp
export CPUS=2
export MEMORY=32G
export TIME_LIMIT=04:00:00
export MAX_PARALLEL=4
```

For a new benchmark run, use a new `RUN_NAME` so previous results are retained.
