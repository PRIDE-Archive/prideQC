# Codon PRIDE full-dataset execution: file-level Slurm array

## Why file-level parallelism

`prideqc --workers N` uses independent Python processes for independent files. Native OpenMS/vendor readers may hold substantial per-file state, so a 32 GiB allocation is not a safe place to run 16 file processes concurrently. On Codon, use one PRIDE data file per Slurm task with `--workers 1`; let the Slurm array provide the parallelism.

This also minimizes scratch use: each task holds only the SIF, one downloaded PRIDE data file, temporary data, and QC output.

## Persistent layout

```text
/nfs/research/juan/DIA/singj/prideQC/
  containers/
  scripts/
  benchmarks/
    data/prideqc_ground_truth.tsv
    sdrf/
    manifests/
  results/
  logs/
```

SDRFs and the generated file manifest are tiny and are intentionally persisted. Raw/mzML/vendor files remain under node-local `/tmp` and are deleted after the task.

## Deploy the scripts

From the workstation repository:

```bash
export CLUSTER_HOST=sing@codon-slurm-login-01
export PERSIST_ROOT=/nfs/research/juan/DIA/singj/prideQC

rsync -avP \
  scripts/slurm/prepare_prideqc_file_manifest.py \
  scripts/slurm/prepare_ground_truth_task_subset.py \
  scripts/slurm/submit_prideqc_files.sh \
  scripts/slurm/prideqc_file_array.sbatch \
  scripts/slurm/inspect_prideqc_array.sh \
  "$CLUSTER_HOST:$PERSIST_ROOT/scripts/"

rsync -avP \
  benchmarks/data/prideqc_ground_truth.tsv \
  "$CLUSTER_HOST:$PERSIST_ROOT/benchmarks/data/"
```

## Submit

On Codon, use the immutable SHA-tagged SIF already pulled from GHCR:

```bash
export PERSIST_ROOT=/nfs/research/juan/DIA/singj/prideQC
export SIF="$PERSIST_ROOT/containers/<your-immutable-image>.sif"

export PARTITION=research
export CPUS=2
export MEMORY=32G
export TIME_LIMIT=04:00:00
export MAX_PARALLEL=4
export PRIDE_PROTOCOL=ftp
export RUN_NAME=ground-truth-files-v1

bash "$PERSIST_ROOT/scripts/submit_prideqc_files.sh"
```

The submit helper downloads the annotated SDRF(s) for every unique PXD in the GT TSV, parses `comment[data file]`, deduplicates repeated SDRF rows, writes a reproducible manifest, and submits `1-N%MAX_PARALLEL` where N is the total number of unique PXD/SDRF/data-file tasks.

Inspect the generated manifest before or after submission:

```bash
column -ts $'\t' "$PERSIST_ROOT/benchmarks/manifests/ground-truth-files-v1.tsv" | less -S
cat "$PERSIST_ROOT/benchmarks/manifests/ground-truth-files-v1.meta.txt"
```


## Reusable ground-truth task subsets

For targeted benchmark experiments, do not rebuild Slurm task lists with ad hoc
CSV snippets. Use the repository helper so task IDs are resolved against the
active manifest and the generated TSV always uses Unix LF line endings:

```bash
python scripts/slurm/prepare_ground_truth_task_subset.py \
  --gt "$PERSIST_ROOT/benchmarks/data/prideqc_ground_truth.tsv" \
  --manifest "$PERSIST_ROOT/benchmarks/manifests/ground-truth-files-v6.tsv" \
  --output "$PERSIST_ROOT/benchmarks/manifests/ground-truth-v6-tasks.tsv"
```

The helper prints validated `task_count=...` and `task_array=...` values that can be passed directly
to Slurm. Prefer the reported `task_count` over counting newlines in the comma-separated array string, because a shell variable normally has no trailing newline:

```bash
HELPER_OUT=$(python scripts/slurm/prepare_ground_truth_task_subset.py \
  --gt "$PERSIST_ROOT/benchmarks/data/prideqc_ground_truth.tsv" \
  --manifest "$PERSIST_ROOT/benchmarks/manifests/ground-truth-files-v6.tsv" \
  --output "$PERSIST_ROOT/benchmarks/manifests/ground-truth-v6-tasks.tsv")

printf '%s\n' "$HELPER_OUT"
TASK_COUNT=$(printf '%s\n' "$HELPER_OUT" | awk -F= '$1 == "task_count" {print $2}')
TASKS=$(printf '%s\n' "$HELPER_OUT" | awk -F= '$1 == "task_array" {print $2}')

printf 'TASK_COUNT=%s\nTASKS=<%s>\n' "$TASK_COUNT" "$TASKS"
```

The helper writes TSV files with `lineterminator="\\n"`. This avoids hidden CR
characters in the final `task_id` column, which otherwise make a visually valid
Slurm array specification fail with `Invalid job array specification`.

### End-to-end targeted experiment rule

For targeted scientific validation, run the experiment from a fresh shell block
instead of reusing variables or report names from a previous version. A complete
experiment should always define, in order: `PERSIST_ROOT`, immutable Git SHA/SIF,
`GT`, `MANIFEST`, task subset/output path, `RESULT_ROOT`, `RUN_NAME`, estimator
flags, Slurm submission, and the version-specific report. Do not rely on an older
`RUN_NAME`, `SIF`, `GT_TASKS`, or report script remaining correct in the shell.

Before submitting, print the immutable container provenance and the resolved array:

```bash
singularity exec "$SIF" cat /opt/prideqc/build-info.txt |
  grep -E '^(vcs_ref|source_state|source_fingerprint_sha256)='

printf 'RUN_NAME=%s\nTASK_COUNT=%s\nTASKS=<%s>\n' \
  "$RUN_NAME" "$TASK_COUNT" "$TASKS"

test "$TASK_COUNT" -eq 10
printf '%s\n' "$TASKS" | grep -Eq '^[0-9]+(,[0-9]+)*$'
```

Use a new run name for every estimator revision even when the same RAW files are
reused. This makes provenance visible in the result path and prevents a correct new
SIF from being mistaken for an older experiment because stale shell names were reused.

For mass-error estimator v6/v19, repeated-spectrum identity is still established with
the 0.2 Da fragment-overlap window. Low-resolution precision is then measured with an
adaptive ladder of 0.5 Da followed by 1.0 Da. The estimator uses the narrowest window
whose robust three-sigma pairwise core stays below 90% of the window and abstains if
1.0 Da remains censored. Do not change the Slurm wrapper or RAW data to pre-centroid
profile spectra; profile peak centers are derived ephemerally inside prideQC.

## Monitor

```bash
squeue -u "$USER"
sacct -j <jobid> --units=G --format=JobID,State,Elapsed,AllocCPUS,ReqMem,MaxRSS,ExitCode
bash "$PERSIST_ROOT/scripts/inspect_prideqc_array.sh" <jobid> "$PERSIST_ROOT/logs"
```

## Results

Each task persists only compact outputs:

```text
results/<run>/<PXD>/<SDRF>/<task>_<data-file>/
  manifest.json
  *.mzQC
  *.summary.json
  metrics.tsv
  annotations.tsv
  refined.sdrf.tsv
  sdrf-changes.tsv
  download-manifest.json
  input.sdrf.tsv
  hpc-analysis.log
  hpc-run-info.txt
  container-build-info.txt
```

Downloaded spectra are not copied back to NFS.

## Memory tuning

Start with `MEMORY=32G`, `MAX_PARALLEL=4`, and one prideQC worker per file. If an individual task still reaches 32 GiB, that is a genuine single-file memory requirement rather than worker multiplication. Re-submit that file/task with 48 or 64 GiB as appropriate. Do not increase `--workers` for file-level jobs.

The old accession-level launcher can still be useful for small mzML projects, but it should not be used with 16 workers under a 32 GiB memory cap for heterogeneous vendor RAW files.
