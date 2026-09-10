# Full-accession prideQC benchmark on Codon / portable Slurm

This profile runs one ProteomeXchange accession per Slurm array task. Inside each task,
prideQC selects all `comment[data file]` entries from the annotated SDRF, downloads those
files from PRIDE into node-local scratch, analyzes them in parallel, and persists only
compact QC/provenance outputs. Downloaded spectra disappear with node-local scratch.

The application/container remains cluster-neutral. Codon paths and resource defaults
live only in the deployment/submission layer.

## Persistent layout

```text
/nfs/research/juan/DIA/singj/prideQC/
  containers/
  scripts/
  benchmarks/data/
  persistent_cache/
  results/
  logs/
```

No source checkout, Cargo build, `uv sync`, or package installation is required on the
cluster.

## 1. Deploy benchmark + launchers from the workstation

```bash
cd /home/sing/Documents/github/prideQC

CLUSTER_HOST=sing@codon-slurm-login-01
PERSIST_ROOT=/nfs/research/juan/DIA/singj/prideQC

ssh "$CLUSTER_HOST" \
  "mkdir -p '$PERSIST_ROOT/containers' '$PERSIST_ROOT/scripts' \
             '$PERSIST_ROOT/benchmarks/data' '$PERSIST_ROOT/persistent_cache' \
             '$PERSIST_ROOT/results' '$PERSIST_ROOT/logs'"

rsync -avP \
  benchmarks/data/prideqc_ground_truth.tsv \
  "$CLUSTER_HOST:$PERSIST_ROOT/benchmarks/data/prideqc_ground_truth.tsv"

rsync -avP \
  scripts/slurm/install_prideqc_sif.sh \
  scripts/slurm/submit_prideqc_accessions.sh \
  scripts/slurm/prideqc_accession_array.sbatch \
  "$CLUSTER_HOST:$PERSIST_ROOT/scripts/"
```

The submit helper derives the number of array tasks from the unique `pxd_accession`
values in the deployed TSV. Nothing is hard-coded to seven or eight accessions.

## 2. Pull the immutable SIF from GHCR on a login node

Prefer the immutable full-SHA tag produced by the GitHub Actions SIF workflow.

```bash
ssh sing@codon-slurm-login-01

export PERSIST_ROOT=/nfs/research/juan/DIA/singj/prideQC
export SIF_REF='oras://ghcr.io/pride-archive/prideqc-sif:sha-<FULL_40_CHAR_GIT_SHA>'

bash "$PERSIST_ROOT/scripts/install_prideqc_sif.sh"
```

The helper smoke-tests prideQC, the exact pyOpenMS environment, native Thermo/Bruker
reader symbols, and .NET runtime availability. It creates local sidecars:

```text
<image>.sif
<image>.sif.sha256
<image>.sif.meta.txt
```

If GHCR requires authentication, log in before the pull. For Apptainer:

```bash
printf '%s' "$GHCR_TOKEN" | apptainer registry login \
  --username '<github-user>' --password-stdin oras://ghcr.io
```

For a public package, no token should be necessary.

## 3. Submit the accession array

Set `SIF` to the path printed by the installation helper.

```bash
export PERSIST_ROOT=/nfs/research/juan/DIA/singj/prideQC
export SIF="$PERSIST_ROOT/containers/prideqc_sha-<FULL_40_CHAR_GIT_SHA>.sif"

export PARTITION=research
export CPUS=17
export PRIDEQC_WORKERS=16
export MEMORY=32G
export TIME_LIMIT=02:00:00
export MAX_PARALLEL=4
export SCRATCH_ROOT=/tmp
export PRIDE_PROTOCOL=ftp
export RUN_NAME=ground-truth-full

bash "$PERSIST_ROOT/scripts/submit_prideqc_accessions.sh"
```

`MAX_PARALLEL=4` means four accessions may run at once, while each accession can analyze
up to 16 independent files concurrently inside its allocation. Tune both values based
on available nodes, PRIDE transfer behavior, file size, and per-file memory use.

For an initial cluster smoke test, use one accession/task at a time by deploying a
small GT TSV or set `MAX_PARALLEL=1`.

## 4. SDRF resolution

Each array task first queries the current `bigbio/sdrf-annotated-datasets` GitHub
directory for its accession and downloads every `*.sdrf.tsv` found there. This handles
both the usual:

```text
datasets/PXD008934/PXD008934.sdrf.tsv
```

and non-canonical/split designs such as:

```text
datasets/PXD018830/PXD018830-DIA.sdrf.tsv
```

If GitHub directory discovery is temporarily unavailable, the launcher falls back to
the `source_sdrf_url` values already stored in the ground-truth TSV.

If an accession has multiple SDRFs, each SDRF is analyzed separately under its own
result/download subdirectory. This is deliberate: it preserves the correct SDRF-to-file
relationship rather than merging experimental designs.

## 5. Scratch/data flow

For every accession:

```text
GHCR SIF on NFS
      -> /tmp/prideqc_<job>_<task>_<PXD>/image.sif

GitHub annotated SDRF
      -> /tmp/.../sdrf/*.sdrf.tsv

PRIDE data files named by SDRF
      -> /tmp/.../downloads/<sdrf>/
      -> prideQC workers
      -> /tmp/.../results/<sdrf>/

compact results/provenance only
      -> /nfs/.../results/<RUN_NAME>/<PXD>/

/tmp/... removed on exit
```

Downloaded RAW/mzML/vendor data are not rsynced to persistent NFS.

## 6. Monitoring

```bash
squeue -u "$USER"
squeue -j <jobid>
sacct -j <jobid> --format=JobID,State,Elapsed,AllocCPUS,MaxRSS,ExitCode
```

Logs are written to:

```text
$PERSIST_ROOT/logs/prideqc-pxd_<jobid>_<array-index>.out
$PERSIST_ROOT/logs/prideqc-pxd_<jobid>_<array-index>.err
```

Follow one task:

```bash
tail -f "$PERSIST_ROOT/logs/prideqc-pxd_<jobid>_1.out"
```

## 7. Inspect results

```bash
find "$PERSIST_ROOT/results/$RUN_NAME" -maxdepth 3 -type f \
  \( -name manifest.json -o -name metrics.tsv -o -name annotations.tsv \
     -o -name refined.sdrf.tsv -o -name hpc-run-info.txt \) -print
```

Per accession you will also retain the exact SDRF(s), PRIDE download manifest(s), GT
rows, SIF checksum, embedded container build information, and Slurm analysis log.

A non-zero array task can still have scientifically useful partial outputs. The launcher
uses `--continue-on-error`, persists compact results first, and then returns non-zero if
one or more files failed. Infrastructure failures before result sync are stored under:

```text
$PERSIST_ROOT/results/$RUN_NAME/_infrastructure_failed/
```

## 8. Important operational notes

- Compute nodes need outbound access to GitHub and the selected PRIDE transfer protocol.
  If FTP is blocked, try another protocol supported by your prideQC/pridepy build, e.g.
  `PRIDE_PROTOCOL=s3` where available.
- Node-local `/tmp` still needs enough free space for the SIF plus the files belonging to
  one accession. This design removes persistent RAW storage pressure; it does not make
  temporary storage requirements disappear.
- `PRIDEQC_WORKERS=16` with `CPUS=17` is the current Codon profile, not an application
  requirement. Reduce workers if `sacct` shows memory pressure.
- BLAS/OpenMP/Rayon thread counts are forced to one per prideQC process to prevent
  oversubscription.
- Existing SDRF assertions are preserved by default. `PRIDEQC_INCLUDE_INFERRED=0` means
  inferred acquisition suggestions are not used to fill SDRF cells unless explicitly
  enabled.

## 9. Commit the infrastructure files

From the workstation repository:

```bash
git add \
  scripts/slurm/install_prideqc_sif.sh \
  scripts/slurm/submit_prideqc_accessions.sh \
  scripts/slurm/prideqc_accession_array.sbatch \
  docs/codon-pride-accessions.md

git diff --cached --check
git diff --cached

git commit -m "Add full-accession Slurm benchmark workflow"
git push origin main
```
