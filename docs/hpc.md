# Reproducible Docker -> Apptainer/Singularity -> Slurm deployment

This deployment keeps prideQC cluster-agnostic. Docker is the canonical software
build; Apptainer/Singularity is a conversion of that image, not a second build.
Codon is one launch profile only.

## Why a Slurm array is the right shape for PRIDE QC

For the ground-truth benchmark, use one RAW file per Slurm array task. Each task
fetches one selected PRIDE file into node-local scratch, runs prideQC, copies only
QC outputs to persistent storage, and deletes scratch on exit. This prevents a
large collection of vendor RAW files from accumulating on NFS. It does **not**
eliminate storage requirements: each node still needs enough `/tmp` space for its
selected RAW file, temporary data, and result files. `INPUT_MODE=pride` also
requires outbound access from compute nodes to PRIDE (and to GitHub if the SDRF
has not been pre-staged). Use `INPUT_MODE=staged` on clusters without outbound
compute-node networking.

The current prideQC repository is Python-only. `containers/Dockerfile` contains a
real Rust builder stage so that a future `Cargo.toml`/`Cargo.lock` is built with
`cargo build --release --locked --bins`; when no Cargo project exists the image
records that no Rust executable was produced.

## Canonical Docker build

Python is fixed to CPython 3.12 because the required OpenMS development wheel is
CPython-3.12/linux-amd64 specific. The Docker build uses the committed `uv.lock`
for the rest of the environment but deliberately omits the locked stable
`pyopenms` package. `containers/pyopenms.requirements.txt` is a second, tiny
hash-locked input containing exactly this required wheel:

```
https://pypi.openms.de/packages/pyopenms-3.6.0.dev20260910-cp312-cp312-manylinux_2_34_x86_64.whl#sha256=9c7cbb35f3a9557c1988b99a8ec51ef3ac2bc961fd818844e4cfc645ad73a7d6
```

The build fails unless the resulting environment exposes both reader classes
that current prideQC actually calls: `ThermoRawFile` and `BrukerTimsFile`. The
OpenMS Thermo reader uses the .NET bridge on Linux, so the runtime image also
contains the pinned .NET 8.0.29 runtime and its Debian 12 runtime libraries. No
.NET SDK is retained, and no dependency resolution occurs when an HPC job runs.

```bash
export PRIDEQC_IMAGE_REPOSITORY=prideqc
export PRIDEQC_IMAGE_TAG=0.2.0-pyopenms-3.6dev20260910
export PRIDEQC_DOCKER_TARGET=runtime
export PRIDEQC_DOCKER_PLATFORM=linux/amd64
export PRIDEQC_REQUIRE_CLEAN=1
export PRIDEQC_RUN_QUALITY=1
# Optional override; default is pinned to .NET 8.0.29 bookworm-slim.
# export PRIDEQC_DOTNET_IMAGE=mcr.microsoft.com/dotnet/runtime:8.0.29-bookworm-slim

./containers/build_docker.sh
```

`PRIDEQC_RUN_QUALITY=1` builds the `quality` stage first, which runs compileall,
the unittest suite, Ruff checks/format check, and mypy. Set it to `0` only when
you intentionally want to skip that gate.

Inspect provenance:

```bash
docker run --rm prideqc:0.2.0-pyopenms-3.6dev20260910 \
  cat /opt/prideqc/build-info.txt
```

Local Docker example:

```bash
mkdir -p "$PWD/local-results" "$PWD/local-data"
docker run --rm --user "$(id -u):$(id -g)" \
  -v "$PWD/local-data:/data" \
  -v "$PWD/local-results:/results" \
  prideqc:0.2.0-pyopenms-3.6dev20260910 \
  prideqc analyze --workers 3 --overwrite \
    --sdrf /data/study.sdrf.tsv -o /results/run /data/*.raw
```

## Convert the already-built image to SIF

```bash
IMAGE="prideqc:0.2.0-pyopenms-3.6dev20260910"
SIF="$PWD/containers/prideqc_0.2.0-pyopenms-3.6dev20260910.sif"

./containers/build_singularity.sh "$IMAGE" "$SIF"
```

This generates:

```
<image>.sif
<image>.sif.sha256
<image>.sif.meta.txt
```

Verify it with:

```bash
cd containers
sha256sum -c prideqc_0.2.0-pyopenms-3.6dev20260910.sif.sha256
```

The conversion script prefers Apptainer, falls back to Singularity, verifies the
SIF checksum, checks prideQC/Python/pyOpenMS, and smoke-tests any Rust executables
listed by the image.

## Codon deployment profile

Codon paths are deployment configuration, not application behavior.

```bash
CLUSTER_HOST=sing@codon-slurm-login-01
PERSIST_ROOT=/nfs/research/juan/DIA/singj/prideQC
IMAGE_TAG=0.2.0-pyopenms-3.6dev20260910
SIF="$PWD/containers/prideqc_${IMAGE_TAG}.sif"

ssh "$CLUSTER_HOST" \
  "mkdir -p '$PERSIST_ROOT/containers' '$PERSIST_ROOT/scripts' \
              '$PERSIST_ROOT/datasets' '$PERSIST_ROOT/datasets/sdrf' \
              '$PERSIST_ROOT/models' '$PERSIST_ROOT/persistent_cache' \
              '$PERSIST_ROOT/results' '$PERSIST_ROOT/logs'"

rsync -avP \
  "$SIF" \
  "$SIF.sha256" \
  "$SIF.meta.txt" \
  "$CLUSTER_HOST:$PERSIST_ROOT/containers/"

# Intentionally flatten scripts/slurm -> remote scripts/.
rsync -avP \
  scripts/slurm/prideqc_ground_truth_array.sbatch \
  "$CLUSTER_HOST:$PERSIST_ROOT/scripts/prideqc_ground_truth_array.sbatch"

rsync -avP \
  benchmarks/prideqc_ground_truth.tsv \
  "$CLUSTER_HOST:$PERSIST_ROOT/datasets/prideqc_ground_truth.tsv"
```

Verify the persistent artifact after transfer:

```bash
ssh "$CLUSTER_HOST" \
  "cd '$PERSIST_ROOT/containers' && \
   sha256sum -c 'prideqc_${IMAGE_TAG}.sif.sha256'"
```

## Submit the ground-truth array

The benchmark TSV has a header, so array tasks are numbered from 1 through the
number of data rows. The launcher resolves column names from the header rather
than relying on fixed column positions.

Codon example, with five RAW files allowed to run concurrently:

```bash
CLUSTER_HOST=sing@codon-slurm-login-01
PERSIST_ROOT=/nfs/research/juan/DIA/singj/prideQC
IMAGE_TAG=0.2.0-pyopenms-3.6dev20260910

ssh "$CLUSTER_HOST" bash -s <<EOF2
set -euo pipefail
cd "$PERSIST_ROOT"

GT_FILE="$PERSIST_ROOT/datasets/prideqc_ground_truth.tsv"
N_TASKS=\$(( \$(wc -l < "\$GT_FILE") - 1 ))

PARTITION=research
CPUS=1
MEMORY=32G
TIME_LIMIT=02:00:00
MAX_PARALLEL=5
SCRATCH_ROOT=/tmp
SIF="$PERSIST_ROOT/containers/prideqc_${IMAGE_TAG}.sif"
DATA_ROOT="$PERSIST_ROOT/datasets"
MODEL_ROOT="$PERSIST_ROOT/models"
CACHE_ROOT="$PERSIST_ROOT/persistent_cache"
RESULT_ROOT="$PERSIST_ROOT/results"

sbatch \
  --partition="\$PARTITION" \
  --cpus-per-task="\$CPUS" \
  --mem="\$MEMORY" \
  --time="\$TIME_LIMIT" \
  --array="1-\${N_TASKS}%\${MAX_PARALLEL}" \
  --export=ALL,PERSIST_ROOT="$PERSIST_ROOT",SIF="\$SIF",GT_FILE="\$GT_FILE",SCRATCH_ROOT="\$SCRATCH_ROOT",DATA_ROOT="\$DATA_ROOT",MODEL_ROOT="\$MODEL_ROOT",CACHE_ROOT="\$CACHE_ROOT",RESULT_ROOT="\$RESULT_ROOT",INPUT_MODE=pride,PRIDEQC_WORKERS=1 \
  "$PERSIST_ROOT/scripts/prideqc_ground_truth_array.sbatch"
EOF2
```

For this one-file-per-task benchmark, `PRIDEQC_WORKERS=1` avoids nested
parallelism and thread oversubscription. Scale with the array concurrency
(`MAX_PARALLEL`) instead. The launcher sets BLAS/OpenMP/NumExpr/Rayon threads to
one. If a later launcher analyzes multiple independent files in one allocation,
set `PRIDEQC_WORKERS` to an appropriate value bounded by `SLURM_CPUS_PER_TASK`.

To use already-staged RAW files instead of PRIDE downloads:

```bash
# Expected input: $DATA_ROOT/<PXD>/<exact RAW filename>
INPUT_MODE=staged
```

The launcher always stages the selected input from NFS to `/tmp` before analysis.

## Optional SDRF pre-staging

If compute nodes cannot access GitHub, place annotated SDRFs under
`$DATA_ROOT/sdrf/`. The launcher first looks there using the basename from each
GT row's `source_sdrf_url`; only if it is absent does it download the SDRF into
scratch.

For example, on a login node with network access:

```bash
mkdir -p "$PERSIST_ROOT/datasets/sdrf"
# Populate this directory with the seven SDRFs referenced in
# benchmarks/prideqc_ground_truth.tsv, preserving their basenames.
```

## Monitoring

```bash
squeue -u "$USER"
squeue -j <jobid>
sacct -j <jobid> --format=JobID,State,Elapsed,AllocCPUS,MaxRSS,ExitCode

tail -f "$PERSIST_ROOT/logs/prideqc-gt_<jobid>_<task>.out"
```

## Results and failure diagnostics

Successful jobs land under:

```
$RESULT_ROOT/<GT_ID>_<PXD>_<RAW filename>/
```

The RAW download remains only in node-local scratch and is deleted. Results
include normal prideQC outputs plus `slurm-run.log` and `hpc-run-info.txt`.

```bash
find "$PERSIST_ROOT/results" -name manifest.json -print
find "$PERSIST_ROOT/results" -name '*.summary.json' -print
find "$PERSIST_ROOT/results" -name refined.sdrf.tsv -print
```

On failure, compact diagnostics and any partial prideQC outputs are copied to
`$RESULT_ROOT/_failed/`; the large downloaded RAW is deliberately not copied to
NFS merely because a job failed.

## Cache behavior

By default, job caches live under node-local scratch (`XDG_CACHE_HOME=/work/cache`)
and vanish with the job. For a future workload where a persistent cache is
valuable:

```
STAGE_CACHE=1
SYNC_CACHE_BACK=1
CACHE_ROOT=/persistent/path
```

The launcher only writes the local cache back after successful scientific work.

## Portability knobs

The image and launcher do not contain Codon paths. Keep these outside the
application and override them per cluster/workload:

```
CLUSTER_HOST
PERSIST_ROOT
PARTITION
CPUS
MEMORY
TIME_LIMIT
SCRATCH_ROOT
SIF
DATA_ROOT
MODEL_ROOT
CACHE_ROOT
RESULT_ROOT
```

The same image can be run under Docker, Apptainer/Singularity, local shell use,
Slurm, or another scheduler launcher.

## Validation before commit

```bash
bash -n containers/build_docker.sh
bash -n containers/build_singularity.sh
bash -n scripts/slurm/prideqc_ground_truth_array.sbatch

uv sync --locked --group dev
uv run --locked python -m compileall -q src tests
uv run --locked python -m unittest discover -s tests -v
uv run --locked ruff check .
uv run --locked ruff format --check .
uv run --locked mypy
uv build
```

Then perform the Docker quality/runtime and SIF smoke tests described above.

## Intended Git commit

Do not use `git add .`. Add only the intended infrastructure and benchmark files:

```bash
git add \
  .dockerignore \
  containers/Dockerfile \
  containers/.gitignore \
  containers/build_docker.sh \
  containers/build_singularity.sh \
  containers/pyopenms.requirements.txt \
  scripts/slurm/prideqc_ground_truth_array.sbatch \
  docs/hpc.md \
  benchmarks/prideqc_ground_truth.tsv

git status --short
git diff --cached --check
git commit -m "Add reproducible Docker Apptainer Slurm deployment"
```

## GitHub Actions: publish Docker and SIF images to GHCR

Two workflows under `.github/workflows/` automate the same canonical build path:

- `container-docker.yml` runs the Docker `quality` stage, builds the `runtime`
  image for `linux/amd64`, and publishes it to `ghcr.io/<owner>/<repo>` on pushes
  to `main`, version tags (`v*`), and manual runs. Pull requests build but do not
  publish. Every published build includes the immutable full-commit tag
  `sha-<40-character-git-sha>`.
- `container-sif.yml` starts only after a successful non-PR Docker workflow (or
  manually), pulls the immutable Docker tag, resolves its registry digest,
  converts that already-built image to SIF, verifies the SIF, uploads the three
  SIF files as a GitHub Actions artifact, and pushes the SIF to the separate
  `ghcr.io/<owner>/<repo>-sif` package via Apptainer's `oras://` transport.

The Docker and SIF packages are deliberately separate. This prevents Docker/OCI
image manifests and native SIF/ORAS artifacts from competing for the same tag
namespace while keeping both artifacts in GitHub Container Registry.

Typical published references for commit `$SHA` are:

```bash
# Docker / OCI image
docker pull ghcr.io/pride-archive/prideqc:sha-$SHA

# Native Apptainer SIF stored as an OCI/ORAS artifact
apptainer pull prideqc.sif \
  oras://ghcr.io/pride-archive/prideqc-sif:sha-$SHA
```

A successful build from `main` also publishes the SIF alias `main`. A stable
`vX.Y.Z` tag publishes the version alias and `latest`; pre-release version tags
publish their version alias but do not replace `latest`.

The SIF metadata records both the local Docker image ID and, after the image has
been pulled from GHCR, its immutable `docker_repo_digest`. The SIF workflow
checks that digest before publishing the artifact. This provides an explicit
chain:

```
Git commit -> Docker image digest -> SIF SHA-256
```

The workflows use the repository-provided `GITHUB_TOKEN`; no personal registry
secret is required for a repository whose Actions workflow token is allowed to
write packages. If an organization policy overrides repository workflow
permissions, enable package write access for Actions or adjust that policy.

GitHub Actions uses the same pinned pyOpenMS requirement file and Dockerfile as
local builds. Python dependencies are therefore resolved only while the Docker
image is built; neither SIF conversion nor Slurm execution runs `uv sync` or
`pip install`.

### Manual runs

From the GitHub Actions UI, `Container - Docker` can be run manually for the
selected ref. After it succeeds, `Container - Apptainer SIF` normally follows
automatically. The SIF workflow can also be dispatched manually. Its optional `source_sha` is
a full Git commit SHA; when omitted it uses the selected ref's commit. The
workflow always converts the corresponding immutable Docker tag
`sha-<source_sha>`. The optional `publish_alias` adds one extra SIF tag.

### HPC pull instead of rsync

If compute/login nodes are permitted to access GHCR, the SIF no longer has to be
built locally or copied from a workstation. On a login node:

```bash
mkdir -p "$PERSIST_ROOT/containers"
cd "$PERSIST_ROOT/containers"
apptainer pull prideqc_${SHA}.sif \
  "oras://ghcr.io/pride-archive/prideqc-sif:sha-${SHA}"
sha256sum prideqc_${SHA}.sif
```

For a locked scientific deployment, use the immutable `sha-...` SIF reference,
record its local SHA-256, and point the Slurm launcher at that file. `main` and
`latest` are convenience aliases and should not be used as the definitive
scientific provenance identifier.
