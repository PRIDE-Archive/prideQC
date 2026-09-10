#!/usr/bin/env bash
#SBATCH --job-name=prideqc-gt
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

set -euo pipefail

: "${PERSIST_ROOT:?Set PERSIST_ROOT to the persistent project root}"
: "${SIF:?Set SIF to the persistent .sif path}"
GT_FILE="${GT_FILE:-$PERSIST_ROOT/datasets/prideqc_ground_truth.tsv}"
RESULT_ROOT="${RESULT_ROOT:-$PERSIST_ROOT/results}"
DATA_ROOT="${DATA_ROOT:-$PERSIST_ROOT/datasets}"
MODEL_ROOT="${MODEL_ROOT:-$PERSIST_ROOT/models}"
CACHE_ROOT="${CACHE_ROOT:-$PERSIST_ROOT/persistent_cache}"
SCRATCH_ROOT="${SCRATCH_ROOT:-/tmp}"
INPUT_MODE="${INPUT_MODE:-pride}"            # pride | staged
PRIDE_PROTOCOL="${PRIDE_PROTOCOL:-ftp}"
PRIDEQC_WORKERS="${PRIDEQC_WORKERS:-1}"
STAGE_CACHE="${STAGE_CACHE:-0}"
SYNC_CACHE_BACK="${SYNC_CACHE_BACK:-0}"

TASK_ID="${SLURM_ARRAY_TASK_ID:-${TASK_ID:-1}}"
JOB_ID="${SLURM_JOB_ID:-local-$$}"
CPUS="${SLURM_CPUS_PER_TASK:-${CPUS:-1}}"
START_EPOCH="$(date +%s)"

if command -v apptainer >/dev/null 2>&1; then
    CONTAINER_BIN=apptainer
elif command -v singularity >/dev/null 2>&1; then
    CONTAINER_BIN=singularity
else
    echo "Neither apptainer nor singularity is available" >&2
    exit 2
fi
command -v rsync >/dev/null 2>&1 || { echo "rsync is required" >&2; exit 2; }
command -v sha256sum >/dev/null 2>&1 || { echo "sha256sum is required" >&2; exit 2; }

test -r "$GT_FILE" || { echo "Ground-truth TSV not readable: $GT_FILE" >&2; exit 2; }
test -r "$SIF" || { echo "SIF not readable: $SIF" >&2; exit 2; }
test -r "${SIF}.sha256" || { echo "SIF checksum file not readable: ${SIF}.sha256" >&2; exit 2; }

get_field() {
    local name="$1"
    awk -F '\t' -v wanted="$name" -v row="$((TASK_ID + 1))" '
        NR == 1 {
            for (i = 1; i <= NF; i++) if ($i == wanted) col = i
            if (!col) exit 3
            next
        }
        NR == row { print $col; found = 1; exit }
        END { if (!found) exit 4 }
    ' "$GT_FILE"
}

GT_ID="$(get_field gt_id)"
PXD="$(get_field pxd_accession)"
RAW_FILE="$(get_field raw_file)"
SDRF_URL="$(get_field source_sdrf_url)"
ACQUISITION="$(get_field acquisition_method)"
[[ -n "$GT_ID" && -n "$PXD" && -n "$RAW_FILE" && -n "$SDRF_URL" ]] || {
    echo "Incomplete GT row $TASK_ID" >&2; exit 2;
}

if [[ "$ACQUISITION" == DIA ]]; then
    SDRF_TEMPLATE=dia-acquisition
else
    SDRF_TEMPLATE=ms-proteomics
fi

safe_raw="$(printf '%s' "$RAW_FILE" | tr -c 'A-Za-z0-9._-' '_')"
RUN_NAME="${RUN_NAME:-${GT_ID}_${PXD}_${safe_raw}}"
WORK_ROOT="$SCRATCH_ROOT/prideqc_${JOB_ID}_${TASK_ID}"
LOCAL_SIF="$WORK_ROOT/$(basename "$SIF")"
LOCAL_SDRF="$WORK_ROOT/input.sdrf.tsv"
LOCAL_RESULTS="$WORK_ROOT/results"
LOCAL_DOWNLOADS="$WORK_ROOT/downloads"
LOCAL_CACHE="$WORK_ROOT/cache"
LOCAL_TMP="$WORK_ROOT/tmp"
FINAL_RESULTS="$RESULT_ROOT/$RUN_NAME"
FAIL_RESULTS="$RESULT_ROOT/_failed/${RUN_NAME}_${JOB_ID}"
RUN_LOG="$WORK_ROOT/run.log"

mkdir -p "$WORK_ROOT" "$LOCAL_RESULTS" "$LOCAL_DOWNLOADS" "$LOCAL_CACHE" "$LOCAL_TMP" "$RESULT_ROOT"

failure_diagnostics() {
    local status=$?
    local elapsed=$(( $(date +%s) - START_EPOCH ))
    if [[ $status -ne 0 ]]; then
        echo "Job failed with status $status after ${elapsed}s" >&2
        mkdir -p "$FAIL_RESULTS"
        {
            echo "status=$status"
            echo "job_id=$JOB_ID"
            echo "task_id=$TASK_ID"
            echo "node=$(hostname)"
            echo "run_name=$RUN_NAME"
            echo "work_root=$WORK_ROOT"
            echo "elapsed_seconds=$elapsed"
            echo
            echo "[disk]"
            df -h "$WORK_ROOT" || true
            echo
            echo "[scratch listing]"
            find "$WORK_ROOT" -maxdepth 3 -type f -printf '%p\t%s bytes\n' 2>/dev/null | sort || true
        } > "$FAIL_RESULTS/failure-diagnostics.txt"
        [[ -f "$RUN_LOG" ]] && cp -a "$RUN_LOG" "$FAIL_RESULTS/" || true
        [[ -d "$LOCAL_RESULTS" ]] && rsync -a "$LOCAL_RESULTS/" "$FAIL_RESULTS/partial-results/" || true
        [[ -f "$LOCAL_DOWNLOADS/download-manifest.json" ]] && cp -a "$LOCAL_DOWNLOADS/download-manifest.json" "$FAIL_RESULTS/" || true
        [[ -f "$LOCAL_SDRF" ]] && cp -a "$LOCAL_SDRF" "$FAIL_RESULTS/input.sdrf.tsv" || true
        echo "Failure diagnostics: $FAIL_RESULTS" >&2
    fi
    rm -rf "$WORK_ROOT"
    exit "$status"
}
trap failure_diagnostics EXIT

printf '%s\n' \
    "job_id=$JOB_ID" \
    "array_task=$TASK_ID" \
    "node=$(hostname)" \
    "slurm_partition=${SLURM_JOB_PARTITION:-unknown}" \
    "cpus=$CPUS" \
    "slurm_mem_per_node=${SLURM_MEM_PER_NODE:-unknown}" \
    "gt_id=$GT_ID" \
    "pxd=$PXD" \
    "raw_file=$RAW_FILE" \
    "acquisition=$ACQUISITION" \
    "sdrf_template=$SDRF_TEMPLATE" \
    "input_mode=$INPUT_MODE" \
    "persist_root=$PERSIST_ROOT" \
    "data_root=$DATA_ROOT" \
    "model_root=$MODEL_ROOT" \
    "cache_root=$CACHE_ROOT" \
    "result_root=$RESULT_ROOT" \
    "scratch_root=$SCRATCH_ROOT" \
    "run_name=$RUN_NAME"

echo "==> Verify persistent SIF"
(
    cd "$(dirname "$SIF")"
    sha256sum -c "$(basename "${SIF}.sha256")"
)
EXPECTED_SIF_SHA="$(awk 'NR==1 {print $1}' "${SIF}.sha256")"

echo "==> Stage SIF to node-local scratch"
cp -a "$SIF" "$LOCAL_SIF"
LOCAL_SIF_SHA="$(sha256sum "$LOCAL_SIF" | awk '{print $1}')"
[[ "$LOCAL_SIF_SHA" == "$EXPECTED_SIF_SHA" ]] || { echo "Staged SIF checksum mismatch" >&2; exit 2; }

if [[ -r "${SIF}.meta.txt" ]]; then
    echo "==> SIF metadata"
    cat "${SIF}.meta.txt"
fi

echo "==> Embedded software provenance"
"$CONTAINER_BIN" exec "$LOCAL_SIF" cat /opt/prideqc/build-info.txt
"$CONTAINER_BIN" exec "$LOCAL_SIF" prideqc --version
"$CONTAINER_BIN" exec "$LOCAL_SIF" python -c \
  "import importlib.metadata as m; print('python pyOpenMS=' + m.version('pyopenms'))"

if [[ "$STAGE_CACHE" == 1 && -d "$CACHE_ROOT" ]]; then
    echo "==> Stage persistent cache"
    rsync -a "$CACHE_ROOT/" "$LOCAL_CACHE/"
fi

SDRF_BASENAME="$(basename "${SDRF_URL%%\?*}")"
if [[ -r "$DATA_ROOT/sdrf/$SDRF_BASENAME" ]]; then
    echo "==> Stage SDRF from persistent storage"
    cp -a "$DATA_ROOT/sdrf/$SDRF_BASENAME" "$LOCAL_SDRF"
else
    echo "==> Fetch SDRF to node-local scratch"
    "$CONTAINER_BIN" exec --cleanenv --bind "$WORK_ROOT:/work" "$LOCAL_SIF" \
      python -c 'import sys, urllib.request; urllib.request.urlretrieve(sys.argv[1], sys.argv[2])' \
      "$SDRF_URL" /work/input.sdrf.tsv
fi

test -s "$LOCAL_SDRF"

container_env=(
    --env "TMPDIR=/work/tmp"
    --env "XDG_CACHE_HOME=/work/cache"
    --env "OMP_NUM_THREADS=1"
    --env "OPENBLAS_NUM_THREADS=1"
    --env "MKL_NUM_THREADS=1"
    --env "NUMEXPR_NUM_THREADS=1"
    --env "RAYON_NUM_THREADS=1"
)

if [[ "$INPUT_MODE" == pride ]]; then
    echo "==> Download selected PRIDE RAW into scratch and analyze"
    set +e
    "$CONTAINER_BIN" exec --cleanenv --bind "$WORK_ROOT:/work" "${container_env[@]}" "$LOCAL_SIF" \
      prideqc analyze \
        --accession "$PXD" \
        --file "$RAW_FILE" \
        --protocol "$PRIDE_PROTOCOL" \
        --download-dir /work/downloads \
        --sdrf /work/input.sdrf.tsv \
        --sdrf-template "$SDRF_TEMPLATE" \
        --workers "$PRIDEQC_WORKERS" \
        --diagnostics \
        --include-inferred \
        --no-progress \
        --overwrite \
        -o /work/results \
      > >(tee "$RUN_LOG") 2> >(tee -a "$RUN_LOG" >&2)
    status=$?
    set -e
    [[ $status -eq 0 ]] || exit "$status"
elif [[ "$INPUT_MODE" == staged ]]; then
    SOURCE_INPUT="$DATA_ROOT/$PXD/$RAW_FILE"
    [[ -e "$SOURCE_INPUT" ]] || { echo "Staged input missing: $SOURCE_INPUT" >&2; exit 2; }
    echo "==> Stage scientific input from shared storage"
    mkdir -p "$WORK_ROOT/inputs"
    cp -a "$SOURCE_INPUT" "$WORK_ROOT/inputs/$RAW_FILE"
    set +e
    "$CONTAINER_BIN" exec --cleanenv --bind "$WORK_ROOT:/work" "${container_env[@]}" "$LOCAL_SIF" \
      prideqc analyze \
        "/work/inputs/$RAW_FILE" \
        --sdrf /work/input.sdrf.tsv \
        --sdrf-template "$SDRF_TEMPLATE" \
        --workers "$PRIDEQC_WORKERS" \
        --diagnostics \
        --include-inferred \
        --no-progress \
        --overwrite \
        -o /work/results \
      > >(tee "$RUN_LOG") 2> >(tee -a "$RUN_LOG" >&2)
    status=$?
    set -e
    [[ $status -eq 0 ]] || exit "$status"
else
    echo "INPUT_MODE must be 'pride' or 'staged'" >&2
    exit 2
fi

# Results are copied atomically-ish through a temporary sibling directory.
echo "==> Sync completed results to persistent storage"
tmp_final="${FINAL_RESULTS}.partial.${JOB_ID}"
rm -rf "$tmp_final"
mkdir -p "$tmp_final"
rsync -a "$LOCAL_RESULTS/" "$tmp_final/"
cp -a "$RUN_LOG" "$tmp_final/slurm-run.log"
printf '%s\n' \
    "gt_id=$GT_ID" \
    "pxd=$PXD" \
    "raw_file=$RAW_FILE" \
    "sif_sha256=$EXPECTED_SIF_SHA" \
    "job_id=$JOB_ID" \
    "task_id=$TASK_ID" \
    > "$tmp_final/hpc-run-info.txt"
rm -rf "$FINAL_RESULTS"
mv "$tmp_final" "$FINAL_RESULTS"

if [[ "$SYNC_CACHE_BACK" == 1 ]]; then
    echo "==> Update persistent cache after successful work"
    mkdir -p "$CACHE_ROOT"
    rsync -a "$LOCAL_CACHE/" "$CACHE_ROOT/"
fi

ELAPSED=$(( $(date +%s) - START_EPOCH ))
echo "elapsed_seconds=$ELAPSED"
echo "final_result_path=$FINAL_RESULTS"

# Disable failure branch in the EXIT trap; it will still remove scratch.
exit 0
