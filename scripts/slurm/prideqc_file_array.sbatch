#!/usr/bin/env bash
set -euo pipefail

: "${PERSIST_ROOT:?Set PERSIST_ROOT}"
: "${SIF:?Set SIF}"
: "${MANIFEST:?Set MANIFEST}"

RESULT_ROOT="${RESULT_ROOT:-$PERSIST_ROOT/results}"
SCRATCH_ROOT="${SCRATCH_ROOT:-/tmp}"
RUN_NAME="${RUN_NAME:-ground-truth-files}"
PRIDE_PROTOCOL="${PRIDE_PROTOCOL:-ftp}"
PRIDEQC_DIAGNOSTICS="${PRIDEQC_DIAGNOSTICS:-1}"
PRIDEQC_INCLUDE_INFERRED="${PRIDEQC_INCLUDE_INFERRED:-0}"
PRIDEQC_ESTIMATE_PEAK_TYPE="${PRIDEQC_ESTIMATE_PEAK_TYPE:-0}"

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
test -r "$MANIFEST" || { echo "Manifest not readable: $MANIFEST" >&2; exit 2; }
test -r "$SIF" || { echo "SIF not readable: $SIF" >&2; exit 2; }
test -r "${SIF}.sha256" || { echo "SIF checksum not readable: ${SIF}.sha256" >&2; exit 2; }
[[ "$TASK_ID" =~ ^[1-9][0-9]*$ ]] || { echo "Invalid task id: $TASK_ID" >&2; exit 2; }

ROW="$(awk -F '\t' -v id="$TASK_ID" 'NR>1 && $1==id {print; exit}' "$MANIFEST")"
[[ -n "$ROW" ]] || { echo "No manifest row for task $TASK_ID" >&2; exit 2; }
IFS=$'\t' read -r ROW_ID PXD SDRF_BASENAME SDRF_PATH SDRF_SHA SDRF_URL SDRF_TEMPLATE SDRF_DATA_FILE ARCHIVE_FILE FILE_URI NEEDS_FILE_MAP <<< "$ROW"
[[ "$ROW_ID" == "$TASK_ID" ]] || { echo "Manifest task mismatch" >&2; exit 2; }
[[ "$PXD" =~ ^PXD[0-9]{6,}$ ]] || { echo "Invalid accession in manifest: $PXD" >&2; exit 2; }
test -r "$SDRF_PATH" || { echo "SDRF not readable: $SDRF_PATH" >&2; exit 2; }
[[ "$(sha256sum "$SDRF_PATH" | awk '{print $1}')" == "$SDRF_SHA" ]] || {
    echo "SDRF checksum mismatch: $SDRF_PATH" >&2
    exit 2
}
[[ "$NEEDS_FILE_MAP" == 0 || "$NEEDS_FILE_MAP" == 1 ]] || {
    echo "Invalid needs_file_map value: $NEEDS_FILE_MAP" >&2
    exit 2
}

SAFE_SDRF="${SDRF_BASENAME%.sdrf.tsv}"
SAFE_SDRF="${SAFE_SDRF//[^A-Za-z0-9._-]/_}"
SAFE_FILE="${ARCHIVE_FILE//[^A-Za-z0-9._-]/_}"
[[ -n "$SAFE_FILE" ]] || SAFE_FILE="file"

WORK_ROOT="$SCRATCH_ROOT/prideqc_${JOB_ID}_${TASK_ID}_${PXD}"
LOCAL_SIF="$WORK_ROOT/$(basename "$SIF")"
LOCAL_SDRF="$WORK_ROOT/input.sdrf.tsv"
LOCAL_FILE_MAP="$WORK_ROOT/file-map.json"
LOCAL_DOWNLOAD="$WORK_ROOT/download"
LOCAL_RESULTS="$WORK_ROOT/results"
LOCAL_TMP="$WORK_ROOT/tmp"
LOCAL_CACHE="$WORK_ROOT/cache"
RUN_LOG="$WORK_ROOT/prideqc.log"
FINAL_RESULTS="$RESULT_ROOT/$RUN_NAME/$PXD/$SAFE_SDRF/${TASK_ID}_${SAFE_FILE}"
FAIL_RESULTS="$RESULT_ROOT/$RUN_NAME/_infrastructure_failed/${PXD}_${TASK_ID}_${JOB_ID}"

mkdir -p "$WORK_ROOT" "$LOCAL_DOWNLOAD" "$LOCAL_RESULTS" "$LOCAL_TMP" "$LOCAL_CACHE"

cleanup() {
    local status=$?
    local elapsed=$(( $(date +%s) - START_EPOCH ))
    if [[ $status -ne 0 && ! -e "$WORK_ROOT/.results-synced" ]]; then
        mkdir -p "$FAIL_RESULTS"
        {
            echo "status=$status"
            echo "job_id=$JOB_ID"
            echo "task_id=$TASK_ID"
            echo "pxd=$PXD"
            echo "sdrf_data_file=$SDRF_DATA_FILE"
            echo "archive_file=$ARCHIVE_FILE"
            echo "file_uri=$FILE_URI"
            echo "node=$(hostname)"
            echo "elapsed_seconds=$elapsed"
            echo "work_root=$WORK_ROOT"
            echo
            echo "[memory]"
            head -n 20 /proc/meminfo 2>/dev/null || true
            echo
            echo "[disk]"
            df -h "$WORK_ROOT" || true
            echo
            echo "[files excluding downloaded spectra]"
            find "$WORK_ROOT" -maxdepth 4 -type f ! -path "$LOCAL_DOWNLOAD/*" -printf '%p\t%s bytes\n' 2>/dev/null | sort || true
        } > "$FAIL_RESULTS/failure-diagnostics.txt"
        [[ -f "$RUN_LOG" ]] && cp -a "$RUN_LOG" "$FAIL_RESULTS/" || true
        [[ -d "$LOCAL_RESULTS" ]] && rsync -a "$LOCAL_RESULTS/" "$FAIL_RESULTS/partial-results/" || true
        echo "Infrastructure failure diagnostics: $FAIL_RESULTS" >&2
    fi
    rm -rf "$WORK_ROOT"
    exit "$status"
}
trap cleanup EXIT

printf '%s\n' \
    "job_id=$JOB_ID" \
    "array_task=$TASK_ID" \
    "pxd=$PXD" \
    "sdrf_data_file=$SDRF_DATA_FILE" \
    "archive_file=$ARCHIVE_FILE" \
    "file_uri=$FILE_URI" \
    "needs_file_map=$NEEDS_FILE_MAP" \
    "sdrf=$SDRF_BASENAME" \
    "sdrf_template=$SDRF_TEMPLATE" \
    "node=$(hostname)" \
    "partition=${SLURM_JOB_PARTITION:-unknown}" \
    "cpus=$CPUS" \
    "workers=1" \
    "slurm_mem_per_node=${SLURM_MEM_PER_NODE:-unknown}" \
    "sif=$SIF" \
    "scratch=$WORK_ROOT" \
    "result=$FINAL_RESULTS"

echo "==> Verify and stage SIF"
(
    cd "$(dirname "$SIF")"
    sha256sum -c "$(basename "${SIF}.sha256")"
)
EXPECTED_SIF_SHA="$(awk 'NR==1 {print $1}' "${SIF}.sha256")"
cp -a "$SIF" "$LOCAL_SIF"
[[ "$(sha256sum "$LOCAL_SIF" | awk '{print $1}')" == "$EXPECTED_SIF_SHA" ]] || {
    echo "Staged SIF checksum mismatch" >&2
    exit 2
}

cp -a "$SDRF_PATH" "$LOCAL_SDRF"

echo "==> Software provenance"
"$CONTAINER_BIN" exec "$LOCAL_SIF" prideqc --version
"$CONTAINER_BIN" exec "$LOCAL_SIF" cat /opt/prideqc/build-info.txt

echo "==> Scratch capacity"
df -h "$WORK_ROOT" || true

container_env=(
    --env "TMPDIR=/work/tmp"
    --env "XDG_CACHE_HOME=/work/cache"
    --env "OMP_NUM_THREADS=1"
    --env "OPENBLAS_NUM_THREADS=1"
    --env "MKL_NUM_THREADS=1"
    --env "NUMEXPR_NUM_THREADS=1"
    --env "RAYON_NUM_THREADS=1"
    --env "MALLOC_ARENA_MAX=2"
)
analysis_flags=()
[[ "$PRIDEQC_DIAGNOSTICS" == 1 ]] && analysis_flags+=(--diagnostics)
[[ "$PRIDEQC_INCLUDE_INFERRED" == 1 ]] && analysis_flags+=(--include-inferred)
[[ "$PRIDEQC_ESTIMATE_PEAK_TYPE" == 1 ]] && analysis_flags+=(--estimate-peak-type)

file_map_flags=()
if [[ "$NEEDS_FILE_MAP" == 1 ]]; then
    "$CONTAINER_BIN" exec --cleanenv --bind "$WORK_ROOT:/work" "$LOCAL_SIF" \
        python -c 'import json,sys; print(json.dumps({sys.argv[1]: sys.argv[2]}))' \
        "$SDRF_DATA_FILE" "$ARCHIVE_FILE" > "$LOCAL_FILE_MAP"
    file_map_flags=(--file-map /work/file-map.json)
    echo "SDRF/archive alias: $SDRF_DATA_FILE -> $ARCHIVE_FILE"
fi

# One archive file, one prideQC worker. Parallelism comes from the Slurm array.
echo "==> Analyze $PXD / $ARCHIVE_FILE"
set +e
"$CONTAINER_BIN" exec --cleanenv --bind "$WORK_ROOT:/work" "${container_env[@]}" "$LOCAL_SIF" \
    prideqc analyze \
      --accession "$PXD" \
      --file "$ARCHIVE_FILE" \
      --protocol "$PRIDE_PROTOCOL" \
      --download-dir /work/download \
      --sdrf /work/input.sdrf.tsv \
      --sdrf-template "$SDRF_TEMPLATE" \
      --workers 1 \
      --no-progress \
      --overwrite \
      "${file_map_flags[@]}" \
      "${analysis_flags[@]}" \
      -o /work/results \
    > >(tee "$RUN_LOG") 2> >(tee -a "$RUN_LOG" >&2)
STATUS=$?
set -e

[[ -f "$LOCAL_DOWNLOAD/download-manifest.json" ]] && cp -a "$LOCAL_DOWNLOAD/download-manifest.json" "$LOCAL_RESULTS/download-manifest.json"
cp -a "$LOCAL_SDRF" "$LOCAL_RESULTS/input.sdrf.tsv"
[[ -f "$LOCAL_FILE_MAP" ]] && cp -a "$LOCAL_FILE_MAP" "$LOCAL_RESULTS/file-map.json"
cp -a "$RUN_LOG" "$LOCAL_RESULTS/hpc-analysis.log"
printf 'task_id\tpxd_accession\tsdrf_basename\tsdrf_path\tsdrf_sha256\tsource_sdrf_url\tsdrf_template\tsdrf_data_file\tarchive_file\tfile_uri\tneeds_file_map\n' > "$LOCAL_RESULTS/manifest-row.tsv"
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$ROW_ID" "$PXD" "$SDRF_BASENAME" "$SDRF_PATH" "$SDRF_SHA" "$SDRF_URL" "$SDRF_TEMPLATE" \
    "$SDRF_DATA_FILE" "$ARCHIVE_FILE" "$FILE_URI" "$NEEDS_FILE_MAP" >> "$LOCAL_RESULTS/manifest-row.tsv"
"$CONTAINER_BIN" exec "$LOCAL_SIF" cat /opt/prideqc/build-info.txt > "$LOCAL_RESULTS/container-build-info.txt"
{
    echo "scientific_status=$STATUS"
    echo "job_id=$JOB_ID"
    echo "task_id=$TASK_ID"
    echo "pxd=$PXD"
    echo "sdrf_data_file=$SDRF_DATA_FILE"
    echo "archive_file=$ARCHIVE_FILE"
    echo "file_uri=$FILE_URI"
    echo "needs_file_map=$NEEDS_FILE_MAP"
    echo "sdrf=$SDRF_BASENAME"
    echo "sdrf_sha256=$SDRF_SHA"
    echo "sdrf_url=$SDRF_URL"
    echo "sdrf_template=$SDRF_TEMPLATE"
    echo "node=$(hostname)"
    echo "cpus=$CPUS"
    echo "workers=1"
    echo "pride_protocol=$PRIDE_PROTOCOL"
    echo "sif=$SIF"
    echo "sif_sha256=$EXPECTED_SIF_SHA"
    echo "started_epoch=$START_EPOCH"
    echo "finished_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "scratch_download_bytes=$(du -sb "$LOCAL_DOWNLOAD" 2>/dev/null | awk '{print $1}' || echo unknown)"
} > "$LOCAL_RESULTS/hpc-run-info.txt"

echo "==> Persist compact results only"
TMP_FINAL="${FINAL_RESULTS}.partial.${JOB_ID}.${TASK_ID}"
rm -rf "$TMP_FINAL"
mkdir -p "$TMP_FINAL"
rsync -a "$LOCAL_RESULTS/" "$TMP_FINAL/"
mkdir -p "$(dirname "$FINAL_RESULTS")"
rm -rf "$FINAL_RESULTS"
mv "$TMP_FINAL" "$FINAL_RESULTS"
touch "$WORK_ROOT/.results-synced"

ELAPSED=$(( $(date +%s) - START_EPOCH ))
echo "elapsed_seconds=$ELAPSED"
echo "final_result_path=$FINAL_RESULTS"
if [[ $STATUS -ne 0 ]]; then
    echo "prideQC returned status $STATUS; compact diagnostics were retained at $FINAL_RESULTS" >&2
fi
exit "$STATUS"
