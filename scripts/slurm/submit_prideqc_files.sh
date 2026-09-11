#!/usr/bin/env bash
set -euo pipefail

: "${PERSIST_ROOT:?Set PERSIST_ROOT to the persistent project root}"
: "${SIF:?Set SIF to the persistent prideQC .sif path}"

GT_FILE="${GT_FILE:-$PERSIST_ROOT/benchmarks/data/prideqc_ground_truth.tsv}"
PREPARE_SCRIPT="${PREPARE_SCRIPT:-$PERSIST_ROOT/scripts/prepare_prideqc_file_manifest.py}"
LAUNCHER="${LAUNCHER:-$PERSIST_ROOT/scripts/prideqc_file_array.sbatch}"
LOG_ROOT="${LOG_ROOT:-$PERSIST_ROOT/logs}"
RESULT_ROOT="${RESULT_ROOT:-$PERSIST_ROOT/results}"
SCRATCH_ROOT="${SCRATCH_ROOT:-/tmp}"
SDRF_ROOT="${SDRF_ROOT:-$PERSIST_ROOT/benchmarks/sdrf}"
MANIFEST_ROOT="${MANIFEST_ROOT:-$PERSIST_ROOT/benchmarks/manifests}"

# Safe defaults for one OpenMS/vendor file per array task.
PARTITION="${PARTITION:-}"
CPUS="${CPUS:-2}"
MEMORY="${MEMORY:-32G}"
TIME_LIMIT="${TIME_LIMIT:-04:00:00}"
MAX_PARALLEL="${MAX_PARALLEL:-4}"
PRIDE_PROTOCOL="${PRIDE_PROTOCOL:-ftp}"
RUN_NAME="${RUN_NAME:-ground-truth-files}"
PRIDEQC_DIAGNOSTICS="${PRIDEQC_DIAGNOSTICS:-1}"
PRIDEQC_INCLUDE_INFERRED="${PRIDEQC_INCLUDE_INFERRED:-0}"
PRIDEQC_ESTIMATE_PEAK_TYPE="${PRIDEQC_ESTIMATE_PEAK_TYPE:-0}"
SDRF_GITHUB_REPO="${SDRF_GITHUB_REPO:-bigbio/sdrf-annotated-datasets}"
SDRF_GITHUB_REF="${SDRF_GITHUB_REF:-main}"

if command -v apptainer >/dev/null 2>&1; then
    CONTAINER_BIN=apptainer
elif command -v singularity >/dev/null 2>&1; then
    CONTAINER_BIN=singularity
else
    echo "Neither apptainer nor singularity is available" >&2
    exit 2
fi
command -v sbatch >/dev/null 2>&1 || { echo "sbatch is required" >&2; exit 2; }
test -r "$GT_FILE" || { echo "Ground-truth TSV not readable: $GT_FILE" >&2; exit 2; }
test -r "$SIF" || { echo "SIF not readable: $SIF" >&2; exit 2; }
test -r "${SIF}.sha256" || { echo "SIF checksum not readable: ${SIF}.sha256" >&2; exit 2; }
test -r "$PREPARE_SCRIPT" || { echo "Manifest helper not readable: $PREPARE_SCRIPT" >&2; exit 2; }
test -r "$LAUNCHER" || { echo "Launcher not readable: $LAUNCHER" >&2; exit 2; }
[[ "$CPUS" =~ ^[1-9][0-9]*$ ]] || { echo "CPUS must be a positive integer" >&2; exit 2; }
[[ "$MAX_PARALLEL" =~ ^[1-9][0-9]*$ ]] || { echo "MAX_PARALLEL must be a positive integer" >&2; exit 2; }

mkdir -p "$LOG_ROOT" "$RESULT_ROOT" "$SDRF_ROOT" "$MANIFEST_ROOT"
MANIFEST="$MANIFEST_ROOT/${RUN_NAME}.tsv"
MANIFEST_META="$MANIFEST_ROOT/${RUN_NAME}.meta.txt"

# Run manifest preparation with the exact Python environment baked into the SIF.
# The persistent root is bound explicitly so this does not depend on site auto-binds.
echo "==> Resolve annotated SDRFs and expand to one file per Slurm task"
"$CONTAINER_BIN" exec --cleanenv --bind "$PERSIST_ROOT:/persist" "$SIF" \
    python "/persist/scripts/$(basename "$PREPARE_SCRIPT")" \
      --gt "/persist/${GT_FILE#"$PERSIST_ROOT"/}" \
      --sdrf-root "/persist/${SDRF_ROOT#"$PERSIST_ROOT"/}/$SDRF_GITHUB_REF" \
      --manifest "/persist/${MANIFEST#"$PERSIST_ROOT"/}" \
      --meta "/persist/${MANIFEST_META#"$PERSIST_ROOT"/}" \
      --manifest-path-base /persist \
      --repo "$SDRF_GITHUB_REPO" \
      --ref "$SDRF_GITHUB_REF"

TASK_COUNT="$(awk 'END {print (NR > 0 ? NR - 1 : 0)}' "$MANIFEST")"
(( TASK_COUNT > 0 )) || { echo "No file tasks generated: $MANIFEST" >&2; exit 2; }
TASK_RANGE="${TASK_RANGE:-1-${TASK_COUNT}}"
[[ "$TASK_RANGE" =~ ^[0-9,:-]+$ ]] || { echo "Invalid TASK_RANGE: $TASK_RANGE" >&2; exit 2; }

printf '%s\n' \
    "manifest=$MANIFEST" \
    "file_tasks=$TASK_COUNT" \
    "task_range=$TASK_RANGE" \
    "partition=${PARTITION:-scheduler-default}" \
    "cpus_per_file=$CPUS" \
    "memory_per_file=$MEMORY" \
    "time_limit_per_file=$TIME_LIMIT" \
    "max_parallel=$MAX_PARALLEL" \
    "pride_protocol=$PRIDE_PROTOCOL" \
    "run_name=$RUN_NAME" \
    "sif=$SIF" \
    "result_root=$RESULT_ROOT"

echo
echo "Tasks per accession:"
awk -F '\t' 'NR>1 {n[$2]++} END {for (pxd in n) print pxd, n[pxd]}' "$MANIFEST" | sort

echo
export_list="ALL"
append_export() {
    local key="$1" value="$2"
    [[ "$value" != *,* ]] || { echo "$key contains a comma, unsupported by sbatch --export: $value" >&2; exit 2; }
    export_list+=",${key}=${value}"
}
append_export PERSIST_ROOT "$PERSIST_ROOT"
append_export SIF "$SIF"
append_export MANIFEST "$MANIFEST"
append_export RESULT_ROOT "$RESULT_ROOT"
append_export SCRATCH_ROOT "$SCRATCH_ROOT"
append_export PRIDE_PROTOCOL "$PRIDE_PROTOCOL"
append_export RUN_NAME "$RUN_NAME"
append_export PRIDEQC_DIAGNOSTICS "$PRIDEQC_DIAGNOSTICS"
append_export PRIDEQC_INCLUDE_INFERRED "$PRIDEQC_INCLUDE_INFERRED"
append_export PRIDEQC_ESTIMATE_PEAK_TYPE "$PRIDEQC_ESTIMATE_PEAK_TYPE"

sbatch_args=(
    --parsable
    --job-name=prideqc-file
    --cpus-per-task="$CPUS"
    --mem="$MEMORY"
    --time="$TIME_LIMIT"
    --array="${TASK_RANGE}%${MAX_PARALLEL}"
    --output="$LOG_ROOT/prideqc-file_%A_%a.out"
    --error="$LOG_ROOT/prideqc-file_%A_%a.err"
    --export="$export_list"
)
[[ -n "$PARTITION" ]] && sbatch_args+=(--partition="$PARTITION")
JOB_ID="$(sbatch "${sbatch_args[@]}" "$LAUNCHER")"

JOB_BASE="${JOB_ID%%;*}"
echo "Submitted file-level Slurm array: $JOB_ID"
echo "Monitor: squeue -j $JOB_BASE"
echo "Accounting: sacct -j $JOB_BASE --format=JobID,State,Elapsed,AllocCPUS,MaxRSS,ExitCode"
echo "Results: $RESULT_ROOT/$RUN_NAME"
