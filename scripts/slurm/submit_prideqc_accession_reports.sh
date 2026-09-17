#!/usr/bin/env bash
set -euo pipefail

: "${PERSIST_ROOT:?Set PERSIST_ROOT to the persistent project root}"
: "${SIF:?Set SIF to the persistent prideQC .sif path}"
: "${RUN_NAME:?Set RUN_NAME to the completed or submitted prideQC run name}"

RESULT_ROOT="${RESULT_ROOT:-$PERSIST_ROOT/results}"
LOG_ROOT="${LOG_ROOT:-$PERSIST_ROOT/logs}"
SCRATCH_ROOT="${SCRATCH_ROOT:-/tmp}"
MANIFEST_ROOT="${MANIFEST_ROOT:-$PERSIST_ROOT/benchmarks/manifests}"
MANIFEST="${MANIFEST:-$MANIFEST_ROOT/${RUN_NAME}.tsv}"
REPORT_MANIFEST="${REPORT_MANIFEST:-$MANIFEST_ROOT/${RUN_NAME}.accession-reports.tsv}"
REPORT_LAUNCHER="${REPORT_LAUNCHER:-$PERSIST_ROOT/scripts/prideqc_accession_report.sbatch}"
DEPENDENCY_JOB_ID="${DEPENDENCY_JOB_ID:-}"

PARTITION="${REPORT_PARTITION:-${PARTITION:-}}"
REPORT_CPUS="${REPORT_CPUS:-1}"
REPORT_MEMORY="${REPORT_MEMORY:-8G}"
REPORT_TIME_LIMIT="${REPORT_TIME_LIMIT:-01:00:00}"
REPORT_MAX_PARALLEL="${REPORT_MAX_PARALLEL:-2}"
PRIDEQC_REPORT_REQUIRE_COMPLETE="${PRIDEQC_REPORT_REQUIRE_COMPLETE:-1}"

command -v sbatch >/dev/null 2>&1 || { echo "sbatch is required" >&2; exit 2; }
test -r "$MANIFEST" || { echo "File manifest not readable: $MANIFEST" >&2; exit 2; }
test -r "$SIF" || { echo "SIF not readable: $SIF" >&2; exit 2; }
test -r "${SIF}.sha256" || { echo "SIF checksum not readable: ${SIF}.sha256" >&2; exit 2; }
test -r "$REPORT_LAUNCHER" || { echo "Report launcher not readable: $REPORT_LAUNCHER" >&2; exit 2; }
[[ "$REPORT_CPUS" =~ ^[1-9][0-9]*$ ]] || { echo "REPORT_CPUS must be a positive integer" >&2; exit 2; }
[[ "$REPORT_MAX_PARALLEL" =~ ^[1-9][0-9]*$ ]] || { echo "REPORT_MAX_PARALLEL must be a positive integer" >&2; exit 2; }
[[ "$PRIDEQC_REPORT_REQUIRE_COMPLETE" == 0 || "$PRIDEQC_REPORT_REQUIRE_COMPLETE" == 1 ]] || {
    echo "PRIDEQC_REPORT_REQUIRE_COMPLETE must be 0 or 1" >&2
    exit 2
}
if [[ -n "$DEPENDENCY_JOB_ID" ]]; then
    [[ "$DEPENDENCY_JOB_ID" =~ ^[0-9]+$ ]] || {
        echo "DEPENDENCY_JOB_ID must be a numeric Slurm job id" >&2
        exit 2
    }
fi

mkdir -p "$LOG_ROOT" "$MANIFEST_ROOT"

tmp_manifest="${REPORT_MANIFEST}.tmp.$$"
{
    printf 'report_id\tpxd_accession\texpected_file_tasks\n'
    awk -F '\t' '
        NR > 1 && $2 ~ /^PXD[0-9]+$/ { count[$2]++ }
        END {
            for (pxd in count) {
                print pxd "\t" count[pxd]
            }
        }
    ' "$MANIFEST" \
      | LC_ALL=C sort -t $'\t' -k1,1 \
      | awk -F '\t' 'BEGIN {OFS="\t"} {print NR, $1, $2}'
} > "$tmp_manifest"
mv "$tmp_manifest" "$REPORT_MANIFEST"

REPORT_COUNT="$(awk 'END {print (NR > 0 ? NR - 1 : 0)}' "$REPORT_MANIFEST")"
(( REPORT_COUNT > 0 )) || { echo "No accession report tasks generated: $REPORT_MANIFEST" >&2; exit 2; }

printf '%s\n' \
    "file_manifest=$MANIFEST" \
    "report_manifest=$REPORT_MANIFEST" \
    "accessions=$REPORT_COUNT" \
    "run_name=$RUN_NAME" \
    "result_root=$RESULT_ROOT" \
    "sif=$SIF" \
    "dependency_job_id=${DEPENDENCY_JOB_ID:-none}" \
    "report_cpus=$REPORT_CPUS" \
    "report_memory=$REPORT_MEMORY" \
    "report_time_limit=$REPORT_TIME_LIMIT" \
    "report_max_parallel=$REPORT_MAX_PARALLEL" \
    "require_complete=$PRIDEQC_REPORT_REQUIRE_COMPLETE"

echo
column -ts $'\t' "$REPORT_MANIFEST" 2>/dev/null || cat "$REPORT_MANIFEST"

export_list="ALL"
append_export() {
    local key="$1" value="$2"
    [[ "$value" != *,* ]] || {
        echo "$key contains a comma, unsupported by sbatch --export: $value" >&2
        exit 2
    }
    export_list+=",${key}=${value}"
}
append_export PERSIST_ROOT "$PERSIST_ROOT"
append_export SIF "$SIF"
append_export REPORT_MANIFEST "$REPORT_MANIFEST"
append_export RESULT_ROOT "$RESULT_ROOT"
append_export SCRATCH_ROOT "$SCRATCH_ROOT"
append_export RUN_NAME "$RUN_NAME"
append_export PRIDEQC_REPORT_REQUIRE_COMPLETE "$PRIDEQC_REPORT_REQUIRE_COMPLETE"

sbatch_args=(
    --parsable
    --job-name=prideqc-report
    --cpus-per-task="$REPORT_CPUS"
    --mem="$REPORT_MEMORY"
    --time="$REPORT_TIME_LIMIT"
    --array="1-${REPORT_COUNT}%${REPORT_MAX_PARALLEL}"
    --output="$LOG_ROOT/prideqc-report_%A_%a.out"
    --error="$LOG_ROOT/prideqc-report_%A_%a.err"
    --export="$export_list"
)
[[ -n "$PARTITION" ]] && sbatch_args+=(--partition="$PARTITION")
[[ -n "$DEPENDENCY_JOB_ID" ]] && sbatch_args+=(--dependency="afterany:${DEPENDENCY_JOB_ID}")

JOB_ID="$(sbatch "${sbatch_args[@]}" "$REPORT_LAUNCHER")"
JOB_BASE="${JOB_ID%%;*}"

echo "Submitted accession-report Slurm array: $JOB_ID"
echo "Monitor: squeue -j $JOB_BASE"
echo "Accounting: sacct -j $JOB_BASE --format=JobID,State,Elapsed,AllocCPUS,MaxRSS,ExitCode"
echo "Reports: $RESULT_ROOT/$RUN_NAME/<PXD>/_multiqc/multiqc_report.html"
