#!/usr/bin/env bash
set -euo pipefail

: "${PERSIST_ROOT:?Set PERSIST_ROOT}"
: "${SIF:?Set SIF}"
: "${REPORT_MANIFEST:?Set REPORT_MANIFEST}"

RESULT_ROOT="${RESULT_ROOT:-$PERSIST_ROOT/results}"
SCRATCH_ROOT="${SCRATCH_ROOT:-/tmp}"
RUN_NAME="${RUN_NAME:-ground-truth-files}"
PRIDEQC_REPORT_REQUIRE_COMPLETE="${PRIDEQC_REPORT_REQUIRE_COMPLETE:-1}"

TASK_ID="${SLURM_ARRAY_TASK_ID:-${TASK_ID:-1}}"
JOB_ID="${SLURM_JOB_ID:-local-$$}"
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
test -r "$REPORT_MANIFEST" || { echo "Report manifest not readable: $REPORT_MANIFEST" >&2; exit 2; }
test -r "$SIF" || { echo "SIF not readable: $SIF" >&2; exit 2; }
test -r "${SIF}.sha256" || { echo "SIF checksum not readable: ${SIF}.sha256" >&2; exit 2; }
[[ "$TASK_ID" =~ ^[1-9][0-9]*$ ]] || { echo "Invalid report task id: $TASK_ID" >&2; exit 2; }
[[ "$PRIDEQC_REPORT_REQUIRE_COMPLETE" == 0 || "$PRIDEQC_REPORT_REQUIRE_COMPLETE" == 1 ]] || {
    echo "PRIDEQC_REPORT_REQUIRE_COMPLETE must be 0 or 1" >&2
    exit 2
}

ROW="$(awk -F '\t' -v id="$TASK_ID" 'NR>1 && $1==id {print; exit}' "$REPORT_MANIFEST")"
[[ -n "$ROW" ]] || { echo "No report manifest row for task $TASK_ID" >&2; exit 2; }
ROW_PARSE="${ROW//$'\t'/$'\x1f'}"
IFS=$'\x1f' read -r ROW_ID PXD EXPECTED_TASKS <<< "$ROW_PARSE"
[[ "$ROW_ID" == "$TASK_ID" ]] || { echo "Report manifest task mismatch" >&2; exit 2; }
[[ "$PXD" =~ ^PXD[0-9]{6,}$ ]] || { echo "Invalid accession in report manifest: $PXD" >&2; exit 2; }
[[ "$EXPECTED_TASKS" =~ ^[1-9][0-9]*$ ]] || { echo "Invalid expected task count: $EXPECTED_TASKS" >&2; exit 2; }

ACCESSION_ROOT="$RESULT_ROOT/$RUN_NAME/$PXD"
FINAL_REPORT="$ACCESSION_ROOT/_multiqc"
WORK_ROOT="$SCRATCH_ROOT/prideqc_report_${JOB_ID}_${TASK_ID}_${PXD}"
LOCAL_SIF="$WORK_ROOT/$(basename "$SIF")"
LOCAL_REPORT="$WORK_ROOT/report"
LOCAL_TMP="$WORK_ROOT/tmp"
LOCAL_CACHE="$WORK_ROOT/cache"
RUN_LOG="$WORK_ROOT/multiqc.log"

mkdir -p "$WORK_ROOT" "$LOCAL_REPORT" "$LOCAL_TMP" "$LOCAL_CACHE"
trap 'status=$?; rm -rf "$WORK_ROOT"; exit "$status"' EXIT

test -d "$ACCESSION_ROOT" || {
    echo "Accession results not found: $ACCESSION_ROOT" >&2
    exit 2
}

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

mapfile -t MZQC_FILES < <(
    find "$ACCESSION_ROOT" \
      -type f \
      -name '*.mzQC' \
      ! -path "$FINAL_REPORT/*" \
      -printf '%P\n' \
      | LC_ALL=C sort
)
FOUND_MZQC="${#MZQC_FILES[@]}"

printf '%s\n' \
    "job_id=$JOB_ID" \
    "array_task=$TASK_ID" \
    "pxd=$PXD" \
    "run_name=$RUN_NAME" \
    "accession_root=$ACCESSION_ROOT" \
    "expected_file_tasks=$EXPECTED_TASKS" \
    "found_mzqc=$FOUND_MZQC" \
    "require_complete=$PRIDEQC_REPORT_REQUIRE_COMPLETE" \
    "sif=$SIF"

if (( FOUND_MZQC == 0 )); then
    echo "No mzQC files found under $ACCESSION_ROOT" >&2
    exit 1
fi

printf '%s\n' "${MZQC_FILES[@]}" > "$LOCAL_REPORT/mzqc-files.txt"

echo "==> Reporting software provenance"
"$CONTAINER_BIN" exec "$LOCAL_SIF" prideqc --version
"$CONTAINER_BIN" exec "$LOCAL_SIF" /opt/pmultiqc/.venv/bin/multiqc --version
"$CONTAINER_BIN" exec "$LOCAL_SIF" cat /opt/prideqc/build-info.txt

echo "==> Aggregate mzQC for $PXD"
set +e
"$CONTAINER_BIN" exec --cleanenv \
    --bind "$ACCESSION_ROOT:/input:ro" \
    --bind "$WORK_ROOT:/work" \
    --env "TMPDIR=/work/tmp" \
    --env "XDG_CACHE_HOME=/work/cache" \
    "$LOCAL_SIF" \
    /opt/pmultiqc/.venv/bin/multiqc \
      --module mzqc \
      --strict \
      --interactive \
      --force \
      --title "$PXD — mzQC Quality Control" \
      --filename multiqc_report.html \
      --outdir /work/report \
      /input \
    > >(tee "$RUN_LOG") 2> >(tee -a "$RUN_LOG" >&2)
MULTIQC_STATUS=$?
set -e

cp -a "$RUN_LOG" "$LOCAL_REPORT/multiqc.log"
"$CONTAINER_BIN" exec "$LOCAL_SIF" cat /opt/prideqc/build-info.txt > "$LOCAL_REPORT/container-build-info.txt"

REPORT_OUTCOME=complete
REPORT_REASON=ok
if (( FOUND_MZQC != EXPECTED_TASKS )); then
    REPORT_OUTCOME=partial
    REPORT_REASON="expected_${EXPECTED_TASKS}_mzqc_found_${FOUND_MZQC}"
fi
if (( MULTIQC_STATUS != 0 )); then
    REPORT_OUTCOME=failure
    REPORT_REASON="multiqc_exit_${MULTIQC_STATUS}"
fi
if [[ ! -s "$LOCAL_REPORT/multiqc_report.html" ]]; then
    REPORT_OUTCOME=failure
    REPORT_REASON=missing_multiqc_report
fi

{
    echo "report_outcome=$REPORT_OUTCOME"
    echo "report_reason=$REPORT_REASON"
    echo "multiqc_exit_code=$MULTIQC_STATUS"
    echo "job_id=$JOB_ID"
    echo "task_id=$TASK_ID"
    echo "pxd=$PXD"
    echo "run_name=$RUN_NAME"
    echo "expected_file_tasks=$EXPECTED_TASKS"
    echo "found_mzqc=$FOUND_MZQC"
    echo "require_complete=$PRIDEQC_REPORT_REQUIRE_COMPLETE"
    echo "sif=$SIF"
    echo "sif_sha256=$EXPECTED_SIF_SHA"
    echo "started_epoch=$START_EPOCH"
    echo "finished_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
} > "$LOCAL_REPORT/report-info.txt"

echo "==> Persist accession report"
TMP_FINAL="${FINAL_REPORT}.partial.${JOB_ID}.${TASK_ID}"
rm -rf "$TMP_FINAL"
mkdir -p "$TMP_FINAL"
rsync -a "$LOCAL_REPORT/" "$TMP_FINAL/"
rm -rf "$FINAL_REPORT"
mv "$TMP_FINAL" "$FINAL_REPORT"

ELAPSED=$(( $(date +%s) - START_EPOCH ))
echo "elapsed_seconds=$ELAPSED"
echo "report_path=$FINAL_REPORT/multiqc_report.html"
echo "report_outcome=$REPORT_OUTCOME"

if (( MULTIQC_STATUS != 0 )); then
    exit "$MULTIQC_STATUS"
fi
if [[ "$REPORT_OUTCOME" == failure ]]; then
    exit 1
fi
if [[ "$REPORT_OUTCOME" == partial && "$PRIDEQC_REPORT_REQUIRE_COMPLETE" == 1 ]]; then
    echo "Report generated but accession is incomplete: $REPORT_REASON" >&2
    exit 1
fi
exit 0
