#!/usr/bin/env bash
set -euo pipefail

JOB_ID="${1:?Usage: inspect_prideqc_array.sh JOB_ID [LOG_ROOT]}"
LOG_ROOT="${2:-${PERSIST_ROOT:-.}/logs}"
command -v sacct >/dev/null 2>&1 || { echo "sacct is required" >&2; exit 2; }

echo "==> Accounting"
sacct -j "$JOB_ID" --units=G --format=JobID,State,Elapsed,AllocCPUS,ReqMem,MaxRSS,ExitCode

echo
echo "==> Error summaries"
shopt -s nullglob
files=("$LOG_ROOT"/prideqc-file_"$JOB_ID"_*.err "$LOG_ROOT"/prideqc-pxd_"$JOB_ID"_*.err)
if (( ${#files[@]} == 0 )); then
    echo "No matching .err files under $LOG_ROOT"
    exit 0
fi
for file in "${files[@]}"; do
    echo "--- $file ---"
    grep -Ei 'out.of.memory|oom|killed|traceback|error|failed|exception|prideqc:' "$file" | tail -n 25 || tail -n 25 "$file"
done
