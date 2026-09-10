#!/usr/bin/env bash
set -euo pipefail

: "${PERSIST_ROOT:?Set PERSIST_ROOT, e.g. /nfs/research/juan/DIA/singj/prideQC}"
: "${SIF_REF:?Set SIF_REF, preferably oras://ghcr.io/pride-archive/prideqc-sif:sha-<full-git-sha>}"

if command -v apptainer >/dev/null 2>&1; then
    CONTAINER_BIN=apptainer
elif command -v singularity >/dev/null 2>&1; then
    CONTAINER_BIN=singularity
else
    echo "Neither apptainer nor singularity is available" >&2
    exit 2
fi
command -v sha256sum >/dev/null 2>&1 || { echo "sha256sum is required" >&2; exit 2; }

CONTAINER_DIR="${CONTAINER_DIR:-$PERSIST_ROOT/containers}"
if [[ -n "${SIF_NAME:-}" ]]; then
    sif_name="$SIF_NAME"
else
    ref_tag="${SIF_REF##*:}"
    ref_tag="${ref_tag//[^A-Za-z0-9._-]/_}"
    sif_name="prideqc_${ref_tag}.sif"
fi
[[ "$sif_name" == *.sif ]] || sif_name="${sif_name}.sif"

mkdir -p "$CONTAINER_DIR"
SIF="$CONTAINER_DIR/$sif_name"
TMP_SIF="${SIF}.partial.$$"
trap 'rm -f "$TMP_SIF"' EXIT

echo "container_runtime=$CONTAINER_BIN"
echo "source_ref=$SIF_REF"
echo "destination=$SIF"

"$CONTAINER_BIN" pull --force "$TMP_SIF" "$SIF_REF"

echo "==> Smoke-test pulled SIF"
"$CONTAINER_BIN" exec "$TMP_SIF" prideqc --version
"$CONTAINER_BIN" exec "$TMP_SIF" python -c \
  "import importlib.metadata as m, pyopenms as oms; print('pyOpenMS=' + m.version('pyopenms')); assert hasattr(oms, 'ThermoRawFile'); assert hasattr(oms, 'BrukerTimsFile')"
"$CONTAINER_BIN" exec "$TMP_SIF" dotnet --list-runtimes

mv -f "$TMP_SIF" "$SIF"
trap - EXIT

(
    cd "$CONTAINER_DIR"
    sha256sum "$(basename "$SIF")" > "$(basename "$SIF").sha256"
)

SIF_SHA="$(sha256sum "$SIF" | awk '{print $1}')"
{
    echo "source_ref=$SIF_REF"
    echo "installed_at_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "container_runtime=$CONTAINER_BIN"
    echo "container_runtime_version=$($CONTAINER_BIN --version 2>&1 | head -n 1)"
    echo "sif_sha256=$SIF_SHA"
    echo
    echo "[embedded build-info]"
    "$CONTAINER_BIN" exec "$SIF" cat /opt/prideqc/build-info.txt
} > "${SIF}.meta.txt"

echo "==> Installed"
echo "SIF=$SIF"
echo "SHA256=${SIF}.sha256"
echo "META=${SIF}.meta.txt"
echo
printf 'export SIF=%q\n' "$SIF"
