#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 <docker-image:tag> <output.sif>" >&2
    exit 2
fi
IMAGE="$1"
SIF="$2"

if command -v apptainer >/dev/null 2>&1; then
    CONTAINER_BIN=apptainer
elif command -v singularity >/dev/null 2>&1; then
    CONTAINER_BIN=singularity
else
    echo "Neither apptainer nor singularity is available" >&2
    exit 2
fi
command -v docker >/dev/null 2>&1 || { echo "docker is required to access the canonical local image" >&2; exit 2; }

mkdir -p "$(dirname "$SIF")"
SIF_DIR="$(cd "$(dirname "$SIF")" && pwd)"
SIF_BASE="$(basename "$SIF")"
SIF="$SIF_DIR/$SIF_BASE"
SHA_FILE="${SIF}.sha256"
META_FILE="${SIF}.meta.txt"

tmpdir="$(mktemp -d)"
cleanup() { rm -rf "$tmpdir"; }
trap cleanup EXIT

rm -f "$SIF" "$SHA_FILE" "$META_FILE"
echo "==> Converting canonical Docker image $IMAGE with $CONTAINER_BIN"
if ! "$CONTAINER_BIN" build --force "$SIF" "docker-daemon://$IMAGE"; then
    echo "docker-daemon transport failed; falling back to docker-archive" >&2
    archive="$tmpdir/image.tar"
    docker save "$IMAGE" -o "$archive"
    "$CONTAINER_BIN" build --force "$SIF" "docker-archive://$archive"
fi

(
    cd "$SIF_DIR"
    sha256sum "$SIF_BASE" > "${SIF_BASE}.sha256"
)

# The checksum is deliberately checked before any smoke test.
(
    cd "$SIF_DIR"
    sha256sum -c "${SIF_BASE}.sha256"
)

echo "==> SIF smoke tests"
"$CONTAINER_BIN" exec "$SIF" prideqc --version
"$CONTAINER_BIN" exec "$SIF" python -c \
  "import importlib.metadata as m, pyopenms as oms; print('pyOpenMS', m.version('pyopenms')); assert hasattr(oms, 'ThermoRawFile'); assert hasattr(oms, 'BrukerTimsFile')"
"$CONTAINER_BIN" exec "$SIF" dotnet --info >/dev/null
"$CONTAINER_BIN" exec "$SIF" test -r /opt/prideqc/build-info.txt

rust_manifest="$("$CONTAINER_BIN" exec "$SIF" cat /opt/prideqc/rust-bin-manifest.txt)"
if [[ "$rust_manifest" != none ]]; then
    while IFS= read -r bin; do
        [[ -n "$bin" ]] || continue
        echo "Rust smoke: $bin"
        "$CONTAINER_BIN" exec "$SIF" test -x "/opt/prideqc/bin/$bin"
        if ! "$CONTAINER_BIN" exec "$SIF" "/opt/prideqc/bin/$bin" --version >/dev/null 2>&1; then
            "$CONTAINER_BIN" exec "$SIF" "/opt/prideqc/bin/$bin" --help >/dev/null
        fi
    done <<< "$rust_manifest"
else
    echo "No Rust binaries in this source revision (current prideQC is Python-only)."
fi

{
    echo "docker_image=$IMAGE"
    echo "docker_image_id=$(docker image inspect "$IMAGE" --format '{{.Id}}')"
    echo "sif_sha256=$(sha256sum "$SIF" | awk '{print $1}')"
    echo "converter=$CONTAINER_BIN"
    printf 'converter_version='; "$CONTAINER_BIN" --version | head -n 1
    echo "converted_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo
    echo "[embedded build-info]"
    "$CONTAINER_BIN" exec "$SIF" cat /opt/prideqc/build-info.txt
} > "$META_FILE"

echo "Created:"
echo "  $SIF"
echo "  $SHA_FILE"
echo "  $META_FILE"
echo "Verify later with: (cd '$SIF_DIR' && sha256sum -c '${SIF_BASE}.sha256')"
