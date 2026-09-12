#!/usr/bin/env bash
set -euo pipefail

: "${PERSIST_ROOT:?Set PERSIST_ROOT, e.g. /nfs/research/juan/DIA/singj/prideQC}"

if [[ $# -ne 1 ]]; then
    echo "Usage: $0 <full-git-sha>" >&2
    echo "Example: $0 489757a400161146f77234a02a8fac137e03cc5a" >&2
    exit 2
fi

GIT_SHA="$1"

if [[ ! "$GIT_SHA" =~ ^[0-9a-fA-F]{40}$ ]]; then
    echo "Git SHA must be exactly 40 hexadecimal characters: $GIT_SHA" >&2
    exit 2
fi

command -v singularity >/dev/null 2>&1 || {
    echo "singularity is required on Codon" >&2
    exit 2
}
command -v sha256sum >/dev/null 2>&1 || {
    echo "sha256sum is required" >&2
    exit 2
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export SIF_REF="oras://ghcr.io/pride-archive/prideqc-sif:sha-${GIT_SHA,,}"

exec "$SCRIPT_DIR/install_prideqc_sif.sh"
