#!/usr/bin/env bash
set -euo pipefail

PROJECT=prideqc
IMAGE_REPOSITORY="${PRIDEQC_IMAGE_REPOSITORY:-prideqc}"
IMAGE_TAG="${PRIDEQC_IMAGE_TAG:-0.2.0-local}"
TARGET="${PRIDEQC_DOCKER_TARGET:-runtime}"
PLATFORM="${PRIDEQC_DOCKER_PLATFORM:-linux/amd64}"
RUN_QUALITY="${PRIDEQC_RUN_QUALITY:-1}"
REQUIRE_CLEAN="${PRIDEQC_REQUIRE_CLEAN:-0}"
PYTHON_IMAGE="${PRIDEQC_PYTHON_IMAGE:-python:3.12.14-slim-bookworm}"
RUST_IMAGE="${PRIDEQC_RUST_IMAGE:-rust:1.98.1-slim-bookworm}"
UV_IMAGE="${PRIDEQC_UV_IMAGE:-ghcr.io/astral-sh/uv:0.12.11}"
DOTNET_IMAGE="${PRIDEQC_DOTNET_IMAGE:-mcr.microsoft.com/dotnet/runtime:8.0.29-bookworm-slim}"
PMULTIQC_GIT_URL="${PRIDEQC_PMULTIQC_GIT_URL:-https://github.com/singjc/pmultiqc.git}"
PMULTIQC_GIT_REF="${PRIDEQC_PMULTIQC_GIT_REF:-fae0c95ff6dac996a7e66138ab0b5992ca7c9b23}"
OPENMS_GIT_URL="${PRIDEQC_OPENMS_GIT_URL:-https://github.com/OpenMS/OpenMS.git}"
OPENMS_GIT_REQUESTED_REF="${PRIDEQC_OPENMS_GIT_REF:-develop}"
OPENMS_BUILD_JOBS="${PRIDEQC_OPENMS_BUILD_JOBS:-8}"
NANOBIND_BACKEND_VERSION="${PRIDEQC_NANOBIND_BACKEND_VERSION:-1.0.0}"

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

command -v docker >/dev/null 2>&1 || { echo "docker is required" >&2; exit 2; }
docker buildx version >/dev/null 2>&1 || { echo "docker buildx is required" >&2; exit 2; }
test -f pyproject.toml
test -f uv.lock
command -v git >/dev/null 2>&1 || { echo "git is required to resolve the OpenMS source ref" >&2; exit 2; }
command -v uv >/dev/null 2>&1 || { echo "uv is required to resolve the Hatch-VCS project version" >&2; exit 2; }

if [[ "$OPENMS_GIT_REQUESTED_REF" =~ ^[0-9a-fA-F]{40}$ ]]; then
    OPENMS_GIT_REF="${OPENMS_GIT_REQUESTED_REF,,}"
else
    OPENMS_GIT_REF="$(
        git ls-remote \
          "$OPENMS_GIT_URL" \
          "$OPENMS_GIT_REQUESTED_REF" \
          "refs/heads/$OPENMS_GIT_REQUESTED_REF" \
          "refs/tags/$OPENMS_GIT_REQUESTED_REF" \
        | awk 'NR == 1 {print tolower($1)}'
    )"
fi
[[ "$OPENMS_GIT_REF" =~ ^[0-9a-f]{40}$ ]] || {
    echo "Could not resolve OpenMS ref '$OPENMS_GIT_REQUESTED_REF' from $OPENMS_GIT_URL" >&2
    exit 2
}
echo "OpenMS source: $OPENMS_GIT_URL $OPENMS_GIT_REQUESTED_REF -> $OPENMS_GIT_REF"

PROJECT_VERSION="$(uvx --from hatch --with hatch-vcs hatch version)"
[[ -n "$PROJECT_VERSION" ]] || {
    echo "Could not resolve prideQC version with uv" >&2
    exit 2
}
PYTHON_LOCK_HASH="$(sha256sum uv.lock | awk '{print $1}')"
if [[ -f Cargo.toml ]]; then
    test -f Cargo.lock || { echo "Cargo.toml exists but Cargo.lock is missing" >&2; exit 2; }
    CARGO_LOCK_HASH="$(sha256sum Cargo.lock | awk '{print $1}')"
else
    CARGO_LOCK_HASH="not-present"
fi

if git rev-parse --is-inside-work-tree >/dev/null 2>&1; then
    VCS_REF="$(git rev-parse HEAD)"
    if [[ -n "$(git status --porcelain --untracked-files=normal)" ]]; then
        SOURCE_STATE=dirty
    else
        SOURCE_STATE=clean
    fi
    SOURCE_FINGERPRINT="$(git ls-files -co --exclude-standard -z | LC_ALL=C sort -z | xargs -0 sha256sum | sha256sum | awk '{print $1}')"
else
    VCS_REF=not-a-git-checkout
    SOURCE_STATE=unknown
    SOURCE_FINGERPRINT="$(find . -type f -not -path './.git/*' -print0 | LC_ALL=C sort -z | xargs -0 sha256sum | sha256sum | awk '{print $1}')"
fi

if [[ "$SOURCE_STATE" != clean ]]; then
    echo "WARNING: building from source state: $SOURCE_STATE" >&2
    if [[ "$REQUIRE_CLEAN" == 1 ]]; then
        echo "PRIDEQC_REQUIRE_CLEAN=1: refusing a non-clean scientific build" >&2
        exit 2
    fi
fi

BUILD_ID="${IMAGE_REPOSITORY}:${IMAGE_TAG}@${VCS_REF:0:12}-${SOURCE_FINGERPRINT:0:12}"
IMAGE="${IMAGE_REPOSITORY}:${IMAGE_TAG}"

build_args=(
    --platform "$PLATFORM"
    -f containers/Dockerfile
    --load
    --build-arg "PROJECT_VERSION=$PROJECT_VERSION"
    --build-arg "VCS_REF=$VCS_REF"
    --build-arg "SOURCE_STATE=$SOURCE_STATE"
    --build-arg "SOURCE_FINGERPRINT=$SOURCE_FINGERPRINT"
    --build-arg "PYTHON_LOCK_HASH=$PYTHON_LOCK_HASH"
    --build-arg "CARGO_LOCK_HASH=$CARGO_LOCK_HASH"
    --build-arg "BUILD_ID=$BUILD_ID"
    --build-arg "OPENMS_GIT_URL=$OPENMS_GIT_URL"
    --build-arg "OPENMS_GIT_REF=$OPENMS_GIT_REF"
    --build-arg "OPENMS_GIT_REQUESTED_REF=$OPENMS_GIT_REQUESTED_REF"
    --build-arg "OPENMS_BUILD_JOBS=$OPENMS_BUILD_JOBS"
    --build-arg "NANOBIND_BACKEND_VERSION=$NANOBIND_BACKEND_VERSION"
    --build-arg "PYTHON_IMAGE=$PYTHON_IMAGE"
    --build-arg "RUST_IMAGE=$RUST_IMAGE"
    --build-arg "UV_IMAGE=$UV_IMAGE"
    --build-arg "DOTNET_IMAGE=$DOTNET_IMAGE"
    --build-arg "PMULTIQC_GIT_URL=$PMULTIQC_GIT_URL"
    --build-arg "PMULTIQC_GIT_REF=$PMULTIQC_GIT_REF"
)

if [[ "${PRIDEQC_DOCKER_NO_CACHE:-0}" == 1 ]]; then
    build_args+=(--no-cache)
fi

if [[ "$RUN_QUALITY" == 1 && "$TARGET" != quality ]]; then
    echo "==> Building quality gate"
    docker buildx build "${build_args[@]}" --target quality -t "${IMAGE}-quality" .
fi

echo "==> Building $IMAGE (target=$TARGET platform=$PLATFORM)"
docker buildx build "${build_args[@]}" --target "$TARGET" -t "$IMAGE" .

echo "==> Runtime smoke tests"
docker run --rm "$IMAGE" prideqc --version
docker run --rm "$IMAGE" python -c \
  "import importlib.metadata as m, pyopenms as oms; print('Python/pyOpenMS:', m.version('pyopenms')); assert hasattr(oms, 'ThermoRawFile'); assert hasattr(oms, 'BrukerTimsFile'); assert hasattr(oms, 'IDFreeMassErrorEstimator'); assert oms.BrukerTimsFile.Config.CentroidAlgo.HILL_BASED is not None"
docker run --rm "$IMAGE" dotnet --info >/dev/null
docker run --rm "$IMAGE" multiqc --version
docker run --rm "$IMAGE" /opt/pmultiqc/.venv/bin/python -c \
  "from pmultiqc.modules.mzqc import MzQCModule, parse_mzqc_document; print('pmultiqc mzQC module: OK')"
docker run --rm "$IMAGE" /opt/prideqc/.venv/bin/python -c \
  "import importlib.metadata as m; print('core pyOpenMS', m.version('pyopenms'))"
docker run --rm "$IMAGE" /opt/pmultiqc/.venv/bin/python -c \
  "import importlib.metadata as m; print('reporting pyOpenMS', m.version('pyopenms'))"
docker run --rm "$IMAGE" cat /opt/prideqc/build-info.txt

echo "==> Docker image ID"
docker image inspect "$IMAGE" --format '{{.Id}}'
echo "Built $IMAGE"
