# syntax=docker/dockerfile:1.7

ARG PYTHON_IMAGE=python:3.12.14-slim-bookworm
ARG RUST_IMAGE=rust:1.98.1-slim-bookworm
ARG UV_IMAGE=ghcr.io/astral-sh/uv:0.12.11
ARG DOTNET_IMAGE=mcr.microsoft.com/dotnet/runtime:8.0.29-bookworm-slim
ARG PMULTIQC_GIT_URL=https://github.com/singjc/pmultiqc.git
ARG PMULTIQC_GIT_REF=fae0c95ff6dac996a7e66138ab0b5992ca7c9b23
ARG OPENMS_GIT_URL=https://github.com/OpenMS/OpenMS.git
ARG OPENMS_GIT_REF=develop
ARG OPENMS_GIT_REQUESTED_REF=develop
ARG OPENMS_MANYLINUX_IMAGE=quay.io/pypa/manylinux_2_34_x86_64
ARG OPENMS_BUILD_JOBS=8
ARG NANOBIND_BACKEND_VERSION=1.0.0

FROM ${UV_IMAGE} AS uv-bin

# ThermoRawFile on Linux uses the OpenMS .NET bridge. Keep the runtime pinned
# and copy it into the final Python image; no SDK is needed at runtime.
FROM ${DOTNET_IMAGE} AS dotnet-runtime

# Build the exact pyOpenMS revision used by the scientific runtime. This is a
# temporary source-build lane for validating OpenMS changes before a release
# wheel containing them exists. build_docker.sh resolves a branch/tag to one
# immutable commit SHA before invoking Docker.
FROM ${OPENMS_MANYLINUX_IMAGE} AS pyopenms-builder
ARG OPENMS_GIT_URL
ARG OPENMS_GIT_REF
ARG OPENMS_GIT_REQUESTED_REF
ARG OPENMS_BUILD_JOBS
ENV CMAKE_BUILD_PARALLEL_LEVEL=${OPENMS_BUILD_JOBS} \
    DOTNET_CLI_TELEMETRY_OPTOUT=1 \
    DOTNET_NOLOGO=1

RUN dnf install -y --disablerepo='epel*' git && dnf clean all
WORKDIR /src
RUN set -eux; \
    git init OpenMS; \
    cd OpenMS; \
    git remote add origin "${OPENMS_GIT_URL}"; \
    # Fetch the immutable SHA already resolved by build_docker.sh. The requested ref
    # is retained only as provenance; fetching a moving branch here creates a race.
    git fetch --depth=1 origin "${OPENMS_GIT_REF}"; \
    git checkout --detach FETCH_HEAD; \
    resolved_sha="$(git rev-parse HEAD)"; \
    test "$resolved_sha" = "${OPENMS_GIT_REF}"; \
    git submodule update --init --depth=1 vcpkg; \
    printf '%s\n' \
      "openms_git_url=${OPENMS_GIT_URL}" \
      "openms_requested_ref=${OPENMS_GIT_REQUESTED_REF}" \
      "openms_resolved_sha=${resolved_sha}" \
      > /tmp/openms-build-info.txt

# Reuse OpenMS's own manylinux build recipe. It installs the C++ libraries into
# /src/OpenMS/install with the same vcpkg/static-dependency policy used by the
# upstream pyOpenMS wheel workflow.
RUN cd /src/OpenMS && bash tools/ci/build_openms_manylinux.sh /src/OpenMS

ENV OpenMS_ROOT=/src/OpenMS/install \
    CMAKE_PREFIX_PATH=/src/OpenMS/vcpkg_installed/x64-linux:/src/OpenMS/install \
    LD_LIBRARY_PATH=/src/OpenMS/install/lib

# Build the abi3 pyOpenMS wheel against that exact OpenMS install, then repair
# it so the prideQC runtime does not need the OpenMS build tree.
RUN set -eux; \
    py=/opt/python/cp312-cp312/bin/python; \
    "$py" -m pip install --no-cache-dir build auditwheel; \
    mkdir -p /tmp/unrepaired /out/wheels; \
    "$py" -m pip wheel \
      --no-deps \
      --wheel-dir /tmp/unrepaired \
      --config-settings=cmake.options.OpenMS_ROOT=/src/OpenMS/install \
      --config-settings=cmake.options.NO_DEPENDENCIES=ON \
      --config-settings=cmake.options.PYOPENMS_PREPARE_WHEEL_REPAIR=ON \
      --config-settings=cmake.options.ARROW_USE_STATIC=ON \
      /src/OpenMS/src/pyOpenMS; \
    auditwheel repair -w /out/wheels /tmp/unrepaired/pyopenms-*.whl; \
    repaired="$(find /out/wheels -maxdepth 1 -type f -name 'pyopenms-*.whl' -print -quit)"; \
    test -n "$repaired"; \
    wheel_sha="$(sha256sum "$repaired" | awk '{print $1}')"; \
    { \
      cat /tmp/openms-build-info.txt; \
      echo "pyopenms_wheel=$(basename "$repaired")"; \
      echo "pyopenms_wheel_sha256=${wheel_sha}"; \
    } > /out/openms-build-info.txt

# prideQC is currently Python-only. This stage is intentionally future-proof:
# if Cargo.toml is added later, release binaries are built with the committed
# Cargo.lock and copied into the runtime image. No Rust toolchain is retained.
FROM ${RUST_IMAGE} AS rust-builder
WORKDIR /src
COPY . .
RUN set -eux; \
    mkdir -p /out/bin; \
    if [ -f Cargo.toml ]; then \
        test -f Cargo.lock; \
        cargo build --release --locked --bins; \
        find target/release -maxdepth 1 -type f -perm /111 -exec cp -a {} /out/bin/ \;; \
        find /out/bin -maxdepth 1 -type f -printf '%f\n' | LC_ALL=C sort > /out/rust-bin-manifest.txt; \
        test -s /out/rust-bin-manifest.txt; \
    else \
        printf '%s\n' none > /out/rust-bin-manifest.txt; \
    fi; \
    rustc --version > /out/rustc-version.txt; \
    cargo --version > /out/cargo-version.txt

FROM ${PYTHON_IMAGE} AS python-builder
ARG NANOBIND_BACKEND_VERSION
ENV UV_PYTHON_DOWNLOADS=0 \
    UV_PROJECT_ENVIRONMENT=/opt/prideqc/.venv \
    UV_LINK_MODE=copy \
    UV_COMPILE_BYTECODE=1 \
    PATH=/opt/prideqc/.venv/bin:/usr/local/bin:/usr/bin:/bin
COPY --from=uv-bin /uv /usr/local/bin/uv
WORKDIR /src/prideqc

# Cache the locked transitive environment separately from project source.
COPY pyproject.toml uv.lock ./
RUN --mount=type=cache,target=/root/.cache/uv \
    uv sync --locked --no-dev --no-install-project --no-install-package pyopenms --no-editable

COPY --from=pyopenms-builder /out/wheels/ /opt/prideqc/wheels/
COPY --from=pyopenms-builder /out/openms-build-info.txt /opt/prideqc/openms-build-info.txt
RUN --mount=type=cache,target=/root/.cache/uv \
    uv pip install --python /opt/prideqc/.venv/bin/python --no-deps "nanobind-backend==${NANOBIND_BACKEND_VERSION}" && \
    uv pip install --python /opt/prideqc/.venv/bin/python --no-deps /opt/prideqc/wheels/*.whl && \
    /opt/prideqc/.venv/bin/python -c \
      "import pyopenms as oms; assert hasattr(oms, 'ThermoRawFile'); assert hasattr(oms, 'BrukerTimsFile'); assert hasattr(oms, 'IDFreeMassErrorEstimator'); assert oms.BrukerTimsFile.Config.CentroidAlgo.HILL_BASED is not None"

COPY . .
RUN --mount=type=cache,target=/root/.cache/uv \
    uv sync --locked --no-dev --no-install-package pyopenms --no-editable && \
    uv pip install --python /opt/prideqc/.venv/bin/python --no-deps \
      "nanobind-backend==${NANOBIND_BACKEND_VERSION}" /opt/prideqc/wheels/*.whl && \
    /opt/prideqc/.venv/bin/python -c \
      "import nanobind_backend, pyopenms as oms; assert hasattr(oms, 'IDFreeMassErrorEstimator')" && \
    /opt/prideqc/.venv/bin/prideqc --version

# Keep reporting isolated from the scientific runtime. The current pmultiqc
# development branch constrains pyOpenMS independently of prideQC, so installing
# both into one environment could downgrade the vendor-reader build used above.
FROM ${PYTHON_IMAGE} AS reporting-builder
ARG PMULTIQC_GIT_URL
ARG PMULTIQC_GIT_REF
ENV UV_PYTHON_DOWNLOADS=0 \
    UV_LINK_MODE=copy \
    UV_COMPILE_BYTECODE=1
COPY --from=uv-bin /uv /usr/local/bin/uv
RUN apt-get update && \
    apt-get install -y --no-install-recommends ca-certificates git && \
    rm -rf /var/lib/apt/lists/*
RUN set -eux; \
    resolved_ref="$(git ls-remote "$PMULTIQC_GIT_URL" "$PMULTIQC_GIT_REF" "refs/heads/$PMULTIQC_GIT_REF" "refs/tags/$PMULTIQC_GIT_REF" | awk 'NR==1 {print $1}')"; \
    if [ -z "$resolved_ref" ] && printf '%s' "$PMULTIQC_GIT_REF" | grep -Eq '^[0-9a-fA-F]{40}$'; then resolved_ref="$PMULTIQC_GIT_REF"; fi; \
    test -n "$resolved_ref"; \
    uv venv --python /usr/local/bin/python /opt/pmultiqc/.venv; \
    uv pip install --python /opt/pmultiqc/.venv/bin/python "git+${PMULTIQC_GIT_URL}@${resolved_ref}"; \
    pmultiqc_version="$(/opt/pmultiqc/.venv/bin/python -c 'import importlib.metadata as m; print(m.version("pmultiqc"))')"; \
    multiqc_version="$(/opt/pmultiqc/.venv/bin/python -c 'import importlib.metadata as m; print(m.version("multiqc"))')"; \
    reporting_pyopenms_version="$(/opt/pmultiqc/.venv/bin/python -c 'import importlib.metadata as m; print(m.version("pyopenms"))')"; \
    mkdir -p /out; \
    { \
      echo "pmultiqc_git_url=${PMULTIQC_GIT_URL}"; \
      echo "pmultiqc_requested_ref=${PMULTIQC_GIT_REF}"; \
      echo "pmultiqc_resolved_sha=${resolved_ref}"; \
      echo "pmultiqc_version=${pmultiqc_version}"; \
      echo "multiqc_version=${multiqc_version}"; \
      echo "reporting_pyopenms_version=${reporting_pyopenms_version}"; \
    } > /out/pmultiqc-build-info.txt; \
    /opt/pmultiqc/.venv/bin/multiqc --version; \
    /opt/pmultiqc/.venv/bin/python -c 'from pmultiqc.modules.mzqc import MzQCModule, parse_mzqc_document; print("pmultiqc mzQC module: OK")'

# Optional quality target. build_docker.sh can build this before runtime.
FROM python-builder AS quality
RUN --mount=type=cache,target=/root/.cache/uv \
    uv sync --locked --group dev --no-install-package pyopenms --no-editable && \
    uv pip install --python /opt/prideqc/.venv/bin/python --no-deps \
      "nanobind-backend==${NANOBIND_BACKEND_VERSION}" /opt/prideqc/wheels/*.whl && \
    /opt/prideqc/.venv/bin/python -c \
      "import nanobind_backend, pyopenms as oms; assert hasattr(oms, 'ThermoRawFile'); assert hasattr(oms, 'BrukerTimsFile'); assert hasattr(oms, 'IDFreeMassErrorEstimator'); assert oms.BrukerTimsFile.Config.CentroidAlgo.HILL_BASED is not None" && \
    /opt/prideqc/.venv/bin/python -m compileall -q src tests && \
    PRIDE_QC_REQUIRE_INTEGRATION=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 /opt/prideqc/.venv/bin/python -m unittest discover -s tests -v && \
    /opt/prideqc/.venv/bin/ruff check . && \
    /opt/prideqc/.venv/bin/mypy

FROM ${PYTHON_IMAGE} AS runtime
ARG PROJECT_VERSION=unknown
ARG VCS_REF=unknown
ARG SOURCE_STATE=unknown
ARG SOURCE_FINGERPRINT=unknown
ARG PYTHON_LOCK_HASH=unknown
ARG CARGO_LOCK_HASH=not-present
ARG BUILD_ID=unknown
ARG OPENMS_GIT_URL=https://github.com/OpenMS/OpenMS.git
ARG OPENMS_GIT_REF=unknown
ARG OPENMS_GIT_REQUESTED_REF=unknown
ARG NANOBIND_BACKEND_VERSION=1.0.0
ARG PYTHON_IMAGE=python:3.12.14-slim-bookworm
ARG RUST_IMAGE=rust:1.98.1-slim-bookworm
ARG UV_IMAGE=ghcr.io/astral-sh/uv:0.12.11
ARG DOTNET_IMAGE=mcr.microsoft.com/dotnet/runtime:8.0.29-bookworm-slim
ARG PMULTIQC_GIT_URL=https://github.com/singjc/pmultiqc.git
ARG PMULTIQC_GIT_REF=fae0c95ff6dac996a7e66138ab0b5992ca7c9b23

RUN apt-get update && \
    apt-get install -y --no-install-recommends ca-certificates libgcc-s1 libgomp1 libgssapi-krb5-2 libicu72 libssl3 libstdc++6 procps tzdata zlib1g && \
    rm -rf /var/lib/apt/lists/* && \
    command -v ps >/dev/null && \
    ps aux >/dev/null && \
    mkdir -p /opt/prideqc/bin /data /results /cache /work && \
    chmod 1777 /data /results /cache /work

COPY --from=python-builder /opt/prideqc/.venv /opt/prideqc/.venv
COPY --from=python-builder /opt/prideqc/openms-build-info.txt /opt/prideqc/openms-build-info.txt
COPY --from=reporting-builder /opt/pmultiqc/.venv /opt/pmultiqc/.venv
COPY --from=reporting-builder /out/pmultiqc-build-info.txt /opt/prideqc/pmultiqc-build-info.txt
COPY --from=dotnet-runtime /usr/share/dotnet /usr/share/dotnet
COPY --from=rust-builder /out/bin/ /opt/prideqc/bin/
COPY --from=rust-builder /out/rust-bin-manifest.txt /opt/prideqc/rust-bin-manifest.txt
COPY --from=rust-builder /out/rustc-version.txt /tmp/rustc-version.txt
COPY --from=rust-builder /out/cargo-version.txt /tmp/cargo-version.txt

ENV PATH=/opt/prideqc/.venv/bin:/opt/pmultiqc/.venv/bin:/opt/prideqc/bin:/usr/share/dotnet:/usr/local/bin:/usr/bin:/bin \
    DOTNET_ROOT=/usr/share/dotnet \
    DOTNET_RUNNING_IN_CONTAINER=true \
    PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    XDG_CACHE_HOME=/cache \
    TMPDIR=/work/tmp \
    OMP_NUM_THREADS=1 \
    OPENBLAS_NUM_THREADS=1 \
    MKL_NUM_THREADS=1 \
    NUMEXPR_NUM_THREADS=1 \
    RAYON_NUM_THREADS=1

RUN set -eux; \
    mkdir -p /work/tmp; \
    PYTHON_VERSION="$(python --version 2>&1)"; \
    PYOPENMS_VERSION="$(python -c 'import importlib.metadata as m; print(m.version("pyopenms"))')"; \
    DOTNET_RUNTIME_VERSION="$(dotnet --list-runtimes | awk '$1 == "Microsoft.NETCore.App" {print $2; exit}')"; \
    test -n "$DOTNET_RUNTIME_VERSION"; \
    { \
      echo "project=prideQC"; \
      echo "project_version=${PROJECT_VERSION}"; \
      echo "vcs_ref=${VCS_REF}"; \
      echo "source_state=${SOURCE_STATE}"; \
      echo "source_fingerprint_sha256=${SOURCE_FINGERPRINT}"; \
      echo "python_version=${PYTHON_VERSION}"; \
      echo "python_package_lock=uv.lock"; \
      echo "python_package_lock_sha256=${PYTHON_LOCK_HASH}"; \
      echo "pyopenms_version=${PYOPENMS_VERSION}"; \
      echo "openms_git_url=${OPENMS_GIT_URL}"; \
      echo "openms_requested_ref=${OPENMS_GIT_REQUESTED_REF}"; \
      echo "openms_build_ref=${OPENMS_GIT_REF}"; \
      echo "nanobind_backend_version=${NANOBIND_BACKEND_VERSION}"; \
      cat /opt/prideqc/openms-build-info.txt; \
      echo "cargo_lock_sha256=${CARGO_LOCK_HASH}"; \
      printf 'rustc_version='; cat /tmp/rustc-version.txt; \
      printf 'cargo_version='; cat /tmp/cargo-version.txt; \
      echo "python_base_image=${PYTHON_IMAGE}"; \
      echo "rust_base_image=${RUST_IMAGE}"; \
      echo "uv_image=${UV_IMAGE}"; \
      echo "dotnet_base_image=${DOTNET_IMAGE}"; \
      echo "dotnet_runtime_version=${DOTNET_RUNTIME_VERSION}"; \
      echo "build_id=${BUILD_ID}"; \
      echo; \
      echo "[optional-reporting]"; \
      cat /opt/prideqc/pmultiqc-build-info.txt; \
    } > /opt/prideqc/build-info.txt; \
    rm -f /tmp/rustc-version.txt /tmp/cargo-version.txt; \
    dotnet --info >/dev/null; \
    prideqc --version; \
    python -c "import pyopenms as oms; assert hasattr(oms, 'ThermoRawFile'); assert hasattr(oms, 'BrukerTimsFile'); assert hasattr(oms, 'IDFreeMassErrorEstimator'); assert oms.BrukerTimsFile.Config.CentroidAlgo.HILL_BASED is not None"; \
    /opt/pmultiqc/.venv/bin/multiqc --version; \
    /opt/pmultiqc/.venv/bin/python -c "from pmultiqc.modules.mzqc import MzQCModule, parse_mzqc_document; print('pmultiqc mzQC module: OK')"; \
    mkdir -p /tmp/mzqc-smoke /tmp/multiqc-smoke; \
    printf '%s\n' '{"mzQC":{"version":"1.0.0","runQualities":[{"metadata":{"label":"container-smoke","inputFiles":[{"name":"smoke.raw","location":"file:///smoke.raw"}]},"qualityMetrics":[{"accession":"MS:4000059","name":"number of MS1 spectra","value":1}]}],"controlledVocabularies":[]}}' > /tmp/mzqc-smoke/smoke.mzQC; \
    /opt/pmultiqc/.venv/bin/multiqc --module mzqc --force --outdir /tmp/multiqc-smoke /tmp/mzqc-smoke; \
    test -s /tmp/multiqc-smoke/multiqc_report.html; \
    rm -rf /tmp/mzqc-smoke /tmp/multiqc-smoke

LABEL org.opencontainers.image.title="prideQC" \
      org.opencontainers.image.description="prideQC scientific runtime with isolated optional pmultiqc/MultiQC reporting"

WORKDIR /work
CMD ["prideqc", "--help"]
