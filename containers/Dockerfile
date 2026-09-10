# syntax=docker/dockerfile:1.7

ARG PYTHON_IMAGE=python:3.12.14-slim-bookworm
ARG RUST_IMAGE=rust:1.98.1-slim-bookworm
ARG UV_IMAGE=ghcr.io/astral-sh/uv:0.12.11
ARG DOTNET_IMAGE=mcr.microsoft.com/dotnet/runtime:8.0.29-bookworm-slim

FROM ${UV_IMAGE} AS uv-bin

# ThermoRawFile on Linux uses the OpenMS .NET bridge. Keep the runtime pinned
# and copy it into the final Python image; no SDK is needed at runtime.
FROM ${DOTNET_IMAGE} AS dotnet-runtime

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
ENV UV_PYTHON_DOWNLOADS=0 \
    UV_PROJECT_ENVIRONMENT=/opt/prideqc/.venv \
    UV_LINK_MODE=copy \
    UV_COMPILE_BYTECODE=1 \
    PATH=/opt/prideqc/.venv/bin:/usr/local/bin:/usr/bin:/bin
COPY --from=uv-bin /uv /usr/local/bin/uv
WORKDIR /src/prideqc

# Cache the locked transitive environment separately from project source.
COPY pyproject.toml uv.lock ./
COPY containers/pyopenms.requirements.txt /tmp/pyopenms.requirements.txt
RUN --mount=type=cache,target=/root/.cache/uv \
    uv sync --locked --no-dev --no-install-project --no-install-package pyopenms --no-editable

COPY . .
RUN --mount=type=cache,target=/root/.cache/uv \
    uv sync --locked --no-dev --no-install-package pyopenms --no-editable && \
    uv pip install --python /opt/prideqc/.venv/bin/python --no-deps --require-hashes \
        -r /tmp/pyopenms.requirements.txt && \
    /opt/prideqc/.venv/bin/python -c \
      "import importlib.metadata as m, pyopenms as oms; v=m.version('pyopenms'); print('pyOpenMS', v); assert v == '3.6.0.dev20260910'; assert hasattr(oms, 'ThermoRawFile'); assert hasattr(oms, 'BrukerTimsFile')" && \
    /opt/prideqc/.venv/bin/prideqc --version

# Optional quality target. build_docker.sh can build this before runtime.
FROM python-builder AS quality
RUN --mount=type=cache,target=/root/.cache/uv \
    uv sync --locked --group dev --no-install-package pyopenms --no-editable && \
    uv pip install --python /opt/prideqc/.venv/bin/python --no-deps --require-hashes \
        -r /tmp/pyopenms.requirements.txt && \
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
ARG PYOPENMS_PIN_HASH=unknown
ARG PYTHON_IMAGE=python:3.12.14-slim-bookworm
ARG RUST_IMAGE=rust:1.98.1-slim-bookworm
ARG UV_IMAGE=ghcr.io/astral-sh/uv:0.12.11
ARG DOTNET_IMAGE=mcr.microsoft.com/dotnet/runtime:8.0.29-bookworm-slim

RUN apt-get update && \
    apt-get install -y --no-install-recommends ca-certificates libgcc-s1 libgomp1 libgssapi-krb5-2 libicu72 libssl3 libstdc++6 tzdata zlib1g && \
    rm -rf /var/lib/apt/lists/* && \
    mkdir -p /opt/prideqc/bin /data /results /cache /work && \
    chmod 1777 /data /results /cache /work

COPY --from=python-builder /opt/prideqc/.venv /opt/prideqc/.venv
COPY --from=python-builder /tmp/pyopenms.requirements.txt /opt/prideqc/pyopenms.requirements.txt
COPY --from=dotnet-runtime /usr/share/dotnet /usr/share/dotnet
COPY --from=rust-builder /out/bin/ /opt/prideqc/bin/
COPY --from=rust-builder /out/rust-bin-manifest.txt /opt/prideqc/rust-bin-manifest.txt
COPY --from=rust-builder /out/rustc-version.txt /tmp/rustc-version.txt
COPY --from=rust-builder /out/cargo-version.txt /tmp/cargo-version.txt

ENV PATH=/opt/prideqc/.venv/bin:/opt/prideqc/bin:/usr/share/dotnet:/usr/local/bin:/usr/bin:/bin \
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
      echo "pyopenms_pin_file_sha256=${PYOPENMS_PIN_HASH}"; \
      printf 'pyopenms_requirement='; cat /opt/prideqc/pyopenms.requirements.txt; \
      echo "cargo_lock_sha256=${CARGO_LOCK_HASH}"; \
      printf 'rustc_version='; cat /tmp/rustc-version.txt; \
      printf 'cargo_version='; cat /tmp/cargo-version.txt; \
      echo "python_base_image=${PYTHON_IMAGE}"; \
      echo "rust_base_image=${RUST_IMAGE}"; \
      echo "uv_image=${UV_IMAGE}"; \
      echo "dotnet_base_image=${DOTNET_IMAGE}"; \
      echo "dotnet_runtime_version=${DOTNET_RUNTIME_VERSION}"; \
      echo "build_id=${BUILD_ID}"; \
    } > /opt/prideqc/build-info.txt; \
    rm -f /tmp/rustc-version.txt /tmp/cargo-version.txt; \
    dotnet --info >/dev/null; \
    prideqc --version; \
    python -c "import pyopenms as oms; assert hasattr(oms, 'ThermoRawFile'); assert hasattr(oms, 'BrukerTimsFile')"

LABEL org.opencontainers.image.title="prideQC" \
      org.opencontainers.image.description="prideQC scientific runtime with pinned pyOpenMS vendor readers"

WORKDIR /work
CMD ["prideqc", "--help"]
