"""Managed local llama.cpp runtime and default GGUF model for prideQC adjudication.

The default installation path is entirely user-space. prideQC downloads pinned,
prebuilt llama.cpp binaries and a pinned quantized model, verifies their SHA-256
checksums, and caches them without requiring a compiler toolchain.
"""

from __future__ import annotations

import hashlib
import json
import os
import platform
import shutil
import stat
import subprocess
import sys
import tarfile
import urllib.request
import zipfile
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path
from typing import Any, BinaryIO

DEFAULT_LLAMA_BUILD = "b11046"
DEFAULT_MODEL_REPOSITORY = "Qwen/Qwen3-4B-GGUF"
DEFAULT_MODEL_REVISION = "a9a60d009fa7ff9606305047c2bf77ac25dbec49"
DEFAULT_MODEL_FILE = "Qwen3-4B-Q4_K_M.gguf"
DEFAULT_MODEL_SHA256 = "7485fe6f11af29433bc51cab58009521f205840f5b4ae3a32fa7f92e8534fdf5"
DEFAULT_MODEL_SIZE = 2_497_280_256
DEFAULT_MODEL_LICENSE = "Apache-2.0"
DEFAULT_MODEL_URL = (
    "https://huggingface.co/Qwen/Qwen3-4B-GGUF/resolve/"
    f"{DEFAULT_MODEL_REVISION}/{DEFAULT_MODEL_FILE}"
)
CACHE_ENVIRONMENT_VARIABLE = "PRIDEQC_LLM_CACHE"


@dataclass(frozen=True)
class RuntimeAsset:
    """One pinned llama.cpp prebuilt runtime archive."""

    system: str
    machine: str
    archive_name: str
    sha256: str
    server_name: str

    @property
    def url(self) -> str:
        return (
            "https://github.com/ggml-org/llama.cpp/releases/download/"
            f"{DEFAULT_LLAMA_BUILD}/{self.archive_name}"
        )


# SHA-256 values are from official ggml-org/llama.cpp GitHub attestations for b11046.
_RUNTIME_ASSETS: dict[tuple[str, str], RuntimeAsset] = {
    ("Linux", "x86_64"): RuntimeAsset(
        system="Linux",
        machine="x86_64",
        archive_name="llama-b11046-bin-ubuntu-x64.tar.gz",
        sha256="ca14dec04b4c5725b6373681cc3685c5dbb914301f7f8be6c89649c0a68734a9",
        server_name="llama-server",
    ),
    ("Linux", "arm64"): RuntimeAsset(
        system="Linux",
        machine="arm64",
        archive_name="llama-b11046-bin-ubuntu-arm64.tar.gz",
        sha256="7920918ce4d5a61ff8b3654102060fe0fc4bf695057d79949e3a4e3b1fb50663",
        server_name="llama-server",
    ),
    ("Darwin", "arm64"): RuntimeAsset(
        system="Darwin",
        machine="arm64",
        archive_name="llama-b11046-bin-macos-arm64.tar.gz",
        sha256="96623092033e83545cd92316ff63d0378a976a4979735f80a1934cc128efb0ad",
        server_name="llama-server",
    ),
    ("Windows", "x86_64"): RuntimeAsset(
        system="Windows",
        machine="x86_64",
        archive_name="llama-b11046-bin-win-cpu-x64.zip",
        sha256="dece5ebc80fbc48063540aa4c48fdc91c59328b7ea44c2bdf6e85ec4383b5fc7",
        server_name="llama-server.exe",
    ),
}


def _normalise_machine(value: str) -> str:
    lowered = value.casefold()
    if lowered in {"x86_64", "amd64"}:
        return "x86_64"
    if lowered in {"aarch64", "arm64"}:
        return "arm64"
    return lowered


def runtime_asset(
    *,
    system: str | None = None,
    machine: str | None = None,
) -> RuntimeAsset:
    """Return the pinned CPU runtime asset for this platform."""
    detected_system = system or platform.system()
    detected_machine = _normalise_machine(machine or platform.machine())
    asset = _RUNTIME_ASSETS.get((detected_system, detected_machine))
    if asset is None:
        raise RuntimeError(
            "No managed llama.cpp CPU binary is currently pinned for "
            f"{detected_system}/{detected_machine}. Use a supported Linux x86_64/arm64, "
            "macOS Apple Silicon, or Windows x86_64 host, or supply a custom runtime later."
        )
    return asset


def runtime_environment(server: Path) -> dict[str, str]:
    """Return an environment that can resolve shared libraries beside llama-server."""
    environment = os.environ.copy()
    system = platform.system()
    if system == "Windows":
        return environment
    variable = "DYLD_LIBRARY_PATH" if system == "Darwin" else "LD_LIBRARY_PATH"
    existing = environment.get(variable)
    environment[variable] = (
        str(server.parent) if not existing else str(server.parent) + os.pathsep + existing
    )
    return environment


def _check_server_executable(server: Path) -> None:
    """Fail setup early when the pinned prebuilt runtime cannot execute on this host."""
    try:
        result = subprocess.run(
            [str(server), "--help"],
            cwd=server.parent,
            env=runtime_environment(server),
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,
            text=True,
            timeout=30,
            check=False,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        raise RuntimeError(f"Unable to execute managed llama-server: {exc}") from exc
    if result.returncode != 0:
        detail = result.stderr.strip()[-2000:]
        raise RuntimeError(
            "Managed llama-server failed its runtime smoke check"
            + (f": {detail}" if detail else "")
        )


def default_cache_dir() -> Path:
    """Return prideQC's user-space LLM cache directory."""
    explicit = os.environ.get(CACHE_ENVIRONMENT_VARIABLE)
    if explicit:
        return Path(explicit).expanduser()
    system = platform.system()
    if system == "Windows" and os.environ.get("LOCALAPPDATA"):
        return Path(os.environ["LOCALAPPDATA"]) / "prideqc" / "llm"
    if system == "Darwin":
        return Path.home() / "Library" / "Caches" / "prideqc" / "llm"
    xdg = os.environ.get("XDG_CACHE_HOME")
    root = Path(xdg).expanduser() if xdg else Path.home() / ".cache"
    return root / "prideqc" / "llm"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def _parse_content_range(value: str | None) -> tuple[int, int, int] | None:
    if not value:
        return None
    units, _, remainder = value.partition(" ")
    if units.casefold() != "bytes":
        return None
    byte_range, _, total_text = remainder.partition("/")
    start_text, separator, end_text = byte_range.partition("-")
    if separator != "-" or not start_text or not end_text or not total_text:
        return None
    if total_text == "*":
        return None
    try:
        start = int(start_text)
        end = int(end_text)
        total = int(total_text)
    except ValueError:
        return None
    if start < 0 or end < start or total <= end:
        return None
    return start, end, total


def _stream_download(url: str, output: BinaryIO) -> int:
    """Stream one ordinary HTTP object into ``output``."""
    request = urllib.request.Request(
        url,
        headers={
            "User-Agent": "prideQC-local-llm",
            "Accept-Encoding": "identity",
        },
    )
    count = 0
    with urllib.request.urlopen(request, timeout=120) as response:
        while chunk := response.read(1024 * 1024):
            output.write(chunk)
            count += len(chunk)
    return count


def _download_verified(
    url: str,
    destination: Path,
    sha256: str,
    *,
    expected_size: int | None = None,
) -> Path:
    """Download a normal-sized asset atomically and verify its checksum."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(destination.name + ".part")
    temporary.unlink(missing_ok=True)
    try:
        with temporary.open("wb") as handle:
            size = _stream_download(url, handle)
        if expected_size is not None and size != expected_size:
            raise RuntimeError(
                f"Downloaded size mismatch for {destination.name}: {size} != {expected_size}"
            )
        observed = _sha256(temporary)
        if observed != sha256:
            raise RuntimeError(
                f"SHA-256 mismatch for {destination.name}: {observed} != {sha256}"
            )
        temporary.replace(destination)
    except Exception:
        temporary.unlink(missing_ok=True)
        raise
    return destination


def _download_verified_ranges(
    url: str,
    destination: Path,
    sha256: str,
    *,
    expected_size: int,
    chunk_size: int = 64 * 1024 * 1024,
    progress: Callable[[int, int], None] | None = None,
) -> Path:
    """Download a large Hub file with explicit bounded HTTP ranges and resume support.

    Hugging Face's ``resolve`` endpoint accepts standard byte-range requests for the
    legacy-compatible download path even when the repository is Xet-backed.  Requesting
    bounded chunks avoids relying on redirect/CDN partial-response behaviour and leaves
    the ``.part`` file in place after transport failures so a later setup can resume.
    Integrity is still enforced over the complete file before installation.
    """
    if expected_size <= 0:
        raise ValueError("expected_size must be positive")
    if chunk_size <= 0:
        raise ValueError("chunk_size must be positive")

    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(destination.name + ".part")
    if temporary.exists() and temporary.stat().st_size > expected_size:
        temporary.unlink()

    current = temporary.stat().st_size if temporary.exists() else 0
    if current == expected_size:
        observed = _sha256(temporary)
        if observed == sha256:
            temporary.replace(destination)
            if progress is not None:
                progress(expected_size, expected_size)
            return destination
        temporary.unlink()
        current = 0

    mode = "ab" if current else "wb"
    try:
        with temporary.open(mode) as handle:
            while current < expected_size:
                requested_end = min(current + chunk_size - 1, expected_size - 1)
                request = urllib.request.Request(
                    url,
                    headers={
                        "User-Agent": "prideQC-local-llm",
                        "Accept-Encoding": "identity",
                        "Range": f"bytes={current}-{requested_end}",
                    },
                )
                with urllib.request.urlopen(request, timeout=120) as response:
                    status = getattr(response, "status", None)
                    content_range = _parse_content_range(response.headers.get("Content-Range"))
                    if status != 206 or content_range is None:
                        raise RuntimeError(
                            "Hugging Face did not honor the requested byte range for the "
                            "local LLM model download"
                        )
                    start, end, total = content_range
                    if start != current or end != requested_end or total != expected_size:
                        raise RuntimeError(
                            "Unexpected Content-Range while downloading local LLM model: "
                            f"requested bytes={current}-{requested_end}, received "
                            f"bytes={start}-{end}/{total}"
                        )

                    segment_size = 0
                    while chunk := response.read(1024 * 1024):
                        handle.write(chunk)
                        segment_size += len(chunk)
                    expected_segment_size = end - start + 1
                    if segment_size != expected_segment_size:
                        raise RuntimeError(
                            "Incomplete byte range while downloading local LLM model: "
                            f"received {segment_size} != {expected_segment_size} bytes"
                        )
                    handle.flush()
                    current += segment_size
                    if progress is not None:
                        progress(current, expected_size)
    except Exception:
        # Keep a valid prefix so the next setup invocation can resume rather than
        # re-downloading multiple gigabytes. Final integrity is always checked below.
        raise

    if temporary.stat().st_size != expected_size:
        raise RuntimeError(
            f"Downloaded size mismatch for {destination.name}: "
            f"{temporary.stat().st_size} != {expected_size}"
        )
    observed = _sha256(temporary)
    if observed != sha256:
        temporary.unlink(missing_ok=True)
        raise RuntimeError(
            f"SHA-256 mismatch for {destination.name}: {observed} != {sha256}"
        )
    temporary.replace(destination)
    return destination


def _model_progress(downloaded: int, total: int) -> None:
    percent = downloaded * 100.0 / total
    downloaded_gib = downloaded / (1024**3)
    total_gib = total / (1024**3)
    print(
        f"\rDownloading {DEFAULT_MODEL_FILE}: "
        f"{downloaded_gib:.2f}/{total_gib:.2f} GiB ({percent:5.1f}%)",
        end="",
        file=sys.stderr,
        flush=True,
    )
    if downloaded >= total:
        print(file=sys.stderr)

def _safe_member_path(root: Path, name: str) -> Path:
    root_resolved = root.resolve()
    target = (root / name).resolve()
    if target != root_resolved and root_resolved not in target.parents:
        raise RuntimeError(f"Unsafe path in llama.cpp archive: {name}")
    return target


def _safe_extract_tar(archive_path: Path, destination: Path) -> None:
    with tarfile.open(archive_path, "r:*") as archive:
        for member in archive.getmembers():
            _safe_member_path(destination, member.name)
            if member.issym():
                relative_target = str(Path(member.name).parent / member.linkname)
                _safe_member_path(destination, relative_target)
            elif member.islnk():
                _safe_member_path(destination, member.linkname)
        if sys.version_info >= (3, 12):
            archive.extractall(destination, filter="fully_trusted")
        else:
            archive.extractall(destination)


def _safe_extract_zip(archive_path: Path, destination: Path) -> None:
    with zipfile.ZipFile(archive_path) as archive:
        for member in archive.infolist():
            _safe_member_path(destination, member.filename)
        archive.extractall(destination)


def _extract_runtime(archive_path: Path, destination: Path) -> None:
    temporary = destination.with_name(destination.name + ".extracting")
    shutil.rmtree(temporary, ignore_errors=True)
    temporary.mkdir(parents=True, exist_ok=True)
    try:
        if archive_path.name.endswith(".zip"):
            _safe_extract_zip(archive_path, temporary)
        else:
            _safe_extract_tar(archive_path, temporary)
        shutil.rmtree(destination, ignore_errors=True)
        temporary.replace(destination)
    except Exception:
        shutil.rmtree(temporary, ignore_errors=True)
        raise


def _find_server(runtime_dir: Path, server_name: str) -> Path:
    matches = sorted(runtime_dir.rglob(server_name))
    if len(matches) != 1:
        raise RuntimeError(
            f"Expected exactly one {server_name} in llama.cpp runtime, found {len(matches)}"
        )
    server = matches[0]
    if os.name != "nt":
        server.chmod(server.stat().st_mode | stat.S_IXUSR)
    return server


def _paths(cache_dir: Path | None = None) -> dict[str, Path]:
    root = (cache_dir or default_cache_dir()).expanduser().resolve()
    asset = runtime_asset()
    platform_key = f"{asset.system.casefold()}-{asset.machine}"
    runtime_dir = root / "runtime" / DEFAULT_LLAMA_BUILD / platform_key
    return {
        "root": root,
        "runtime_dir": runtime_dir,
        "server": runtime_dir / asset.server_name,
        "model": root / "models" / DEFAULT_MODEL_FILE,
        "manifest": root / "local-llm-manifest.json",
        "downloads": root / "downloads",
    }


def _resolve_server(runtime_dir: Path, asset: RuntimeAsset) -> Path | None:
    direct = runtime_dir / asset.server_name
    if direct.is_file():
        return direct
    matches = sorted(runtime_dir.rglob(asset.server_name)) if runtime_dir.exists() else []
    return matches[0] if len(matches) == 1 else None


def local_llm_status(cache_dir: Path | None = None) -> dict[str, Any]:
    """Return the managed local runtime/model installation state."""
    asset = runtime_asset()
    paths = _paths(cache_dir)
    server = _resolve_server(paths["runtime_dir"], asset)
    model = paths["model"]
    manifest = paths["manifest"]
    return {
        "cache_dir": str(paths["root"]),
        "runtime": {
            "build": DEFAULT_LLAMA_BUILD,
            "asset": asset.archive_name,
            "sha256": asset.sha256,
            "server": str(server) if server else None,
            "ready": bool(server and server.is_file()),
        },
        "model": {
            "repository": DEFAULT_MODEL_REPOSITORY,
            "revision": DEFAULT_MODEL_REVISION,
            "filename": DEFAULT_MODEL_FILE,
            "sha256": DEFAULT_MODEL_SHA256,
            "size_bytes": DEFAULT_MODEL_SIZE,
            "license": DEFAULT_MODEL_LICENSE,
            "path": str(model),
            "ready": model.is_file() and model.stat().st_size == DEFAULT_MODEL_SIZE,
        },
        "manifest": str(manifest) if manifest.is_file() else None,
    }


def setup_local_llm(
    cache_dir: Path | None = None,
    *,
    force: bool = False,
    show_progress: bool = False,
) -> dict[str, Any]:
    """Install the pinned CPU llama.cpp runtime and default Qwen GGUF model."""
    asset = runtime_asset()
    paths = _paths(cache_dir)
    paths["root"].mkdir(parents=True, exist_ok=True)
    downloads = paths["downloads"]
    downloads.mkdir(parents=True, exist_ok=True)

    server = _resolve_server(paths["runtime_dir"], asset)
    if force or server is None:
        archive_path = downloads / asset.archive_name
        _download_verified(asset.url, archive_path, asset.sha256)
        _extract_runtime(archive_path, paths["runtime_dir"])
        archive_path.unlink(missing_ok=True)
        server = _find_server(paths["runtime_dir"], asset.server_name)

    assert server is not None
    _check_server_executable(server)

    model = paths["model"]
    model_ready = model.is_file() and model.stat().st_size == DEFAULT_MODEL_SIZE
    if model_ready and not force:
        model_ready = _sha256(model) == DEFAULT_MODEL_SHA256
    if force or not model_ready:
        _download_verified_ranges(
            DEFAULT_MODEL_URL,
            model,
            DEFAULT_MODEL_SHA256,
            expected_size=DEFAULT_MODEL_SIZE,
            progress=_model_progress if show_progress else None,
        )

    manifest = {
        "schema_version": "prideqc-local-llm-install-v1",
        "runtime": {
            "project": "ggml-org/llama.cpp",
            "build": DEFAULT_LLAMA_BUILD,
            "asset": asset.archive_name,
            "asset_sha256": asset.sha256,
            "server": str(server),
        },
        "model": {
            "repository": DEFAULT_MODEL_REPOSITORY,
            "revision": DEFAULT_MODEL_REVISION,
            "filename": DEFAULT_MODEL_FILE,
            "sha256": DEFAULT_MODEL_SHA256,
            "size_bytes": DEFAULT_MODEL_SIZE,
            "license": DEFAULT_MODEL_LICENSE,
            "path": str(model),
        },
        "compiler_required": False,
    }
    paths["manifest"].write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    status = local_llm_status(paths["root"])
    if not status["runtime"]["ready"] or not status["model"]["ready"]:
        raise RuntimeError("Local LLM setup completed without a usable runtime/model")
    return status


def managed_server_path(cache_dir: Path | None = None) -> Path:
    """Return the installed managed llama-server path or fail with setup guidance."""
    asset = runtime_asset()
    paths = _paths(cache_dir)
    server = _resolve_server(paths["runtime_dir"], asset)
    if server is None:
        raise RuntimeError("Local llama.cpp runtime is not installed; run `prideqc llm setup`.")
    return server


def managed_model_path(cache_dir: Path | None = None) -> Path:
    """Return the installed default GGUF path or fail with setup guidance."""
    path = _paths(cache_dir)["model"]
    if not path.is_file() or path.stat().st_size != DEFAULT_MODEL_SIZE:
        raise RuntimeError("Local Qwen model is not installed; run `prideqc llm setup`.")
    return path
