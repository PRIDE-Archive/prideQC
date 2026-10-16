#!/usr/bin/env python3
"""Bounded, traversal-safe materialization of logical PRIDE MS runs.

One input archive may contain many RAW/mzML runs or multiple complete .d folders.
This module is intentionally independent of pridepy and pyOpenMS for offline tests.
"""
from __future__ import annotations

import gzip
import hashlib
import json
import re
import shutil
import stat
import tarfile
import zipfile
from pathlib import Path, PurePosixPath

RAW_SUFFIXES = (".raw", ".mzml", ".mzml.gz")
ARCHIVE_SUFFIXES = (".zip", ".tar.gz", ".tgz", ".raw.gz")
MAX_FILES = 100_000
MAX_BYTES = 500 * 1024**3
MAX_DEPTH = 2


class ArchiveError(ValueError):
    """Archive is unsafe, corrupt, ambiguous, or contains no supported MS runs."""


def _safe_name(raw: str) -> PurePosixPath:
    if not raw or "\\" in raw or "\x00" in raw or re.search(r"[\x00-\x1f]", raw):
        raise ArchiveError(f"Unsafe archive member name: {raw!r}")
    path = PurePosixPath(raw)
    if path.is_absolute() or any(part in {".", ".."} for part in raw.rstrip("/").split("/")):
        raise ArchiveError(f"Unsafe archive member path: {raw!r}")
    if any(part.endswith(":") for part in path.parts):
        raise ArchiveError("Windows drive path in archive")
    return path


def _safe_copy(src, dst: Path, budget: list[int], limit: int) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists():
        raise ArchiveError(f"Duplicate member path: {dst.name}")
    with dst.open("xb") as handle:
        while chunk := src.read(1024 * 1024):
            budget[0] += len(chunk)
            if budget[0] > limit:
                raise ArchiveError("Archive expansion exceeds maximum byte budget")
            handle.write(chunk)


def _extract_archive(archive: Path, dest: Path, budget: list[int], limit: int) -> None:
    """Extract only regular files/dirs, rejecting links and unsafe paths."""
    dest.mkdir(parents=True, exist_ok=False)
    if archive.name.casefold().endswith(".zip"):
        with zipfile.ZipFile(archive) as zf:
            entries = zf.infolist()
            if len(entries) > MAX_FILES:
                raise ArchiveError("ZIP contains too many members")
            names: set[str] = set()
            for item in entries:
                rel = _safe_name(item.filename)
                key = rel.as_posix().rstrip("/").casefold()
                if key in names:
                    raise ArchiveError(f"Duplicate case-insensitive ZIP member: {item.filename}")
                names.add(key)
                mode = item.external_attr >> 16
                if stat.S_ISLNK(mode):
                    raise ArchiveError("ZIP symlink not allowed")
                if item.is_dir():
                    (dest / rel).mkdir(parents=True, exist_ok=True)
                else:
                    if item.file_size > limit or (item.compress_size == 0 and item.file_size > 0):
                        raise ArchiveError("Invalid ZIP member size")
                    with zf.open(item) as stream:
                        _safe_copy(stream, dest / rel, budget, limit)
    elif archive.name.casefold().endswith((".tar.gz", ".tgz")):
        with tarfile.open(archive, "r:gz") as tf:
            entries = tf.getmembers()
            if len(entries) > MAX_FILES:
                raise ArchiveError("TAR contains too many members")
            names: set[str] = set()
            for item in entries:
                rel = _safe_name(item.name)
                key = rel.as_posix().rstrip("/").casefold()
                if key in names:
                    raise ArchiveError(f"Duplicate case-insensitive TAR member: {item.name}")
                names.add(key)
                if item.isdir():
                    (dest / rel).mkdir(parents=True, exist_ok=True)
                elif item.isfile():
                    if item.size > limit:
                        raise ArchiveError("TAR member exceeds maximum size")
                    stream = tf.extractfile(item)
                    if stream is None:
                        raise ArchiveError("Missing TAR member stream")
                    with stream:
                        _safe_copy(stream, dest / rel, budget, limit)
                else:
                    raise ArchiveError("TAR links and special files are forbidden")
    else:
        raise ArchiveError(f"Unsupported archive: {archive.name}")


def _verify_gzip(source: Path, budget: list[int], limit: int, target: Path | None = None) -> None:
    """Consume gzip to EOF to verify CRC and expansion bounds, optionally unpack."""
    sink = None
    try:
        if target is not None:
            sink = target.open("xb")
        with gzip.open(source, "rb") as reader:
            while chunk := reader.read(1024 * 1024):
                budget[0] += len(chunk)
                if budget[0] > limit:
                    raise ArchiveError("Gzip expansion exceeds maximum byte budget")
                if sink is not None:
                    sink.write(chunk)
    finally:
        if sink is not None:
            sink.close()


def _logical_inputs(root: Path) -> list[Path]:
    candidates: list[Path] = []
    for path in sorted(root.rglob("*"), key=lambda p: p.as_posix().casefold()):
        if path.is_symlink():
            raise ArchiveError("Symlink encountered after extraction")
        name = path.name.casefold()
        if path.is_dir() and name.endswith(".d"):
            candidates.append(path)
        elif path.is_file() and name.endswith(RAW_SUFFIXES):
            if not any(parent.name.casefold().endswith(".d") for parent in path.parents):
                candidates.append(path)
    return [p for p in candidates if not any(parent in candidates for parent in p.parents)]


def expand(source: Path, outdir: Path, *, max_bytes: int = MAX_BYTES,
           depth: int = MAX_DEPTH) -> dict:
    """Materialize self-contained unit folders, output member-level provenance."""
    if outdir.exists() and any(outdir.iterdir()):
        raise ArchiveError("Refusing nonempty output directory")
    if source.is_symlink():
        raise ArchiveError("Symlink source is forbidden")
    source = source.resolve(strict=True)
    if not (source.is_file() or source.is_dir()):
        raise ArchiveError("Source must be a regular file or directory")
    if max_bytes < 1 or depth < 0:
        raise ArchiveError("Invalid expansion limits")
    outdir.mkdir(parents=True, exist_ok=True)
    staging = outdir / "_stage"
    staging.mkdir()
    budget = [0]
    try:
        if source.is_dir():
            if not source.name.casefold().endswith(".d"):
                raise ArchiveError("Only vendor .d directory inputs are supported")
            if any(p.is_symlink() for p in source.rglob("*")):
                raise ArchiveError("Symlink inside vendor folder is forbidden")
            shutil.copytree(source, staging / source.name)
            root = staging
        elif source.name.casefold().endswith((".zip", ".tar.gz", ".tgz")):
            root = staging / "expanded"
            _extract_archive(source, root, budget, max_bytes)
        elif source.name.casefold().endswith(".raw.gz"):
            root = staging
            _verify_gzip(source, budget, max_bytes, root / source.name[:-3])
        elif source.name.casefold().endswith(RAW_SUFFIXES):
            root = staging
            if source.name.casefold().endswith(".mzml.gz"):
                _verify_gzip(source, budget, max_bytes)
            shutil.copyfile(source, root / source.name)
        else:
            raise ArchiveError(f"Unsupported source extension: {source.name}")
        # Support nested ZIPs, e.g. bundle.zip containing run.raw.zip, to bounded depth.
        for level in range(depth):
            nested = sorted(p for p in root.rglob("*") if p.is_file() and
                            p.name.casefold().endswith((".zip", ".tar.gz", ".tgz")) and
                            not any(parent.name.casefold().endswith(".d") for parent in p.parents))
            if not nested:
                break
            for index, item in enumerate(nested):
                dst = item.parent / ("_extracted_%s_%d" % (level, index))
                _extract_archive(item, dst, budget, max_bytes)
                item.unlink()
        # ZIPs sometimes contain individually gzipped Thermo runs.
        if source.name.casefold().endswith((".zip", ".tar.gz", ".tgz")):
            for candidate in sorted(root.rglob("*")):
                if candidate.is_file() and candidate.name.casefold().endswith(".raw.gz"):
                    target = candidate.with_name(candidate.name[:-3])
                    _verify_gzip(candidate, budget, max_bytes, target)
                    candidate.unlink()
        # Validate gzip streams even when they were contained in a ZIP.
        if source.name.casefold().endswith((".zip", ".tar.gz", ".tgz")):
            for candidate in root.rglob("*"):
                if candidate.is_file() and candidate.name.casefold().endswith(".mzml.gz"):
                    _verify_gzip(candidate, budget, max_bytes)
        if any(p.is_file() and p.name.casefold().endswith(ARCHIVE_SUFFIXES)
               for p in root.rglob("*") if not any(t.name.casefold().endswith(".d") for t in p.parents)):
            raise ArchiveError("Nested compression depth exceeded or unsupported .raw.gz")
        items = _logical_inputs(root)
        if not items:
            raise ArchiveError("No supported .raw/.mzML/.mzML.gz/.d inputs found")
        folds: set[str] = set()
        output: list[dict] = []
        for item in items:
            relative = item.relative_to(root).as_posix()
            fingerprint = hashlib.sha256(relative.encode()).hexdigest()[:12]
            unit_id = f"unit-{len(output):05d}-{fingerprint}"
            unit_dir = outdir / unit_id
            unit_dir.mkdir()
            target = unit_dir / item.name
            if item.is_dir():
                shutil.move(str(item), str(target))
                kind = "vendor_d"
            else:
                shutil.move(str(item), str(target))
                kind = "mzml_gz" if item.name.casefold().endswith(".mzml.gz") else (
                    "mzml" if item.name.casefold().endswith(".mzml") else "raw")
            if target.name.casefold() in folds:
                raise ArchiveError(f"Duplicate logical RAW basename: {target.name}")
            folds.add(target.name.casefold())
            record = {"unit_id": unit_id, "source_member": relative,
                      "input_basename": target.name, "format": kind}
            (unit_dir / "unit.json").write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
            output.append(record)
        inventory = {"schema_version": 1, "source_name": source.name,
                     "logical_run_count": len(output), "units": output}
        (outdir / "inventory.json").write_text(json.dumps(inventory, indent=2, sort_keys=True) + "\n")
        return inventory
    finally:
        shutil.rmtree(staging, ignore_errors=True)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--max-uncompressed-bytes", type=int, default=MAX_BYTES)
    arguments = parser.parse_args()
    print(json.dumps(expand(arguments.source, arguments.outdir,
                            max_bytes=arguments.max_uncompressed_bytes), sort_keys=True))
