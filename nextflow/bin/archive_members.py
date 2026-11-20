#!/usr/bin/env python3
"""Run members of a PRIDE ZIP archive: classify, list remotely, extract one at a time.

A PRIDE file such as `runs.zip` or `run01.raw.zip` can hold one or many MS runs. This module
reads only the ZIP central directory to list the members (locally, or remotely with HTTP
range reads of the last bytes of the file), and extracts one run member at a time to
node-local scratch, so the scratch use is the size of the largest run, not of the archive.

The same rules decide which members are runs, both when the archive is inspected before a
batch and when it is unpacked on the compute node:
  - Thermo .raw files, .mzML and .mzML.gz files are run members (QC inputs);
  - a .raw.gz member, a nested archive, a Bruker .d or Waters .raw folder, a .wiff file and other vendor
    files are "unsupported" (the archive is still processed for its supported runs, but the
    project cannot be complete);
  - results, peak lists, text and other files are ignored;
  - an mzML member is not a run when a .raw member of the same run stem exists (vendor first).
Python 3.9 compatible (Codon host python3). No network library beyond urllib.
"""
from __future__ import annotations

import hashlib
import http.client
import io
import re
import stat
import time
import urllib.error
import urllib.request
import zipfile
from collections.abc import Callable, Iterator
from pathlib import Path, PurePosixPath
from typing import IO

RUN_SUFFIXES = (".mzml.gz", ".mzml", ".raw")
UNSUPPORTED_SUFFIXES = (".raw.gz", ".zip", ".tar.gz", ".tgz", ".tar", ".7z", ".rar", ".wiff", ".wiff2",
                        ".wiff.scan", ".baf", ".tdf", ".tdf_bin", ".yep", ".lcd", ".mzxml", ".mzxml.gz")
FOLDER_RUN_SUFFIXES = (".d", ".raw")
MAX_MEMBERS = 1_000_000
MAX_MEMBER_BYTES = 200 * 1024**3
SAFE_BASENAME = re.compile(r"[A-Za-z0-9_.+\-]{1,255}\Z")


class ArchiveMemberError(ValueError):
    """Archive is unsafe, corrupt or ambiguous; the caller writes a fail-closed receipt."""


def run_stem(name: str) -> str:
    lower = PurePosixPath(name).name.casefold()
    for suffix in RUN_SUFFIXES:
        if lower.endswith(suffix):
            return lower[: -len(suffix)]
    return lower


def member_kind(name: str) -> str:
    """'run', 'unsupported' or 'ignored' for one ZIP member name.

    A folder run (Bruker `x.d/`, Waters `x.raw/`) is a run without a reader: 'unsupported'.
    """
    path = PurePosixPath(name)
    lower = name.casefold().rstrip("/")
    if name.endswith("/"):
        return "unsupported" if lower.endswith(FOLDER_RUN_SUFFIXES) else "ignored"
    if folder_run(path) is not None:
        return "unsupported"
    if path.name.startswith(".") or "__macosx" in lower:
        return "ignored"
    if lower.endswith(UNSUPPORTED_SUFFIXES):
        return "unsupported"
    if lower.endswith(RUN_SUFFIXES):
        return "run"
    return "ignored"


def folder_run(path: PurePosixPath) -> str | None:
    """The first parent folder that is a vendor run folder (x.d, x.raw), or None."""
    for i, part in enumerate(path.parts[:-1]):
        if part.casefold().endswith(FOLDER_RUN_SUFFIXES):
            return "/".join(path.parts[: i + 1]) + "/"
    return None


def _safe_member(name: str) -> PurePosixPath:
    if not name or "\\" in name or "\x00" in name or re.search(r"[\x00-\x1f]", name):
        raise ArchiveMemberError(f"Unsafe archive member name: {name!r}")
    path = PurePosixPath(name)
    if path.is_absolute() or any(part in {".", ".."} for part in name.rstrip("/").split("/")):
        raise ArchiveMemberError(f"Unsafe archive member path: {name!r}")
    if any(part.endswith(":") for part in path.parts):
        raise ArchiveMemberError("Windows drive path in archive")
    return path


def safe_basename(name: str) -> str:
    """The local file name of a run member: other characters become '_'.

    PXD011935 has 'hog1_exp1a_5min_GPF(1,2).RAW'. The member path in the receipt keeps the
    original name (the SDRF link); only the scratch file and the mzQC name use this form.
    """
    return re.sub(r"[^A-Za-z0-9_.+\-]", "_", name)


def plan_members(infos: list[zipfile.ZipInfo]) -> dict:
    """Classify the members of one ZIP. Raises ArchiveMemberError on an unsafe or ambiguous archive.

    Returns {"runs": [{member, basename, size, format}], "unsupported": [names], "ignored": n}.
    """
    if len(infos) > MAX_MEMBERS:
        raise ArchiveMemberError("ZIP contains too many members")
    names: set[str] = set()
    raws: dict[str, dict] = {}
    mzmls: dict[str, dict] = {}
    unsupported: list[str] = []
    ignored = 0
    unsafe_ignored = 0
    for info in infos:
        kind = member_kind(info.filename)
        try:
            path = _safe_member(info.filename)
        except ArchiveMemberError:
            # Only run members are extracted (by basename). An unsafe name of another member
            # (PXD010595: '../mqpar.xml' next to the runs) is not a risk: ignore it.
            if kind != "ignored":
                raise
            unsafe_ignored += 1
            continue
        key = path.as_posix().rstrip("/").casefold()
        if key in names:
            raise ArchiveMemberError(f"Duplicate case-insensitive ZIP member: {info.filename}")
        names.add(key)
        if stat.S_ISLNK(info.external_attr >> 16):
            raise ArchiveMemberError("ZIP symlink not allowed")
        kind = member_kind(info.filename)
        if kind == "ignored":
            ignored += 1
            continue
        if kind == "unsupported":
            # one entry for a vendor run folder, not one for each file inside it
            label = folder_run(path) or info.filename
            if info.filename.endswith("/"):
                label = info.filename
            if label not in unsupported:
                unsupported.append(label)
            continue
        if info.file_size > MAX_MEMBER_BYTES or (info.compress_size == 0 and info.file_size > 0):
            raise ArchiveMemberError(f"Invalid ZIP member size: {info.filename}")
        basename = safe_basename(path.name)
        if not SAFE_BASENAME.fullmatch(basename):
            raise ArchiveMemberError(f"Unsafe run file name in archive: {info.filename!r}")
        lower = basename.casefold()
        record = {"member": info.filename, "basename": basename, "size": info.file_size,
                  "format": "raw" if lower.endswith(".raw") else ("mzml_gz" if lower.endswith(".gz") else "mzml")}
        target = raws if record["format"] == "raw" else mzmls
        stem = run_stem(basename)
        if stem in target:
            raise ArchiveMemberError(f"Duplicate run in archive: {target[stem]['member']} and {info.filename}")
        target[stem] = record
    runs = list(raws.values()) + [r for stem, r in mzmls.items() if stem not in raws]
    seen: set[str] = set()
    for r in runs:
        if r["basename"].casefold() in seen:
            raise ArchiveMemberError(f"Duplicate run file name in archive: {r['basename']}")
        seen.add(r["basename"].casefold())
    runs.sort(key=lambda r: r["member"].casefold())
    return {"runs": runs, "unsupported": sorted(unsupported, key=str.casefold), "ignored": ignored + unsafe_ignored,
            "unsafe_ignored": unsafe_ignored,
            "skipped_mzml_with_raw": sum(1 for stem in mzmls if stem in raws)}


class HttpRangeFile(io.RawIOBase):
    """Read-only, seekable file over HTTP range requests (enough for zipfile to read the directory).

    The first request asks for the last `tail` bytes (a suffix range). The response gives the
    true file size (Content-Range) and, for most archives, the whole central directory, so
    one request lists the archive.
    """

    def __init__(self, url: str, *, opener: Callable | None = None, timeout: int = 60,
                 tail: int = 1 << 20, block: int = 1 << 20, retries: int = 4, backoff: float = 2.0) -> None:
        super().__init__()
        self.url, self.pos, self.timeout, self.block = url, 0, timeout, block
        self.retries, self.backoff = retries, backoff
        self.opener = opener or urllib.request.urlopen
        self.requests = 0
        data, total = self._fetch(f"-{tail}")
        self.size = total
        self._cache: tuple[int, bytes] = (total - len(data), data)

    def seekable(self) -> bool:
        return True

    def readable(self) -> bool:
        return True

    def tell(self) -> int:
        return self.pos

    def seek(self, offset: int, whence: int = 0) -> int:
        base = {0: 0, 1: self.pos, 2: self.size}[whence]
        self.pos = max(0, base + offset)
        return self.pos

    def _fetch(self, spec: str) -> tuple[bytes, int]:
        """One range request, retried on a network error (trial 2 with 32 threads: 1,202 of 1,321
        archives had a URLError; one of them, tried again alone, listed without error)."""
        for attempt in range(self.retries + 1):
            try:
                return self._fetch_once(spec)
            except (urllib.error.URLError, OSError, http.client.HTTPException) as err:
                code = getattr(err, "code", None)
                # 403 is throttling on ftp.pride.ebi.ac.uk (trial 6: 244 x 403 for PXD010595 under
                # load, the same file gave 206 alone); 404 and other 4xx are final
                if attempt == self.retries or (code is not None and code < 500 and code not in (403, 429)):
                    raise
                time.sleep(self.backoff * (2 ** attempt))
        raise ArchiveMemberError("unreachable")

    def _fetch_once(self, spec: str) -> tuple[bytes, int]:
        self.requests += 1
        request = urllib.request.Request(self.url, headers={"Range": f"bytes={spec}",
                                                           "User-Agent": "prideqc-archive-inspect/1"})
        with self.opener(request, timeout=self.timeout) as response:
            status = getattr(response, "status", None) or response.getcode()
            if status != 206:
                raise ArchiveMemberError(f"Server does not support range reads (HTTP {status})")
            match = re.fullmatch(r"bytes (\d+)-(\d+)/(\d+)", str(response.headers.get("Content-Range", "")).strip())
            if not match:
                raise ArchiveMemberError("Range response without Content-Range")
            data = response.read()
        start, end, total = (int(x) for x in match.groups())
        if len(data) != end - start + 1:
            raise ArchiveMemberError("Short range read")
        return data, total

    def read(self, n: int = -1) -> bytes:
        if self.pos >= self.size:
            return b""
        if n is None or n < 0:
            n = self.size - self.pos
        n = min(n, self.size - self.pos)
        start, end = self.pos, self.pos + n - 1
        cached_start, data = self._cache
        if not (cached_start <= start and end < cached_start + len(data)):
            # read at least one block, so the many small reads of the directory parser are cheap
            fetch_end = min(self.size - 1, max(end, start + self.block - 1))
            data, _ = self._fetch(f"{start}-{fetch_end}")
            self._cache = (start, data)
            cached_start = start
        out = data[start - cached_start:end - cached_start + 1]
        self.pos += len(out)
        return out

    def readinto(self, b) -> int:  # type: ignore[no-untyped-def]
        data = self.read(len(b))
        b[:len(data)] = data
        return len(data)


def list_zip(source: str | Path | IO[bytes], *, opener: Callable | None = None, attempts: int = 2) -> dict:
    """Classify the members of a local ZIP path, a file object, or a remote https URL.

    A remote listing that is not a valid ZIP is read again once (a bad response is not a
    bad archive). Raises ArchiveMemberError.
    """
    remote = isinstance(source, str) and source.startswith(("http://", "https://"))
    for attempt in range(attempts if remote else 1):
        if remote:
            handle: IO[bytes] = io.BufferedReader(HttpRangeFile(str(source), opener=opener), buffer_size=1 << 16)
        elif isinstance(source, (str, Path)):
            handle = open(source, "rb")
        else:
            handle = source  # type: ignore[assignment]
        try:
            with zipfile.ZipFile(handle) as zf:
                return plan_members(zf.infolist())
        except zipfile.BadZipFile as err:
            if attempt + 1 >= (attempts if remote else 1):
                raise ArchiveMemberError(f"Not a valid ZIP archive: {err}") from None
        finally:
            if handle is not source:
                handle.close()
    raise ArchiveMemberError("Not a valid ZIP archive")


def extract_member(archive: Path, member: str, dest_dir: Path, *, max_bytes: int = MAX_MEMBER_BYTES) -> dict:
    """Extract one run member to dest_dir/<basename>; return {path, size, sha256} of the member bytes."""
    with zipfile.ZipFile(archive) as zf:
        info = zf.getinfo(member)
        plan = plan_members([info])
        if len(plan["runs"]) != 1:
            raise ArchiveMemberError(f"Not a run member: {member}")
        dest_dir.mkdir(parents=True, exist_ok=True)
        target = dest_dir / plan["runs"][0]["basename"]
        digest = hashlib.sha256()
        written = 0
        with zf.open(info) as src, target.open("xb") as out:
            while chunk := src.read(1 << 20):
                written += len(chunk)
                if written > max_bytes:
                    raise ArchiveMemberError("Member expansion exceeds the byte limit")
                digest.update(chunk)
                out.write(chunk)
        # zipfile checks the CRC-32 at the end of the stream; a size mismatch is corruption
        if written != info.file_size:
            raise ArchiveMemberError(f"Member size mismatch: {member}")
    return {"path": target, "size": written, "sha256": digest.hexdigest()}


def iter_runs(archive: Path) -> Iterator[dict]:
    """The run members of a local ZIP in a fixed order (member name), as plan_members gives them."""
    yield from list_zip(archive)["runs"]
