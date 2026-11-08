"""PRIDE acquisition through the maintained pridepy download client."""

from __future__ import annotations

import hashlib
import re
from collections.abc import Iterable
from dataclasses import asdict, dataclass
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Any, Protocol

from prideqc.io import write_json
from prideqc.validation import dependency_version


class DownloadClient(Protocol):
    """The subset of pridepy's public Client API used for acquisition."""

    def download_file_by_name(
        self, *, accession: str, file_name: str, output_folder: str,
        skip_if_downloaded_already: bool, protocol: str, username: str | None,
        password: str | None, aspera_maximum_bandwidth: str, checksum_check: bool,
    ) -> Any: ...


def _expected_size_bytes(client: Any, accession: str, name: str) -> int | None:
    """Return a unique positive PRIDE API file size when available.

    pridepy validates only existence/non-emptiness when checksum checking is
    disabled. A timed-out transfer can therefore leave a non-empty partial
    file that looks valid to pridepy. Use the normal project file listing as
    an independent size guard without depending on the checksum endpoint.
    """
    lookup = getattr(client, "get_file_from_api", None)
    if not callable(lookup):
        return None
    try:
        records = lookup(accession, name)
    except Exception:
        return None
    sizes: set[int] = set()
    for record in records or []:
        if not isinstance(record, dict) or record.get("fileName") != name:
            continue
        value = record.get("fileSizeBytes")
        if value is None:
            continue
        try:
            size = int(value)
        except (TypeError, ValueError):
            continue
        if size > 0:
            sizes.add(size)
    if len(sizes) == 1:
        return sizes.pop()
    return None


def _listing_records(client: Any, accession: str, name: str) -> list[dict]:
    lookup = getattr(client, "get_file_from_api", None)
    if not callable(lookup):
        return []
    try:
        records = lookup(accession, name)
    except Exception:
        return []
    return [r for r in records or [] if isinstance(r, dict) and r.get("fileName") == name]


def _expected_sha1(client: Any, accession: str, name: str) -> str | None:
    """Return the unique SHA1 from the PRIDE project file listing, if any.

    2026-11-08: /ws/archive/v3/files/checksum/<PXD> returned HTTP 500 for 38 of 40
    sampled projects, and pridepy aborts the download when that call fails. The file
    listing carries a per-file SHA1 (verified equal to sha1sum on complete downloads),
    so we verify that instead and do not depend on the checksum endpoint.
    """
    values = {str(r.get("checksum", "")).strip().lower() for r in _listing_records(client, accession, name)}
    values = {v for v in values if re.fullmatch(r"[0-9a-f]{40}", v)}
    return values.pop() if len(values) == 1 else None


def _sha1(path: Path) -> str:
    digest = hashlib.sha1()  # noqa: S324 - integrity check against PRIDE's published SHA1
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _gzip_isize(path: Path) -> int | None:
    """Uncompressed size modulo 2**32 from the gzip trailer, or None if not gzip."""
    try:
        with path.open("rb") as handle:
            if handle.read(2) != b"\x1f\x8b":
                return None
            handle.seek(-4, 2)
            return int.from_bytes(handle.read(4), "little")
    except OSError:
        return None


def _size_matches(path: Path, actual: int, expected: int) -> str | None:
    """Return the basis on which the downloaded size matches the PRIDE listing.

    PRIDE lists some gzip files (e.g. PXD006392 *.mzML.gz, 2026-11-08) with the
    UNCOMPRESSED size: the server sent 183,645,619 bytes, the listing says
    460,186,169, and the gzip trailer ISIZE is 460,186,169. Accept that case only
    when the trailer agrees. Gzip CRC and full length are verified later by
    archive_units when the stream is read to EOF.
    """
    if actual == expected:
        return "bytes"
    if path.name.casefold().endswith(".gz") and _gzip_isize(path) == expected % (1 << 32):
        return "gzip_isize"
    return None


@dataclass(frozen=True, slots=True)
class DownloadOptions:
    protocol: str = "ftp"
    checksum_check: bool = True
    aspera_maximum_bandwidth: str = "100M"

    def __post_init__(self) -> None:
        if self.protocol not in {"ftp", "aspera", "s3", "globus", "fire"}:
            raise ValueError("Unsupported pridepy transfer protocol.")
        if not self.aspera_maximum_bandwidth.strip():
            raise ValueError("Aspera bandwidth must be nonempty.")


def selected_files(names: Iterable[str]) -> list[str]:
    """Deduplicate repeated SDRF rows, rejecting unsafe or colliding flat names."""
    result: list[str] = []
    folded: dict[str, str] = {}
    for name in names:
        if (not name or name != name.strip() or name in {".", ".."}
                or any(character in name for character in "/\\:\x00\r\n\t")
                or name.casefold() in {"not available", "not provided", "not applicable"}):
            raise ValueError(f"PRIDE selections must be exact file basenames: {name!r}")
        previous = folded.get(name.casefold())
        if previous is not None and previous != name:
            raise ValueError(
                f"Case-insensitive download filename collision: {previous!r}, {name!r}",
            )
        if previous is None:
            folded[name.casefold()] = name
            result.append(name)
    if not result:
        raise ValueError("Select at least one PRIDE file with --file or --sdrf.")
    if "download-manifest.json" in folded:
        raise ValueError("download-manifest.json is reserved for acquisition provenance.")
    return result


class PrideRepository:
    """Acquire an explicit file subset; pridepy owns transport and checksum logic.

    A fresh destination prevents reuse of partial or unrelated files. Each file
    is staged on the same filesystem, then published only after pridepy returns
    and a nonempty regular file exists. When the PRIDE listing exposes a file
    size, that size must also match before publication. Credentials are never
    serialized.
    """

    def __init__(self, client: DownloadClient | None = None) -> None:
        self._client = client

    def _get_client(self) -> DownloadClient:
        if self._client is None:
            try:
                from pridepy.download.client import Client
            except ImportError as exc:
                raise RuntimeError("PRIDE acquisition requires pridepy>=0.0.16; run `uv sync`.") from exc
            self._client = Client()
        assert self._client is not None
        return self._client

    def download_files(
        self, accession: str, filenames: Iterable[str], destination: Path | str, *,
        options: DownloadOptions | None = None, username: str | None = None,
        password: str | None = None,
    ) -> list[Path]:
        accession = accession.strip().upper()
        if not re.fullmatch(r"PXD\d{6,}", accession):
            raise ValueError("Expected a PRIDE project accession such as PXD008644.")
        names = selected_files(filenames)
        options = options or DownloadOptions()
        if bool(username) != bool(password):
            raise ValueError("Private PRIDE access requires both username and password.")
        target = Path(destination).resolve()
        if target.exists() and any(target.iterdir()):
            raise FileExistsError(f"Download directory is not empty: {target}")
        client = self._get_client()
        target.mkdir(parents=True, exist_ok=True)
        records: list[dict[str, Any]] = []
        paths: list[Path] = []
        try:
            for name in names:
                record: dict[str, Any] = {"filename": name, "status": "failed"}
                records.append(record)
                expected_size = _expected_size_bytes(client, accession, name)
                if expected_size is not None:
                    record["expected_size_bytes"] = expected_size
                # With a listing SHA1, verify it ourselves and skip pridepy's checksum-endpoint path.
                expected_sha1 = _expected_sha1(client, accession, name) if options.checksum_check else None
                pridepy_checksum = options.checksum_check and expected_sha1 is None
                record["checksum_method"] = ("listing_sha1" if expected_sha1 else
                                             "pridepy" if pridepy_checksum else "none")
                try:
                    with TemporaryDirectory(prefix=".prideqc-", dir=target) as staging:
                        transfer_result = client.download_file_by_name(
                            accession=accession, file_name=name, output_folder=staging,
                            skip_if_downloaded_already=False, protocol=options.protocol,
                            username=username, password=password,
                            aspera_maximum_bandwidth=options.aspera_maximum_bandwidth,
                            checksum_check=pridepy_checksum,
                        )
                        if transfer_result is False:
                            raise RuntimeError("pridepy reported an unsuccessful transfer.")
                        downloaded = Path(staging) / name
                        if downloaded.is_symlink() or not downloaded.is_file() or downloaded.stat().st_size == 0:
                            raise RuntimeError("pridepy did not produce a nonempty regular file.")
                        actual_size = downloaded.stat().st_size
                        record["size_bytes"] = actual_size
                        basis = (_size_matches(downloaded, actual_size, expected_size)
                                 if expected_size is not None else None)
                        if expected_size is not None and basis is None:
                            record["size_match"] = False
                            raise RuntimeError(
                                "pridepy produced an incomplete download: "
                                f"got {actual_size} bytes, expected {expected_size}."
                            )
                        if expected_size is not None:
                            record["size_match"] = True
                            record["size_match_basis"] = basis
                        if expected_sha1 is not None:
                            record["sha1_match"] = _sha1(downloaded) == expected_sha1
                            if not record["sha1_match"]:
                                raise RuntimeError("Downloaded file SHA1 does not match the PRIDE listing.")
                        final = target / name
                        downloaded.replace(final)
                        paths.append(final)
                        record["status"] = "downloaded"
                except Exception as exc:
                    # Upstream errors may contain authenticated URLs. Do not copy
                    # their text into manifests or the public CLI error message.
                    record["error_type"] = type(exc).__name__
                    raise RuntimeError(
                        f"PRIDE download failed for {accession}/{name} ({type(exc).__name__}). "
                        f"See {target / 'download-manifest.json'}."
                    ) from None
        finally:
            write_json(target / "download-manifest.json", {
                "engine": "pridepy", "engine_version": dependency_version("pridepy"),
                "accession": accession, "options": asdict(options), "files": records,
                "unprocessed": names[len(records):], "success": len(paths) == len(names),
            })
        return paths
