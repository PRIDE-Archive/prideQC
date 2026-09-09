"""PRIDE acquisition through the maintained pridepy download client."""

from __future__ import annotations

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


@dataclass(frozen=True, slots=True)
class DownloadOptions:
    protocol: str = "ftp"
    checksum_check: bool = True
    aspera_maximum_bandwidth: str = "100M"

    def __post_init__(self) -> None:
        if self.protocol not in {"ftp", "aspera", "s3", "globus"}:
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
    and a nonempty regular file exists. Credentials are never serialized.
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
                try:
                    with TemporaryDirectory(prefix=".prideqc-", dir=target) as staging:
                        transfer_result = client.download_file_by_name(
                            accession=accession, file_name=name, output_folder=staging,
                            skip_if_downloaded_already=False, protocol=options.protocol,
                            username=username, password=password,
                            aspera_maximum_bandwidth=options.aspera_maximum_bandwidth,
                            checksum_check=options.checksum_check,
                        )
                        if transfer_result is False:
                            raise RuntimeError("pridepy reported an unsuccessful transfer.")
                        downloaded = Path(staging) / name
                        if downloaded.is_symlink() or not downloaded.is_file() or downloaded.stat().st_size == 0:
                            raise RuntimeError("pridepy did not produce a nonempty regular file.")
                        record["size_bytes"] = downloaded.stat().st_size
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
