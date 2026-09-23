from __future__ import annotations

import hashlib
import io
import json
import tarfile
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from prideqc import local_llm
from prideqc.cli import parser


class LocalLLMTests(unittest.TestCase):
    def test_runtime_assets_cover_reference_laptop_platforms(self) -> None:
        linux = local_llm.runtime_asset(system="Linux", machine="AMD64")
        mac = local_llm.runtime_asset(system="Darwin", machine="arm64")
        windows = local_llm.runtime_asset(system="Windows", machine="x86_64")
        self.assertEqual(linux.archive_name, "llama-b11046-bin-ubuntu-x64.tar.gz")
        self.assertEqual(mac.server_name, "llama-server")
        self.assertEqual(windows.server_name, "llama-server.exe")
        with self.assertRaisesRegex(RuntimeError, "No managed llama.cpp"):
            local_llm.runtime_asset(system="Plan9", machine="mips")

    def test_llm_cli_has_setup_status_and_adjudicate(self) -> None:
        setup = parser().parse_args(["llm", "setup", "--force"])
        status = parser().parse_args(["llm", "status"])
        adjudicate = parser().parse_args(
            ["llm", "adjudicate", "--request", "request.json"]
        )
        self.assertEqual(setup.llm_command, "setup")
        self.assertTrue(setup.force)
        self.assertEqual(status.llm_command, "status")
        self.assertEqual(adjudicate.llm_command, "adjudicate")

    def test_verified_download_is_atomic_and_rejects_wrong_hash(self) -> None:
        data = b"verified payload"
        digest = hashlib.sha256(data).hexdigest()
        with tempfile.TemporaryDirectory() as folder:
            destination = Path(folder) / "artifact.bin"

            def stream(_url: str, output: io.BufferedWriter) -> int:
                output.write(data)
                return len(data)

            with patch("prideqc.local_llm._stream_download", side_effect=stream):
                result = local_llm._download_verified(
                    "https://example.invalid/artifact",
                    destination,
                    digest,
                    expected_size=len(data),
                )
            self.assertEqual(result.read_bytes(), data)
            self.assertFalse(destination.with_name("artifact.bin.part").exists())

            destination.unlink()
            with patch("prideqc.local_llm._stream_download", side_effect=stream):
                with self.assertRaisesRegex(RuntimeError, "SHA-256 mismatch"):
                    local_llm._download_verified(
                        "https://example.invalid/artifact",
                        destination,
                        "0" * 64,
                    )
            self.assertFalse(destination.exists())
            self.assertFalse(destination.with_name("artifact.bin.part").exists())


    def test_stream_download_follows_http_partial_content_ranges(self) -> None:
        payload = b"abcdefghijklmnopqrstuvwxyz"

        class FakeHeaders(dict[str, str]):
            pass

        class FakeResponse(io.BytesIO):
            def __init__(self, data: bytes, *, content_range: str) -> None:
                super().__init__(data)
                self.status = 206
                self.headers = FakeHeaders({"Content-Range": content_range})

            def __enter__(self) -> FakeResponse:
                return self

            def __exit__(self, *_args: object) -> None:
                self.close()

        requests: list[str | None] = []

        def urlopen(request: object, timeout: int) -> FakeResponse:
            self.assertEqual(timeout, 120)
            range_header = request.get_header("Range")  # type: ignore[attr-defined]
            requests.append(range_header)
            if range_header is None:
                return FakeResponse(
                    payload[:10],
                    content_range=f"bytes 0-9/{len(payload)}",
                )
            if range_header == "bytes=10-":
                return FakeResponse(
                    payload[10:20],
                    content_range=f"bytes 10-19/{len(payload)}",
                )
            if range_header == "bytes=20-":
                return FakeResponse(
                    payload[20:],
                    content_range=f"bytes 20-25/{len(payload)}",
                )
            raise AssertionError(range_header)

        output = io.BytesIO()
        with patch("prideqc.local_llm.urllib.request.urlopen", side_effect=urlopen):
            size = local_llm._stream_download("https://example.invalid/model", output)

        self.assertEqual(size, len(payload))
        self.assertEqual(output.getvalue(), payload)
        self.assertEqual(requests, [None, "bytes=10-", "bytes=20-"])

    def test_tar_extraction_rejects_path_traversal(self) -> None:
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            archive_path = root / "bad.tar.gz"
            with tarfile.open(archive_path, "w:gz") as archive:
                payload = b"escape"
                info = tarfile.TarInfo("../escape.txt")
                info.size = len(payload)
                archive.addfile(info, io.BytesIO(payload))
            with self.assertRaisesRegex(RuntimeError, "Unsafe path"):
                local_llm._safe_extract_tar(archive_path, root / "runtime")
            self.assertFalse((root / "escape.txt").exists())

    def test_setup_uses_pinned_downloads_without_compiler(self) -> None:
        runtime_bytes = io.BytesIO()
        with tarfile.open(fileobj=runtime_bytes, mode="w:gz") as archive:
            payload = b"#!/bin/sh\nexit 0\n"
            info = tarfile.TarInfo("llama-server")
            info.mode = 0o755
            info.size = len(payload)
            archive.addfile(info, io.BytesIO(payload))
        runtime_data = runtime_bytes.getvalue()
        model_data = b"tiny-test-model"
        asset = local_llm.RuntimeAsset(
            system="Linux",
            machine="x86_64",
            archive_name="llama-test.tar.gz",
            sha256=hashlib.sha256(runtime_data).hexdigest(),
            server_name="llama-server",
        )

        def stream(url: str, output: io.BufferedWriter) -> int:
            data = runtime_data if "github" in url else model_data
            output.write(data)
            return len(data)

        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            with (
                patch("prideqc.local_llm.runtime_asset", return_value=asset),
                patch("prideqc.local_llm._stream_download", side_effect=stream),
                patch.object(local_llm, "DEFAULT_MODEL_SIZE", len(model_data)),
                patch.object(
                    local_llm,
                    "DEFAULT_MODEL_SHA256",
                    hashlib.sha256(model_data).hexdigest(),
                ),
            ):
                status = local_llm.setup_local_llm(root)
            self.assertTrue(status["runtime"]["ready"])
            self.assertTrue(status["model"]["ready"])
            manifest = json.loads((root / "local-llm-manifest.json").read_text())
            self.assertFalse(manifest["compiler_required"])
            self.assertEqual(manifest["model"]["license"], "Apache-2.0")


if __name__ == "__main__":
    unittest.main()
