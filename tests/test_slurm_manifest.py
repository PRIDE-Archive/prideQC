from __future__ import annotations

import importlib.util
import json
import unittest
from pathlib import Path
from unittest import mock

SCRIPT = (
    Path(__file__).parents[1] / "scripts" / "slurm" / "prepare_prideqc_file_manifest.py"
)
SPEC = importlib.util.spec_from_file_location("prepare_prideqc_file_manifest", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
manifest = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(manifest)


class PrideArchiveAliasTests(unittest.TestCase):
    def test_py_alias_resolves_to_unique_date_prefixed_archive_file(self) -> None:
        live = {
            "20111225_EXQ5_KiSh_SA_LabelFree_HeLa_Proteome_Control_rep4_pH11.raw",
            "other.raw",
        }
        resolved, reason = manifest.resolve_pride_archive_name(
            "pY_EXQ5_KiSh_SA_LabelFree_HeLa_Proteome_Control_rep4_pH11.raw",
            live,
        )
        self.assertEqual(
            resolved,
            "20111225_EXQ5_KiSh_SA_LabelFree_HeLa_Proteome_Control_rep4_pH11.raw",
        )
        self.assertEqual(reason, "py-date-alias")

    def test_resolved_alias_is_deduplicated_against_exact_physical_file(self) -> None:
        physical = (
            "20111225_EXQ5_KiSh_SA_LabelFree_HeLa_Proteome_Control_rep4_pH11.raw"
        )
        alias = "pY_EXQ5_KiSh_SA_LabelFree_HeLa_Proteome_Control_rep4_pH11.raw"
        records = [
            self._record(physical),
            self._record(alias),
        ]
        resolved, aliases = manifest.resolve_and_deduplicate_records(
            records, {physical}
        )
        self.assertEqual(len(resolved), 1)
        self.assertEqual(resolved[0]["sdrf_data_file"], physical)
        self.assertEqual(resolved[0]["archive_file"], physical)
        self.assertEqual(len(aliases), 1)
        self.assertEqual(aliases[0]["sdrf_data_file"], alias)
        self.assertEqual(aliases[0]["archive_file"], physical)
        self.assertEqual(aliases[0]["canonical_sdrf_data_file"], physical)
        self.assertEqual(aliases[0]["resolution"], "py-date-alias")

    def test_exact_physical_row_wins_even_when_alias_appears_first(self) -> None:
        physical = (
            "20111225_EXQ5_KiSh_SA_LabelFree_HeLa_Proteome_Control_rep4_pH11.raw"
        )
        alias = "pY_EXQ5_KiSh_SA_LabelFree_HeLa_Proteome_Control_rep4_pH11.raw"
        resolved, aliases = manifest.resolve_and_deduplicate_records(
            [self._record(alias), self._record(physical)], {physical}
        )
        self.assertEqual(len(resolved), 1)
        self.assertEqual(resolved[0]["sdrf_data_file"], physical)
        self.assertEqual(resolved[0]["archive_file"], physical)
        self.assertEqual(aliases[0]["sdrf_data_file"], alias)
        self.assertEqual(aliases[0]["canonical_sdrf_data_file"], physical)

    def test_pxd000612_style_273_rows_reduce_to_231_physical_plus_42_aliases(self) -> None:
        physical_names = [
            f"20111225_EXQ5_KiSh_alias_core_{index:03d}.raw"
            for index in range(42)
        ]
        physical_names.extend(
            f"20120310_EXQ5_KiSh_unique_{index:03d}.raw"
            for index in range(189)
        )
        aliases = [
            f"pY_EXQ5_KiSh_alias_core_{index:03d}.raw"
            for index in range(42)
        ]
        records = [self._record(name) for name in physical_names + aliases]

        resolved, alias_rows = manifest.resolve_and_deduplicate_records(
            records, set(physical_names)
        )

        self.assertEqual(len(records), 273)
        self.assertEqual(len(resolved), 231)
        self.assertEqual(len(alias_rows), 42)
        self.assertTrue(all(row["resolution"] == "py-date-alias" for row in alias_rows))

    def test_unknown_alias_is_not_fuzzy_matched(self) -> None:
        name = "unexpected_EXQ5_file.raw"
        resolved, reason = manifest.resolve_pride_archive_name(
            name, {"20111225_EXQ5_file.raw"}
        )
        self.assertEqual(resolved, name)
        self.assertEqual(reason, "unresolved")

    def test_ambiguous_date_prefixed_alias_is_rejected(self) -> None:
        alias = "pY_EXQ5_file.raw"
        with self.assertRaisesRegex(ValueError, "Ambiguous PRIDE pY/date alias"):
            manifest.resolve_pride_archive_name(
                alias,
                {"20111225_EXQ5_file.raw", "20111226_EXQ5_file.raw"},
            )

    @staticmethod
    def _record(name: str) -> dict[str, str]:
        return {
            "sdrf_data_file": name,
            "archive_file": name,
            "file_uri": f"https://example.invalid/{name}",
            "needs_file_map": "0",
        }


class PrideApiPagingTests(unittest.TestCase):
    def test_page_file_names_accepts_current_v3_bare_list(self) -> None:
        payload = [
            {"fileName": "a.raw", "fileSizeBytes": 1},
            {"fileName": "b.raw", "fileSizeBytes": 2},
        ]
        names, number, total = manifest._page_file_names(payload, "PXD000612")
        self.assertEqual(names, {"a.raw", "b.raw"})
        self.assertIsNone(number)
        self.assertIsNone(total)

    def test_page_file_names_accepts_legacy_embedded_files(self) -> None:
        payload = {
            "_embedded": {"files": [{"fileName": "a.raw"}, {"fileName": "b.raw"}]},
            "page": {"number": 0, "totalPages": 3},
        }
        names, number, total = manifest._page_file_names(payload, "PXD000612")
        self.assertEqual(names, {"a.raw", "b.raw"})
        self.assertEqual(number, 0)
        self.assertEqual(total, 3)

    def test_page_file_names_rejects_records_without_file_name(self) -> None:
        payload = [{"name": "a.raw", "fileSizeBytes": 1}]
        with self.assertRaisesRegex(RuntimeError, "no fileName values"):
            manifest._page_file_names(payload, "PXD000612")

    def test_project_file_names_uses_total_records_header(self) -> None:
        payloads = [
            ([{"fileName": "a.raw"}], {"total_records": "2"}),
            ([{"fileName": "b.raw"}], {"total_records": "2"}),
        ]
        urls: list[str] = []

        def fake_request(
            url: str, timeout: int = 120, *, accept: str = "application/json"
        ) -> tuple[object, dict[str, str]]:
            self.assertEqual(accept, "application/json")
            urls.append(url)
            return payloads[len(urls) - 1]

        with mock.patch.object(manifest, "request_json_with_headers", fake_request):
            names = manifest.pride_project_file_names("PXD000612")

        self.assertEqual(names, {"a.raw", "b.raw"})
        self.assertEqual(len(urls), 2)
        self.assertIn("/projects/PXD000612/files?pageSize=100&page=0", urls[0])
        self.assertIn("/projects/PXD000612/files?pageSize=100&page=1", urls[1])

    def test_project_file_names_pages_bare_lists_until_empty_without_header(self) -> None:
        payloads = [
            ([{"fileName": "a.raw"}], {}),
            ([{"fileName": "b.raw"}], {}),
            ([], {}),
        ]
        calls = 0

        def fake_request(
            url: str, timeout: int = 120, *, accept: str = "application/json"
        ) -> tuple[object, dict[str, str]]:
            nonlocal calls
            response = payloads[calls]
            calls += 1
            return response

        with mock.patch.object(manifest, "request_json_with_headers", fake_request):
            names = manifest.pride_project_file_names("PXD000612")

        self.assertEqual(names, {"a.raw", "b.raw"})
        self.assertEqual(calls, 3)

    def test_project_file_names_rejects_incomplete_pagination(self) -> None:
        payloads = [
            ([{"fileName": "a.raw"}], {"total_records": "2"}),
            ([], {"total_records": "2"}),
        ]
        calls = 0

        def fake_request(
            url: str, timeout: int = 120, *, accept: str = "application/json"
        ) -> tuple[object, dict[str, str]]:
            nonlocal calls
            response = payloads[calls]
            calls += 1
            return response

        with mock.patch.object(manifest, "request_json_with_headers", fake_request):
            with self.assertRaisesRegex(RuntimeError, "total_records=2"):
                manifest.pride_project_file_names("PXD000612")


class HttpHeaderTests(unittest.TestCase):
    def test_github_discovery_requests_github_json_media_type(self) -> None:
        captured: list[str] = []

        def fake_request(
            url: str, timeout: int = 120, *, accept: str = "*/*"
        ) -> bytes:
            captured.append(accept)
            return json.dumps(
                [
                    {
                        "type": "file",
                        "name": "PXD000612.sdrf.tsv",
                        "download_url": "https://example.invalid/PXD000612.sdrf.tsv",
                    }
                ]
            ).encode()

        with mock.patch.object(manifest, "request_bytes", fake_request):
            urls = manifest.discover_urls(
                "PXD000612", "bigbio/sdrf-annotated-datasets", "main"
            )

        self.assertEqual(captured, ["application/vnd.github+json"])
        self.assertEqual(urls, ["https://example.invalid/PXD000612.sdrf.tsv"])

    def test_pride_json_request_defaults_to_application_json(self) -> None:
        seen: dict[str, str] = {}

        class FakeResponse:
            headers = {"total_records": "1"}

            def __enter__(self) -> FakeResponse:
                return self

            def __exit__(self, *args: object) -> None:
                return None

            def read(self) -> bytes:
                return b'[{"fileName": "a.raw"}]'

        def fake_urlopen(request: object, timeout: int = 120) -> FakeResponse:
            seen["accept"] = request.get_header("Accept")  # type: ignore[attr-defined]
            return FakeResponse()

        with mock.patch.object(manifest.urllib.request, "urlopen", fake_urlopen):
            payload, headers = manifest.request_json_with_headers(
                "https://www.ebi.ac.uk/pride/ws/archive/v3/projects/PXD000612/files"
            )

        self.assertEqual(seen["accept"], "application/json")
        self.assertEqual(payload, [{"fileName": "a.raw"}])
        self.assertEqual(headers["total_records"], "1")


class ManifestFailurePolicyTests(unittest.TestCase):
    def test_pride_api_failure_is_fatal_by_default(self) -> None:
        import tempfile
        import urllib.error

        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            argv = [
                "prepare_prideqc_file_manifest.py",
                "--gt",
                str(root / "gt.tsv"),
                "--sdrf-root",
                str(root / "sdrf"),
                "--manifest",
                str(root / "manifest.tsv"),
                "--meta",
                str(root / "manifest.meta.txt"),
            ]
            with (
                mock.patch.object(manifest.sys, "argv", argv),
                mock.patch.object(
                    manifest,
                    "load_gt",
                    return_value=([{"pxd_accession": "PXD000612"}], ["PXD000612"]),
                ),
                mock.patch.object(
                    manifest,
                    "pride_project_file_names",
                    side_effect=urllib.error.URLError("offline"),
                ),
            ):
                with self.assertRaisesRegex(
                    RuntimeError, "refusing to generate a potentially duplicated manifest"
                ):
                    manifest.main()


if __name__ == "__main__":
    unittest.main()
