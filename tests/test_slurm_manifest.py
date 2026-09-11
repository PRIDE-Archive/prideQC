from __future__ import annotations

import importlib.util
from pathlib import Path
import unittest


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
            {
                "sdrf_data_file": physical,
                "archive_file": physical,
                "file_uri": f"https://example.invalid/{physical}",
                "needs_file_map": "0",
            },
            {
                "sdrf_data_file": alias,
                "archive_file": alias,
                "file_uri": f"https://example.invalid/{alias}",
                "needs_file_map": "0",
            },
        ]
        resolved, aliases = manifest.resolve_and_deduplicate_records(
            records, {physical}
        )
        self.assertEqual(len(resolved), 1)
        self.assertEqual(resolved[0]["archive_file"], physical)
        self.assertEqual(len(aliases), 1)
        self.assertEqual(aliases[0]["sdrf_data_file"], alias)
        self.assertEqual(aliases[0]["archive_file"], physical)
        self.assertEqual(aliases[0]["canonical_sdrf_data_file"], physical)

    def test_unknown_alias_is_not_fuzzy_matched(self) -> None:
        name = "unexpected_EXQ5_file.raw"
        resolved, reason = manifest.resolve_pride_archive_name(
            name, {"20111225_EXQ5_file.raw"}
        )
        self.assertEqual(resolved, name)
        self.assertEqual(reason, "unresolved")


if __name__ == "__main__":
    unittest.main()

class PrideApiPagingTests(unittest.TestCase):
    def test_page_file_names_accepts_embedded_files_and_pagination(self) -> None:
        payload = {
            "_embedded": {"files": [{"fileName": "a.raw"}, {"fileName": "b.raw"}]},
            "page": {"number": 0, "totalPages": 3},
        }
        names, number, total = manifest._page_file_names(payload, "PXD000612")
        self.assertEqual(names, {"a.raw", "b.raw"})
        self.assertEqual(number, 0)
        self.assertEqual(total, 3)

    def test_project_file_names_pages_until_total_pages(self) -> None:
        payloads = [
            {
                "_embedded": {"files": [{"fileName": "a.raw"}]},
                "page": {"number": 0, "totalPages": 2},
            },
            {
                "_embedded": {"files": [{"fileName": "b.raw"}]},
                "page": {"number": 1, "totalPages": 2},
            },
        ]
        urls: list[str] = []

        def fake_request(url: str, timeout: int = 120) -> bytes:
            urls.append(url)
            index = len(urls) - 1
            return __import__("json").dumps(payloads[index]).encode("utf-8")

        original = manifest.request_bytes
        try:
            manifest.request_bytes = fake_request
            names = manifest.pride_project_file_names("PXD000612")
        finally:
            manifest.request_bytes = original

        self.assertEqual(names, {"a.raw", "b.raw"})
        self.assertEqual(len(urls), 2)
        self.assertIn("/projects/PXD000612/files?pageSize=100&page=0", urls[0])
        self.assertIn("/projects/PXD000612/files?pageSize=100&page=1", urls[1])
