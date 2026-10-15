"""Synthetic mzQC 1.0 merger regression fixtures (not actual Codon data)."""

from __future__ import annotations

import copy
import hashlib
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from prideqc.mzqc_merge import merge_mzqc


class MergeMzQCTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.input_dir = self.root / "repository-metadata"
        self.output = self.root / "_study" / "PXD041271.mzQC"

    def make_run(
        self,
        index: int,
        *,
        basename: str | None = None,
        cv_version: str = "1",
        obo_definition: str = "Spectra evidence",
    ) -> Path:
        name = basename or f"run{index:03}.raw"
        path = self.input_dir / f"{index:03}_{name}" / f"{name}.mzQC"
        path.parent.mkdir(parents=True, exist_ok=True)
        data = {
            "mzQC": {
                "version": "1.0.0",
                "creationDate": "2026-10-15T12:00:00+00:00",
                "runQualities": [
                    {
                        "metadata": {
                            "label": name,
                            "inputFiles": [
                                {
                                    "name": name,
                                    "location": f"file:///input/{name}",
                                    "fileFormat": {
                                        "accession": "MS:1000563",
                                        "name": "Thermo RAW format",
                                    },
                                    "fileProperties": [
                                        {
                                            "accession": "MS:1001919",
                                            "name": "ProteomeXchange accession number",
                                            "value": "PXD041271",
                                        }
                                    ],
                                }
                            ],
                            "analysisSoftware": [
                                {
                                    "accession": "QCPRIDE:SOFTWARE",
                                    "name": "prideqc",
                                    "version": "0.2",
                                }
                            ],
                        },
                        "qualityMetrics": [
                            {"accession": "QCPRIDE:METRIC", "name": "signal", "value": index},
                            {
                                "accession": "MS:4000059",
                                "name": "number of MS1 spectra",
                                "value": 1,
                            },
                        ],
                    }
                ],
                "controlledVocabularies": [
                    {"name": "PSI-MS", "uri": "https://example.org/psi.obo"},
                    {
                        "name": "prideqc local vocabulary",
                        "version": cv_version,
                        "uri": f"file:///temporary/{index}.obo",
                    },
                ],
            }
        }
        path.write_text(json.dumps(data), encoding="utf-8")
        path.with_suffix(".obo").write_text(
            "format-version: 1.2\ndata-version: 1\nontology: prideqc\n"
            "\n[Term]\nid: QCPRIDE:SOFTWARE\nname: prideqc\n"
            "\n[Term]\nid: QCPRIDE:SOURCEFILE\nname: original vendor file\n"
            f'\n[Term]\nid: QCPRIDE:METRIC\nname: signal\ndef: "{obo_definition}" []\n',
            encoding="utf-8",
        )
        return path

    def test_one_and_67_runs_are_lossless_valid_and_idempotent(self) -> None:
        from jsonschema import Draft7Validator, FormatChecker

        schema = json.loads((Path(__file__).parent / "data" / "mzqc_schema.json").read_text())
        for n in (1, 67):
            with self.subTest(n=n):
                for index in range(1, n + 1):
                    self.make_run(index)
                member_runs = [
                    json.loads(path.read_text())["mzQC"]["runQualities"][0]
                    for path in sorted(self.input_dir.rglob("*.mzQC"))
                ]
                manifest = merge_mzqc(
                    self.input_dir,
                    self.output,
                    accession="PXD041271",
                    expected_runs=n,
                )
                doc = json.loads(self.output.read_text())
                self.assertEqual(doc["mzQC"]["runQualities"], member_runs)
                self.assertNotIn("setQualities", doc["mzQC"])
                self.assertEqual(manifest["merged_runs"], n)
                self.assertEqual(manifest["status"], "complete")
                errors = list(
                    Draft7Validator(schema, format_checker=FormatChecker()).iter_errors(doc)
                )
                self.assertFalse(errors, errors)
                self.assertEqual(
                    [cv["name"] for cv in doc["mzQC"]["controlledVocabularies"]],
                    ["PSI-MS", "prideqc local vocabulary"],
                )
                self.assertEqual(
                    next(
                        cv["uri"]
                        for cv in doc["mzQC"]["controlledVocabularies"]
                        if cv["name"] == "prideqc local vocabulary"
                    ),
                    self.output.with_suffix(".obo").absolute().as_uri(),
                )
                before = self.output.stat().st_mtime_ns
                self.assertEqual(
                    merge_mzqc(self.input_dir, self.output, accession="PXD041271", expected_runs=n),
                    manifest,
                )
                self.assertEqual(self.output.stat().st_mtime_ns, before)
                self.assertEqual(
                    hashlib.sha256(self.output.read_bytes()).hexdigest(),
                    json.loads(self.output.with_suffix(".manifest.json").read_text())[
                        "mzqc_sha256"
                    ],
                )

    def test_partial_does_not_publish_or_replace(self) -> None:
        self.make_run(1)
        merge_mzqc(self.input_dir, self.output, accession="PXD041271", expected_runs=1)
        original = self.output.read_bytes()
        with self.assertRaisesRegex(ValueError, "Incomplete accession"):
            merge_mzqc(self.input_dir, self.output, accession="PXD041271", expected_runs=2)
        self.assertEqual(self.output.read_bytes(), original)

    def test_missing_companion_obo_is_rejected(self) -> None:
        self.make_run(1).with_suffix(".obo").unlink()
        with self.assertRaisesRegex(ValueError, "Missing companion OBO"):
            merge_mzqc(self.input_dir, self.output, accession="PXD041271")
        self.assertFalse(self.output.exists())

    def test_duplicate_input_names_are_rejected(self) -> None:
        self.make_run(1, basename="duplicate.raw")
        self.make_run(2, basename="duplicate.raw")
        with self.assertRaisesRegex(ValueError, "Duplicate mzQC run label"):
            merge_mzqc(self.input_dir, self.output, accession="PXD041271")

    def test_conflicting_terms_are_rejected(self) -> None:
        self.make_run(1)
        self.make_run(2, obo_definition="Conflicting definition")
        with self.assertRaisesRegex(ValueError, "Conflicting QCPRIDE OBO definition"):
            merge_mzqc(self.input_dir, self.output, accession="PXD041271")

    def test_invalid_json_and_multiple_runqualites_are_rejected(self) -> None:
        path = self.make_run(1)
        path.write_text("{bad", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "Invalid JSON"):
            merge_mzqc(self.input_dir, self.output, accession="PXD041271")
        self.make_run(1)
        data = json.loads(path.read_text())
        data["mzQC"]["runQualities"].append(copy.deepcopy(data["mzQC"]["runQualities"][0]))
        path.write_text(json.dumps(data))
        with self.assertRaisesRegex(ValueError, "exactly one runQuality"):
            merge_mzqc(self.input_dir, self.output, accession="PXD041271")

    def test_mismatched_project_and_cv_are_rejected(self) -> None:
        one = self.make_run(1)
        self.make_run(2)
        data = json.loads(one.read_text())
        data["mzQC"]["runQualities"][0]["metadata"]["inputFiles"][0]["fileProperties"][0][
            "value"
        ] = "PXD999999"
        one.write_text(json.dumps(data))
        with self.assertRaisesRegex(ValueError, "accession mismatch"):
            merge_mzqc(self.input_dir, self.output, accession="PXD041271")
        self.make_run(1)
        data = json.loads(one.read_text())
        data["mzQC"]["controlledVocabularies"][0]["uri"] = "https://other.org/psi.obo"
        one.write_text(json.dumps(data))
        with self.assertRaisesRegex(ValueError, "Conflicting controlled vocabulary"):
            merge_mzqc(self.input_dir, self.output, accession="PXD041271")

    def test_local_vocab_version_setqualities_and_output_inside_inputs_are_rejected(self) -> None:
        path = self.make_run(1, cv_version="2")
        with self.assertRaisesRegex(ValueError, "Unsupported prideqc local vocabulary version"):
            merge_mzqc(self.input_dir, self.output, accession="PXD041271")
        self.make_run(1)
        data = json.loads(path.read_text())
        data["mzQC"]["setQualities"] = []
        path.write_text(json.dumps(data))
        with self.assertRaisesRegex(ValueError, "Unexpected setQualities"):
            merge_mzqc(self.input_dir, self.output, accession="PXD041271")
        self.make_run(1)
        with self.assertRaisesRegex(ValueError, "outside the per-run input"):
            merge_mzqc(
                self.input_dir,
                self.input_dir / "_study" / "PXD041271.mzQC",
                accession="PXD041271",
            )

    def test_cli(self) -> None:
        self.make_run(1)
        done = subprocess.run(
            [
                sys.executable,
                "-m",
                "prideqc",
                "mzqc",
                "merge",
                "--input-dir",
                str(self.input_dir),
                "--output",
                str(self.output),
                "--accession",
                "PXD041271",
                "--expected-runs",
                "1",
            ],
            capture_output=True,
            text=True,
        )
        self.assertEqual(done.returncode, 0, done.stderr)
        self.assertIn("Merged 1 runs", done.stdout)


if __name__ == "__main__":
    unittest.main()
