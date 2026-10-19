"""CV replacements must agree between writer and the tracked metric registry."""

from __future__ import annotations

import csv
import tempfile
import unittest
from pathlib import Path

from prideqc.models import Metric
from prideqc.mzqc import MzQCWriter, definition


class SurvivingCvMappingTests(unittest.TestCase):
    def test_official_charge_terms_and_value_shapes(self) -> None:
        registry = Path(__file__).resolve().parents[1] / "docs" / "metrics.tsv"
        with registry.open(encoding="utf-8", newline="") as handle:
            rows = {row["key"]: row for row in csv.DictReader(handle, delimiter="\t")}
        expected = {
            "ChargeMean": "MS:4000173",
            "ChargeMedian": "MS:4000175",
            "ChargeRatio_3over2": "MS:4000169",
            "ChargeRatio_4over2": "MS:4000171",
        }
        self.assertEqual(len(rows), 30)
        for key, accession in expected.items():
            with self.subTest(key=key):
                term = definition(Metric(key, 2))
                self.assertEqual(term.accession, accession)
                self.assertEqual(term.shape, "scalar")
                self.assertEqual(rows[key]["accession"], accession)
                self.assertEqual(rows[key]["shape"], "scalar")
        # The 76 retired metric keys must remain absent.
        retired = Path(__file__).resolve().parents[1] / "docs" / "retired-metrics.tsv"
        with retired.open(encoding="utf-8", newline="") as handle:
            retired_keys = {row["key"] for row in csv.DictReader(handle, delimiter="\t")}
        self.assertEqual(len(retired_keys), 76)
        self.assertTrue(retired_keys.isdisjoint(rows))
        # Official terms belong to PSI-MS; the companion local OBO must not
        # redefine their former QCPRIDE accessions.
        with tempfile.TemporaryDirectory() as temporary:
            companion = Path(temporary) / "run.obo"
            MzQCWriter().write_vocabulary(
                [Metric(key, 2) for key in expected], companion,
            )
            content = companion.read_text(encoding="utf-8")
            for old_accession in (
                "QCPRIDE:E303A6C9571AC1DD9AD6",
                "QCPRIDE:3E06B5777AC79EC2AAC6",
                "QCPRIDE:3263CF93B552C37EE7C0",
                "QCPRIDE:D916EC99E4BD5371D2AE",
            ):
                self.assertNotIn(old_accession, content)


if __name__ == "__main__":
    unittest.main()
