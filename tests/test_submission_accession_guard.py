"""Regression tests for complete project-accession path matching."""

import unittest
from pathlib import Path

from prideqc.submission import _project_accession_guard


class AccessionGuardBoundaryTests(unittest.TestCase):
    def test_partial_directory_accession_is_ignored(self) -> None:
        path = Path("/tmp/PXD2_diagnostics/PXD123456.sdrf.tsv")
        result = _project_accession_guard(path, {"project_accession": "PXD123456"})
        self.assertEqual(result, "PXD123456")

    def test_different_full_accession_is_rejected(self) -> None:
        path = Path("/tmp/PXD999999/PXD123456.sdrf.tsv")
        with self.assertRaisesRegex(ValueError, "PXD999999"):
            _project_accession_guard(path, {"project_accession": "PXD123456"})


if __name__ == "__main__":
    unittest.main()
