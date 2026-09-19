from __future__ import annotations

import subprocess
import unittest
from pathlib import Path


class SlurmFileArrayTests(unittest.TestCase):
    @staticmethod
    def wrapper() -> Path:
        return Path(__file__).resolve().parents[1] / "scripts" / "slurm" / "prideqc_file_array.sbatch"

    def test_wrapper_has_valid_bash_syntax(self) -> None:
        subprocess.run(["bash", "-n", str(self.wrapper())], check=True)

    def test_checksum_check_can_be_disabled_and_is_recorded(self) -> None:
        text = self.wrapper().read_text(encoding="utf-8")
        self.assertIn('PRIDEQC_CHECKSUM_CHECK="${PRIDEQC_CHECKSUM_CHECK:-1}"', text)
        self.assertIn('download_flags+=(--no-checksum-check)', text)
        self.assertIn('pride_checksum_check=$PRIDEQC_CHECKSUM_CHECK', text)

    def test_wrapper_supports_metadata_only_manifest_without_fake_sdrf(self) -> None:
        text = self.wrapper().read_text(encoding="utf-8")
        self.assertIn('HAS_SDRF=0', text)
        self.assertIn('if [[ -n "$SDRF_MANIFEST_PATH" ]]; then', text)
        self.assertIn('SAFE_SDRF="repository-metadata"', text)
        self.assertIn('sdrf_flags=()', text)
        self.assertIn('sdrf_flags=(--sdrf /work/input.sdrf.tsv --sdrf-template "$SDRF_TEMPLATE")', text)
        self.assertIn('"${sdrf_flags[@]}"', text)
        self.assertIn('Metadata-only manifest rows cannot require an SDRF file map', text)


if __name__ == "__main__":
    unittest.main()
