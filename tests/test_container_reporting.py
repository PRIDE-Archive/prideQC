from __future__ import annotations

import unittest
from pathlib import Path


class ContainerReportingTests(unittest.TestCase):
    @staticmethod
    def dockerfile() -> str:
        path = Path(__file__).resolve().parents[1] / "containers" / "Dockerfile"
        return path.read_text(encoding="utf-8")

    def test_pmultiqc_reporting_environment_is_isolated(self) -> None:
        text = self.dockerfile()
        self.assertIn("FROM ${PYTHON_IMAGE} AS reporting-builder", text)
        self.assertIn("/opt/prideqc/.venv", text)
        self.assertIn("/opt/pmultiqc/.venv", text)
        self.assertIn("COPY --from=reporting-builder /opt/pmultiqc/.venv", text)

    def test_runtime_smokes_direct_mzqc_module(self) -> None:
        text = self.dockerfile()
        self.assertIn("/opt/pmultiqc/.venv/bin/multiqc --version", text)
        self.assertIn("from pmultiqc.modules.mzqc import MzQCModule", text)
        self.assertIn("pmultiqc_resolved_sha=", text)

    def test_runtime_includes_process_listing_for_multiqc_flat_rendering(self) -> None:
        text = self.dockerfile()
        self.assertIn("procps", text)
        self.assertIn("command -v ps >/dev/null", text)
        self.assertIn("ps aux >/dev/null", text)


if __name__ == "__main__":
    unittest.main()
