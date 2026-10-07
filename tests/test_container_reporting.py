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


    @staticmethod
    def docker_workflow() -> str:
        path = (
            Path(__file__).resolve().parents[1]
            / ".github"
            / "workflows"
            / "container-docker.yml"
        )
        return path.read_text(encoding="utf-8")

    def test_docker_workflow_resolves_openms_ref_to_immutable_sha(self) -> None:
        text = self.docker_workflow()
        self.assertIn("OPENMS_GIT_REF: 'develop'", text)
        self.assertIn('openms_sha="$(' , text)
        self.assertIn('git ls-remote', text)
        self.assertIn('OPENMS_GIT_REF=${{ steps.buildinfo.outputs.openms_sha }}', text)
        self.assertIn('OPENMS_GIT_REQUESTED_REF=${{ env.OPENMS_GIT_REF }}', text)
        self.assertIn('openms_resolved_sha=${{ steps.buildinfo.outputs.openms_sha }}', text)

    def test_container_version_provenance_uses_hatch_vcs(self) -> None:
        workflow = self.docker_workflow()
        build_script = (
            Path(__file__).resolve().parents[1] / "containers" / "build_docker.sh"
        ).read_text(encoding="utf-8")
        self.assertIn('project_version="$(uvx --from hatch --with hatch-vcs hatch version)"', workflow)
        self.assertIn('PROJECT_VERSION="$(uvx --from hatch --with hatch-vcs hatch version)"', build_script)
        self.assertIn('astral-sh/setup-uv@v10.2.0', workflow)
        self.assertNotIn("/^version = ", workflow)
        self.assertNotIn("/^version = ", build_script)

    def test_runtime_includes_process_listing_for_multiqc_flat_rendering(self) -> None:
        text = self.dockerfile()
        self.assertIn("procps", text)
        self.assertIn("command -v ps >/dev/null", text)
        self.assertIn("ps aux >/dev/null", text)

    @staticmethod
    def accession_report_script() -> str:
        path = (
            Path(__file__).resolve().parents[1]
            / "scripts"
            / "slurm"
            / "prideqc_accession_report.sbatch"
        )
        return path.read_text(encoding="utf-8")

    def test_accession_report_uses_strict_reporting_environment(self) -> None:
        text = self.accession_report_script()
        self.assertIn("/opt/pmultiqc/.venv/bin/multiqc", text)
        self.assertIn("--strict", text)
        self.assertIn("--interactive", text)
        self.assertIn('--title "$PXD — mzQC Quality Control"', text)


if __name__ == "__main__":
    unittest.main()
