from __future__ import annotations

import os
import subprocess
import tempfile
import unittest
from pathlib import Path


class SlurmAccessionReportTests(unittest.TestCase):
    @staticmethod
    def report_wrapper() -> Path:
        return (
            Path(__file__).resolve().parents[1]
            / "scripts"
            / "slurm"
            / "prideqc_accession_report.sbatch"
        )

    @staticmethod
    def report_submitter() -> Path:
        return (
            Path(__file__).resolve().parents[1]
            / "scripts"
            / "slurm"
            / "submit_prideqc_accession_reports.sh"
        )

    def test_report_scripts_have_valid_bash_syntax(self) -> None:
        subprocess.run(["bash", "-n", str(self.report_wrapper())], check=True)
        subprocess.run(["bash", "-n", str(self.report_submitter())], check=True)

    def test_report_wrapper_uses_only_mzqc_module_and_checks_completeness(self) -> None:
        text = self.report_wrapper().read_text(encoding="utf-8")
        self.assertIn("--module mzqc", text)
        self.assertIn("EXPECTED_TASKS", text)
        self.assertIn("FOUND_MZQC", text)
        self.assertIn("PRIDEQC_REPORT_REQUIRE_COMPLETE", text)
        self.assertIn("multiqc_report.html", text)

    def test_study_merge_never_enters_multiqc_discovery(self) -> None:
        text = self.report_wrapper().read_text(encoding="utf-8")
        self.assertIn('PER_RUN_ROOT="$ACCESSION_ROOT/repository-metadata"', text)
        self.assertIn('STUDY_ROOT="$ACCESSION_ROOT/_study"', text)
        self.assertIn('find "$PER_RUN_ROOT"', text)
        self.assertIn('--bind "$PER_RUN_ROOT:/input:ro"', text)
        self.assertIn("prideqc mzqc merge", text)
        self.assertIn('--expected-runs "$EXPECTED_TASKS"', text)
        self.assertIn("MULTIQC_STATUS == 0 && FOUND_MZQC == EXPECTED_TASKS", text)
        self.assertIn("study_merge_status=$MERGE_STATUS", text)
        self.assertNotIn('find "$ACCESSION_ROOT"', text)

    def test_submitter_builds_one_report_task_per_accession(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            manifest = root / "run.tsv"
            manifest.write_text(
                "task_id\tpxd_accession\n1\tPXD000001\n2\tPXD000001\n3\tPXD000002\n",
                encoding="utf-8",
            )
            sif = root / "prideqc.sif"
            sif.write_bytes(b"sif")
            (root / "prideqc.sif.sha256").write_text("dummy  prideqc.sif\n", encoding="utf-8")
            launcher = root / "report.sbatch"
            launcher.write_text("#!/usr/bin/env bash\n", encoding="utf-8")

            fake_bin = root / "bin"
            fake_bin.mkdir()
            sbatch_log = root / "sbatch-args.txt"
            sbatch = fake_bin / "sbatch"
            sbatch.write_text(
                '#!/usr/bin/env bash\nprintf "%s\\n" "$@" > "$SBATCH_LOG"\nprintf "12345\\n"\n',
                encoding="utf-8",
            )
            sbatch.chmod(0o755)

            report_manifest = root / "reports.tsv"
            env = os.environ.copy()
            env.update(
                {
                    "PATH": f"{fake_bin}:{env['PATH']}",
                    "SBATCH_LOG": str(sbatch_log),
                    "PERSIST_ROOT": str(root),
                    "SIF": str(sif),
                    "RUN_NAME": "run",
                    "MANIFEST": str(manifest),
                    "REPORT_MANIFEST": str(report_manifest),
                    "REPORT_LAUNCHER": str(launcher),
                    "DEPENDENCY_JOB_ID": "999",
                }
            )
            subprocess.run(
                ["bash", str(self.report_submitter())],
                env=env,
                check=True,
                capture_output=True,
                text=True,
            )

            self.assertEqual(
                report_manifest.read_text(encoding="utf-8").splitlines(),
                [
                    "report_id\tpxd_accession\texpected_file_tasks",
                    "1\tPXD000001\t2",
                    "2\tPXD000002\t1",
                ],
            )
            args = sbatch_log.read_text(encoding="utf-8")
            self.assertIn("--dependency=afterany:999", args)
            self.assertIn("--array=1-2%2", args)


if __name__ == "__main__":
    unittest.main()
