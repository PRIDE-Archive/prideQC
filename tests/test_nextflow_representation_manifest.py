"""Fail-closed mixed native-vendor and converted mzML representation tests."""
from __future__ import annotations

import csv
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

BIN = Path(__file__).resolve().parents[1] / "nextflow" / "bin"


def load(name: str):
    spec = importlib.util.spec_from_file_location(name, BIN / f"{name}.py")
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


selector = load("build_representation_manifest")
metadata_probe = load("inspect_pride_raw_metadata")
science_audit = load("audit_science_metrics")
planner = load("plan_manifest")


class RepresentationSelectionTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.inventory = self.root / "inventory.tsv"

    def write(self, rows):
        with self.inventory.open("w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=sorted(selector.REQUIRED), delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)

    def row(self, filename, run="acq-1", *, native="0", evidence="MSRUN:pride-project-MSRun-1", pxd="PXD039992"):
        return {"pxd_accession": pxd, "file_name": filename,
                "file_category": "RAW", "size_bytes": "1000", "logical_run_id": run,
                "equivalence_evidence": evidence, "native_reader_verified": native}

    def test_curated_vendor_and_mzml_select_one_and_preserve_alternate(self):
        self.write([self.row("A.raw", native="1"), self.row("A.mzML")])
        dest = self.root / "out"
        info = selector.build(self.inventory, dest)
        self.assertEqual(info["logical_runs"], 1)
        self.assertEqual(info["suppressed_alternates"], 1)
        with (dest / "file-manifest.tsv").open() as f:
            self.assertEqual(list(csv.DictReader(f, delimiter="\t"))[0]["archive_file"], "A.raw")
        with (dest / "representation-decisions.tsv").open() as f:
            self.assertEqual(len(list(csv.DictReader(f, delimiter="\t"))), 2)
        plan = self.root / "plan"
        self.assertEqual(planner.prepare(dest / "file-manifest.tsv", dest / "accession-reports.tsv", plan, set())["tasks"], 1)
        self.assertFalse(json.loads((dest / "representation-provenance.json").read_text())["study_merge_allowed"])

    def test_curated_metadata_carries_to_planned_task(self):
        self.write([self.row("A.raw", native="1"), self.row("A.mzML")])
        dest = self.root / "out"
        selector.build(self.inventory, dest)
        plan_dir = self.root / "plan"
        planner.prepare(dest / "file-manifest.tsv", dest / "accession-reports.tsv", plan_dir, set())
        task = json.loads(next((plan_dir / "tasks").glob("*.json")).read_text())
        self.assertEqual(task["representation_policy"], "curated_single_primary_v1")
        self.assertEqual(task["logical_run_id"], "acq-1")
        self.assertRegex(task["representation_inventory_sha256"], r"^[a-f0-9]{64}$")

    def test_unverified_native_uses_mzml_if_equivalent(self):
        self.write([self.row("A.raw"), self.row("A.mzML")])
        selector.build(self.inventory, self.root / "out")
        with (self.root / "out/file-manifest.tsv").open() as f:
            self.assertEqual(list(csv.DictReader(f, delimiter="\t"))[0]["archive_file"], "A.mzML")

    def test_unverified_no_equivalence_blocked(self):
        self.write([self.row("A.raw", native="1", evidence=""), self.row("A.mzML", evidence="")])
        with self.assertRaisesRegex(ValueError, "reviewed evidence"):
            selector.build(self.inventory, self.root / "out")
        self.assertFalse((self.root / "out").exists())

    def test_conflicting_equivalence_evidence_blocked(self):
        self.write([self.row("A.raw", native="1"), self.row("A.mzML", evidence="CURATED:other")])
        with self.assertRaisesRegex(ValueError, "Unverified equivalent"):
            selector.build(self.inventory, self.root / "out")

    def test_derived_same_stem_different_run_ids_rejected(self):
        self.write([self.row("A.raw", run="acq-1", native="1"),
                    self.row("A.mzML", run="acq-2")])
        with self.assertRaisesRegex(ValueError, "Same-stem"):
            selector.build(self.inventory, self.root / "out")

    def test_multi_run_project_with_distinct_acquisitions(self):
        self.write([self.row("A.raw", native="1", run="acq-1"),
                    self.row("A.mzML", run="acq-1"),
                    self.row("B.mzML", run="acq-2", evidence="CURATED:other-run")])
        info = selector.build(self.inventory, self.root / "out")
        self.assertEqual(info["logical_runs"], 2)
        plan = planner.prepare(self.root / "out/file-manifest.tsv",
                               self.root / "out/accession-reports.tsv", self.root / "plan", set())
        self.assertEqual(plan["expected_by_accession"]["PXD039992"], 2)

    def test_mixed_manual_manifest_without_review_rejected(self):
        manifest = self.root / "files.tsv"
        reports = self.root / "reports.tsv"
        manifest.write_text("task_id\tpxd_accession\tarchive_file\tsdrf_path\tneeds_file_map\n"
                            "1\tPXD039992\tA.raw\t\t0\n2\tPXD039992\tB.mzML\t\t0\n")
        reports.write_text("report_id\tpxd_accession\texpected_file_tasks\n1\tPXD039992\t2\n")
        with self.assertRaisesRegex(ValueError, "curated manifest"):
            planner.prepare(manifest, reports, self.root / "plan", set())

    def test_same_stem_manual_manifest_rejected(self):
        manifest = self.root / "files.tsv"
        reports = self.root / "reports.tsv"
        manifest.write_text("task_id\tpxd_accession\tarchive_file\tsdrf_path\tneeds_file_map\n"
                            "1\tPXD039992\tA.raw\t\t0\n2\tPXD039992\tA.mzML\t\t0\n")
        reports.write_text("report_id\tpxd_accession\texpected_file_tasks\n1\tPXD039992\t2\n")
        with self.assertRaisesRegex(ValueError, "ambiguous RAW/mzML"):
            planner.prepare(manifest, reports, self.root / "plan", set())

    def test_invalid_inventory_does_not_write_partial_results(self):
        self.write([self.row("A.raw", native="1"), self.row("../escape.mzML")])
        with self.assertRaises(ValueError):
            selector.build(self.inventory, self.root / "out")
        self.assertFalse((self.root / "out").exists())

    def test_scientific_evidence_audit_does_not_claim_metric_validity(self):
        run = self.root / "published"
        unit = run / "PXD039992" / "_run_outputs" / "unit"
        qc = unit / "qc"
        qc.mkdir(parents=True)
        import hashlib
        mz = qc / "one.mzQC"
        mz.write_text(json.dumps({"mzQC": {"version": "1.0.0", "runQualities": [
            {"qualityMetrics": [{"accession": "MS:4000059", "name": "number of MS1 spectra", "value": 1},
                                {"accession": "QCPRIDE:x", "name": "mass error estimator diagnostics", "value": {}}]}
        ]}}))
        (unit / "result.json").write_text(json.dumps({"accession": "PXD039992", "unit_id": "one",
            "status": "success", "mzqc_file": mz.name,
            "mzqc_sha256": hashlib.sha256(mz.read_bytes()).hexdigest(),
            "requested_estimators": {"mass_error": True, "mass_shifts": False}}))
        report = science_audit.audit(run)
        self.assertEqual(report["successful_receipts"], 1)
        self.assertEqual(report["records"][0]["metric_groups_present"], ["mass_error"])
        self.assertFalse(report["scientific_validity_verified"])

    def test_metadata_probe_is_read_only_and_not_claim_complete(self):
        report = metadata_probe.summarize([
            {"fileName": "BNU_TMT_P10_01.raw", "fileSizeBytes": 1000},
            {"fileName": "BNU_TMT_P10_01.mzML", "fileSizeBytes": 900},
        ], "BNU_TMT_P10_01.raw")
        self.assertEqual(report["status"], "metadata_found")
        self.assertFalse(report["open_format_inventory_complete"])
        self.assertEqual(report["raw_category_mzml_same_stem_candidates"],
                         ["BNU_TMT_P10_01.mzML"])

    def test_idempotent_no_overwrites(self):
        self.write([self.row("A.raw", native="1")])
        selector.build(self.inventory, self.root / "out")
        with self.assertRaisesRegex(ValueError, "Preserving"):
            selector.build(self.inventory, self.root / "out")


if __name__ == "__main__":
    unittest.main()

# Verified runtime policy is checked without depending on a Java installation.
class ProductionSafetySourceTests(unittest.TestCase):
    def test_full_study_bool_is_not_sufficient_evidence(self):
        source = (BIN.parent / "workflows" / "prideqc.nf").read_text()
        self.assertIn("if (params.allow_complete_study)", source)
        self.assertIn("Full-study merging is disabled", source)

    def test_active_source_qc_is_throttled_in_codon_config(self):
        config = (BIN.parent / "conf" / "codon_slurm.config").read_text()
        self.assertIn("withName: PRIDEQC_SOURCE_QC", config)
        self.assertIn("maxForks = 2", config)
        self.assertIn("queueSize = 4", config)
