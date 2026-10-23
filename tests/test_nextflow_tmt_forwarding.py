"""Regression checks: both TMT collectors must reach `prideqc analyze`."""
from __future__ import annotations

import importlib.util
import json
import subprocess
import tempfile
import types
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]


def load_analyzer():
    filename = ROOT / "nextflow/bin/analyze_unit.py"
    spec = importlib.util.spec_from_file_location("prideqc_test_analyze_unit", filename)
    if spec is None or spec.loader is None:
        raise RuntimeError("Cannot load source analyzer")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


class TMTForwardingTests(unittest.TestCase):
    def test_nextflow_config_defaults_off(self):
        config = (ROOT / "nextflow/nextflow.config").read_text()
        self.assertIn("tmt_reporter_evidence = false", config)
        self.assertIn("tmt_profile_reporter_evidence = false", config)

    def test_both_nextflow_execution_paths_forward_flags(self):
        for rel in ("modules/local/source_qc/main.nf", "modules/local/analyze/main.nf"):
            with self.subTest(rel=rel):
                source = (ROOT / "nextflow" / rel).read_text()
                self.assertIn("params.tmt_reporter_evidence ? '--tmt-reporter-evidence'", source)
                self.assertIn("params.tmt_profile_reporter_evidence ? '--tmt-profile-reporter-evidence'", source)

    def test_python_orchestrators_accept_and_forward_both_flags(self):
        for rel in ("source_qc.py", "analyze_unit.py", "prideqc_run_file.py"):
            with self.subTest(rel=rel):
                source = (ROOT / "nextflow/bin" / rel).read_text()
                self.assertIn('ap.add_argument("--tmt-reporter-evidence", action="store_true")', source)
                self.assertIn('ap.add_argument("--tmt-profile-reporter-evidence", action="store_true")', source)
                self.assertIn('"--tmt-reporter-evidence"', source)
                self.assertIn('"--tmt-profile-reporter-evidence"', source)

    def test_source_qc_passes_tmt_options_and_records_requested_flags(self):
        import importlib.util
        import sys
        analyzer = types.ModuleType("analyze_unit")
        prepmod = types.ModuleType("prepare_archive")
        for enabled in (False, True):
            with self.subTest(enabled=enabled):
                seen = []
                with tempfile.TemporaryDirectory() as tmp:
                    base = Path(tmp)
                    task_file = base / "task.json"
                    task_file.write_text(json.dumps({"accession": "PXD999999", "task_id": 1,
                                                    "archive_file": "dummy.raw", "expected_runs": 1}))
                    def fake_prepare(task_file, protocol, prepared, **kwargs):
                        unit = prepared / "units" / "one"
                        unit.mkdir(parents=True, exist_ok=True)
                        (unit / "unit.json").write_text(json.dumps({"unit_id": "one"}))
                        return {"status": "success"}
                    def fake_analyze(unit, output, *, _seen=seen, **kwargs):
                        _seen.append(dict(kwargs))
                        output.mkdir(parents=True, exist_ok=True)
                        return {"status": "success"}
                    prepmod.prepare = fake_prepare
                    analyzer.analyze = fake_analyze
                    with patch.dict(sys.modules, {"prepare_archive": prepmod, "analyze_unit": analyzer}):
                        spec = importlib.util.spec_from_file_location(
                            "prideqc_test_source_qc", ROOT / "nextflow/bin/source_qc.py")
                        self.assertIsNotNone(spec)
                        self.assertIsNotNone(spec.loader)
                        mod = importlib.util.module_from_spec(spec)
                        spec.loader.exec_module(mod)
                    with patch.object(Path, "cwd", return_value=base):
                        result = mod.run(task_file, "ftp", base / "scratch", base / "out",
                                         tmt_reporter=enabled, tmt_profile=enabled)
                    self.assertEqual(result["status"], "success", result)
                    self.assertEqual(len(seen), 1)
                    self.assertEqual(seen[0]["tmt_reporter"], enabled)
                    self.assertEqual(seen[0]["tmt_profile"], enabled)
                    self.assertEqual(result["requested_estimators"]["tmt_reporter_evidence"], enabled)
                    self.assertEqual(result["requested_estimators"]["tmt_profile_reporter_evidence"], enabled)
                    self.assertFalse((base / "scratch").exists())

    def test_analyzer_enabled_and_disabled_invocations(self):
        analyzer = load_analyzer()
        for tmt_enabled in (False, True):
            with self.subTest(tmt_enabled=tmt_enabled):
                with tempfile.TemporaryDirectory() as tmp:
                    p = Path(tmp)
                    unit = p / "unit"
                    unit.mkdir()
                    (unit / "sample.raw").write_bytes(b"synthetic-only")
                    (unit / "unit.json").write_text(json.dumps({
                        "accession": "PXD999999", "task_id": 1,
                        "source_archive": "sample.raw", "unit_id": "unit-1",
                        "expected_file_tasks": 1, "discovery_status": "success",
                        "input_basename": "sample.raw",
                    }))
                    captured = []
                    def fake_run(cmd, *, _captured=captured, _base=p, **kwargs):
                        _captured.append(list(cmd))
                        qc = _base / "out" / "qc"
                        qc.mkdir(parents=True, exist_ok=True)
                        (qc / "manifest.json").write_text(json.dumps({
                            "files": [{"status": "success", "error": None}]
                        }))
                        (qc / "sample.raw.mzQC").write_text(json.dumps({"mzQC": {
                            "version": "1.0.0",
                            "runQualities": [{"qualityMetrics": [
                                {"accession": "MS:4000059",
                                 "name": "number of MS1 spectra", "value": 1},
                                {"accession": "MS:4000060",
                                 "name": "number of MS2 spectra", "value": 2},
                            ]}]
                        }}))
                        (qc / "sample.raw.obo").write_text("TEST")
                        return subprocess.CompletedProcess(cmd, 0, "", "")
                    with patch.object(analyzer.subprocess, "run", side_effect=fake_run):
                        result = analyzer.analyze(unit, p / "out",
                                                  tmt_reporter=tmt_enabled,
                                                  tmt_profile=tmt_enabled)
                    self.assertEqual(result["status"], "success", result)
                    self.assertEqual(len(captured), 1)
                    self.assertEqual("--tmt-reporter-evidence" in captured[0], tmt_enabled)
                    self.assertEqual("--tmt-profile-reporter-evidence" in captured[0], tmt_enabled)


if __name__ == "__main__":
    unittest.main()
