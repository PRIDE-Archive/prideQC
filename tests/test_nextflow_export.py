"""Offline Mongo staging: validated, no network or MongoDB writes."""

import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "nextflow/bin/export_mongo_ndjson.py"
SPEC = importlib.util.spec_from_file_location("prideqc_nf_mongo", SCRIPT)
MONGO = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MONGO)


class MongoExportTests(unittest.TestCase):
    def test_complete_revision_to_individual_metrics(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            mzqc = root / "PXD041271.mzQC"
            data = {
                "mzQC": {
                    "runQualities": [
                        {
                            "metadata": {"label": "a.raw", "inputFiles": [{"name": "a.raw"}]},
                            "qualityMetrics": [
                                {"accession": "QCPRIDE:1", "name": "signal", "value": 5}
                            ],
                        }
                    ]
                }
            }
            mzqc.write_text(json.dumps(data))
            obo = root / "PXD041271.obo"
            obo.write_text("format-version: 1.2\n")
            manifest = {
                "accession": "PXD041271",
                "status": "complete",
                "expected_runs": 1,
                "merged_runs": 1,
                "mzqc_sha256": MONGO.sha256(mzqc),
                "obo_sha256": MONGO.sha256(obo),
            }
            sidecar = root / "PXD041271.manifest.json"
            sidecar.write_text(json.dumps(manifest))
            receipt = MONGO.export(mzqc, root / "out")
            self.assertEqual(receipt["metric_documents"], 1)
            self.assertFalse(receipt["mongo_write_performed"])
            record = json.loads((root / "out/PXD041271.metrics.ndjson").read_text())
            self.assertEqual(record["metric_value"], 5)
            manifest["status"] = "count_unverified"
            sidecar.write_text(json.dumps(manifest))
            with self.assertRaisesRegex(ValueError, "incomplete"):
                MONGO.export(mzqc, root / "not-ready")
