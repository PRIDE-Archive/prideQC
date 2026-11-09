from __future__ import annotations

import hashlib
import importlib.util
import json
import tempfile
import threading
import unittest
import unittest.mock
from http.server import BaseHTTPRequestHandler, HTTPServer
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "nextflow/bin/upload_qc_mongo.py"
SPEC = importlib.util.spec_from_file_location("prideqc_upload_qc_mongo", SCRIPT)
UP = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(UP)


def _run(source_file: str = "a.raw", **over) -> dict:
    run = {"_id": hashlib.sha256(f"PXD000001|{source_file}".encode()).hexdigest(),
           "schema_version": 1, "pxd_accession": "PXD000001", "source_file": source_file,
           "source_sha256": "b" * 64, "source_size_bytes": 10, "mzqc_sha256": "a" * 64,
           "mzqc_path": f"PXD000001/{source_file}.mzQC", "prideqc_version": "0.2.0",
           "pyopenms_version": "3.7.0", "mzqc_creation_date": "2026-11-06T00:00:00+00:00",
           "run_tag": "t", "scientific_approval": False,
           "metrics": [{"accession": "MS:4000059", "name": "number of MS1 spectra",
                        "unit": {"accession": "UO:0000189", "name": "count unit"}, "value": 5},
                       {"accession": "QCPRIDE:X", "name": "nested", "unit": None,
                        "value": {"a": [1, 2.5], "b": "x"}}]}
    run.update(over)
    return run


class _Server:
    def __init__(self, status: int = 200):
        self.requests: list = []
        outer = self

        class Handler(BaseHTTPRequestHandler):
            def do_POST(self):  # noqa: N802
                body = json.loads(self.rfile.read(int(self.headers["Content-Length"])))
                outer.requests.append((self.path, self.headers.get("mongo-ws-api-key"), body))
                self.send_response(status)
                self.send_header("Content-Type", "application/json")
                self.end_headers()
                self.wfile.write(json.dumps(body if status == 200 else {"error": "x"}).encode())

            def log_message(self, *args):
                pass

        self.httpd = HTTPServer(("127.0.0.1", 0), Handler)
        self.url = f"http://127.0.0.1:{self.httpd.server_port}/pride/mongo/ws"
        threading.Thread(target=self.httpd.serve_forever, daemon=True).start()

    def close(self):
        self.httpd.shutdown()
        self.httpd.server_close()


class UploadQcMongoTests(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.path = Path(tmp.name) / "runs.ndjson"

    def _write(self, runs):
        self.path.write_text("".join(json.dumps(r) + "\n" for r in runs), encoding="utf-8")

    def test_mapping_keeps_values_and_units(self):
        doc = UP.to_document(_run())
        self.assertEqual(doc["id"], _run()["_id"])
        self.assertEqual(doc["projectAccession"], "PXD000001")
        self.assertEqual(doc["mzqcPath"], "PXD000001/a.raw.mzQC")
        self.assertIs(doc["scientificApproval"], False)
        self.assertEqual(doc["metrics"][0]["unit"]["name"], "count unit")
        self.assertEqual(doc["metrics"][1]["value"], {"a": [1, 2.5], "b": "x"})
        self.assertNotIn("_id", doc)

    def test_rejects_bad_documents(self):
        for bad in (_run(_id="0" * 64), _run(mzqc_sha256="xyz"), _run(metrics=[]),
                    _run(scientific_approval=True)):
            with self.assertRaises(ValueError):
                UP.to_document(bad)

    def test_dry_run_makes_no_network_call(self):
        server = _Server()
        self.addCleanup(server.close)
        self._write([_run("a.raw"), _run("b.raw")])
        self.assertEqual(UP.main([str(self.path), "--base-url", server.url]), 0)
        self.assertEqual(server.requests, [])

    def test_execute_requires_key(self):
        self._write([_run()])
        with unittest.mock.patch.dict("os.environ", {}, clear=True), self.assertRaises(SystemExit):
            UP.main([str(self.path), "--execute"])

    def test_execute_batches_and_filters(self):
        server = _Server()
        self.addCleanup(server.close)
        runs = [_run(f"{i}.raw") for i in range(5)] + [_run("x.raw", pxd_accession="PXD000002",
                                                            _id=UP.run_id("PXD000002", "x.raw"))]
        self._write(runs)
        with unittest.mock.patch.dict("os.environ", {"PRIDE_MONGO_WS_API_KEY": "k"}):
            UP.main([str(self.path), "--base-url", server.url, "--execute", "--batch-size", "2",
                     "--accession", "PXD000001"])
        self.assertEqual([len(b) for _, _, b in server.requests], [2, 2, 1])
        self.assertTrue(all(p == "/pride/mongo/ws/qc-runs" and k == "k" for p, k, _ in server.requests))
        self.assertEqual({d["projectAccession"] for _, _, b in server.requests for d in b}, {"PXD000001"})

    def test_client_error_is_not_retried(self):
        server = _Server(status=400)
        self.addCleanup(server.close)
        self._write([_run()])
        with unittest.mock.patch.dict("os.environ", {"PRIDE_MONGO_WS_API_KEY": "k"}), \
                self.assertRaises(RuntimeError):
            UP.main([str(self.path), "--base-url", server.url, "--execute"])
        self.assertEqual(len(server.requests), 1)

    def test_batch_size_limit(self):
        self._write([_run()])
        with self.assertRaises(SystemExit):
            UP.main([str(self.path), "--batch-size", "21"])


if __name__ == "__main__":
    unittest.main()
