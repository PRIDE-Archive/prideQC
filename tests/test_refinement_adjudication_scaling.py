from __future__ import annotations

import json
import unittest
from pathlib import Path

import jsonschema

from prideqc.refinement_adjudication import (
    build_llm_adjudication_request,
    validate_llm_adjudication_request,
)


class RefinementAdjudicationScalingTests(unittest.TestCase):
    def _empty_packet(self) -> dict[str, object]:
        return {
            "schema_version": "prideqc-llm-refinement-packet-v1",
            "packet_scope": "accession-sdrf",
            "project_accession": "PXDZERO",
            "provenance": {
                "prideqc_version": "0.2.0",
                "input_mode": "sdrf-backed",
            },
            "sdrf": {
                "available": True,
                "source_name": "PXDZERO.sdrf.tsv",
                "sha256": "a" * 64,
                "row_count": 1,
                "target_columns": {},
            },
            "runs": [
                {
                    "run_id": "run.raw",
                    "experiment_group": "Experiment group 1",
                }
            ],
            "experiment_groups": [
                {
                    "group_id": "Experiment group 1",
                    "decision_ids": [],
                }
            ],
            "ptm_context": [],
            "decision_candidates": [],
            "llm_contract": {},
            "limitations": {},
        }

    def test_zero_actionable_decisions_are_valid_request(self) -> None:
        request = build_llm_adjudication_request(self._empty_packet())
        validate_llm_adjudication_request(request)

        self.assertEqual(request["project_accession"], "PXDZERO")
        self.assertEqual(request["input_mode"], "sdrf-backed")
        self.assertEqual(request["decisions"], [])

        schema = json.loads(
            Path("src/prideqc/data/llm-adjudication-request-v1.schema.json").read_text(
                encoding="utf-8"
            )
        )
        jsonschema.validate(request, schema)


if __name__ == "__main__":
    unittest.main()
