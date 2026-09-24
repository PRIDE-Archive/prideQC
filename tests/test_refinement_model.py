from __future__ import annotations

import copy
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from prideqc.refinement_adjudication import build_llm_adjudication_request
from prideqc.refinement_model import (
    LocalLlamaCppAdapter,
    _extract_decision_response,
    build_llama_chat_payload,
)


class RefinementModelTests(unittest.TestCase):
    def _request(self) -> dict[str, object]:
        packet: dict[str, object] = {
            "schema_version": "prideqc-llm-refinement-packet-v1",
            "packet_scope": "accession-sdrf",
            "project_accession": "PXDTEST",
            "provenance": {"prideqc_version": "0.2.0"},
            "sdrf": {
                "source_name": "PXDTEST.sdrf.tsv",
                "sha256": "a" * 64,
                "row_count": 1,
                "target_columns": {},
            },
            "runs": [{"run_id": "run.raw", "experiment_group": "Experiment group 1"}],
            "experiment_groups": [
                {
                    "group_id": "Experiment group 1",
                    "decision_ids": ["Experiment group 1:precursor-mass-tolerance"],
                }
            ],
            "ptm_context": [],
            "decision_candidates": [
                {
                    "decision_id": "Experiment group 1:precursor-mass-tolerance",
                    "decision_type": "mass_tolerance",
                    "experiment_group": "Experiment group 1",
                    "target_field": "comment[precursor mass tolerance]",
                    "write_semantics": "fill_or_replace_canonical_value",
                    "target_rows": [2],
                    "target_runs": ["run.raw"],
                    "original": {"rows": [{"row": 2, "values": ["not available"]}]},
                    "candidate_values": ["8 ppm"],
                    "allowed_values": ["8 ppm"],
                    "allowed_decisions": ["accept", "reject", "abstain"],
                    "evidence": {"supporting_runs": 1, "group_runs": 1},
                }
            ],
            "llm_contract": {},
            "limitations": {},
        }
        return build_llm_adjudication_request(packet)

    def _accepted(self, request: dict[str, object]) -> dict[str, object]:
        decision = request["decisions"][0]  # type: ignore[index]
        return {
            "schema_version": "prideqc-llm-refinement-decision-v1",
            "request_id": request["request_id"],
            "project_accession": "PXDTEST",
            "decisions": [
                {
                    "decision_id": decision["decision_id"],  # type: ignore[index]
                    "decision": "accept",
                    "selected_value": "8 ppm",
                    "reason": "The original value is unavailable and the cohort supports 8 ppm.",
                }
            ],
        }

    def test_payload_is_local_deterministic_and_schema_constrained(self) -> None:
        request = self._request()
        first = build_llama_chat_payload(request)
        second = build_llama_chat_payload(copy.deepcopy(request))
        self.assertEqual(first, second)
        self.assertEqual(first["temperature"], 0.0)
        self.assertFalse(first["stream"])
        self.assertEqual(first["max_tokens"], 384)
        self.assertTrue(first["cache_prompt"])
        self.assertFalse(first["chat_template_kwargs"]["enable_thinking"])
        model_input = json.loads(first["messages"][1]["content"])
        self.assertEqual(
            set(model_input),
            {"request_id", "project_accession", "decision"},
        )
        self.assertNotIn("source_packet", model_input)
        self.assertNotIn("sdrf", model_input)
        self.assertNotIn("contract", model_input)
        self.assertEqual(first["response_format"]["type"], "json_object")
        self.assertEqual(
            first["response_format"]["schema"]["properties"]["decisions"]["type"],
            "array",
        )
        self.assertNotIn("tools", first)

    def test_extract_decision_response_requires_json_content(self) -> None:
        request = self._request()
        accepted = self._accepted(request)
        raw = {"choices": [{"message": {"content": json.dumps(accepted)}}]}
        self.assertEqual(_extract_decision_response(raw), accepted)
        with self.assertRaisesRegex(RuntimeError, "non-JSON"):
            _extract_decision_response({"choices": [{"message": {"content": "not json"}}]})

    def test_local_adapter_validates_model_candidate_before_returning(self) -> None:
        request = self._request()
        accepted = self._accepted(request)
        raw = {"choices": [{"message": {"content": json.dumps(accepted)}}]}

        class FakeServer:
            base_url = "http://127.0.0.1:12345"

            def __init__(self, *args: object, **kwargs: object) -> None:
                pass

            def __enter__(self) -> FakeServer:
                return self

            def __exit__(self, *args: object) -> None:
                pass

        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            server = root / "llama-server"
            model = root / "model.gguf"
            server.touch()
            model.touch()
            adapter = LocalLlamaCppAdapter(server_path=server, model_path=model)
            with (
                patch("prideqc.refinement_model._LocalServer", FakeServer),
                patch("prideqc.refinement_model._json_request", return_value=raw) as request_mock,
            ):
                result = adapter.adjudicate(request)
            self.assertEqual(result.decisions, accepted)
            self.assertEqual(result.audit["adapter"], "local-llama.cpp")
            self.assertFalse(result.audit["runtime"]["managed"])
            self.assertIsNone(result.audit["runtime"]["build"])
            self.assertFalse(result.audit["model"]["managed"])
            self.assertIsNone(result.audit["model"]["repository"])
            self.assertFalse(result.audit["model"]["thinking_enabled"])
            payload = request_mock.call_args.args[1]
            self.assertEqual(payload["model"], "model.gguf")

            bad = copy.deepcopy(accepted)
            bad["decisions"][0]["selected_value"] = "9 ppm"  # type: ignore[index]
            bad_raw = {"choices": [{"message": {"content": json.dumps(bad)}}]}
            with (
                patch("prideqc.refinement_model._LocalServer", FakeServer),
                patch("prideqc.refinement_model._json_request", return_value=bad_raw),
            ):
                with self.assertRaisesRegex(ValueError, "exact supplied candidate"):
                    adapter.adjudicate(request)

    def test_local_adapter_splits_accession_request_into_one_decision_calls(self) -> None:
        request = self._request()
        second = copy.deepcopy(request["decisions"][0])  # type: ignore[index]
        second["decision_id"] = "Experiment group 1:fragment-mass-tolerance"
        second["target_field"] = "comment[fragment mass tolerance]"
        second["candidate_values"] = ["11 ppm"]
        request["decisions"].append(second)  # type: ignore[union-attr]

        class FakeServer:
            base_url = "http://127.0.0.1:12345"

            def __init__(self, *args: object, **kwargs: object) -> None:
                pass

            def __enter__(self) -> FakeServer:
                return self

            def __exit__(self, *args: object) -> None:
                pass

        seen_prompt_decisions: list[list[str]] = []

        def fake_json_request(
            url: str,
            payload: dict[str, object] | None,
            timeout: float,
        ) -> dict[str, object]:
            self.assertIsNotNone(payload)
            assert payload is not None
            messages = payload["messages"]
            assert isinstance(messages, list)
            user_message = messages[1]
            assert isinstance(user_message, dict)
            model_input = json.loads(str(user_message["content"]))
            self.assertEqual(
                set(model_input),
                {"request_id", "project_accession", "decision"},
            )
            decision = model_input["decision"]
            decision_id = decision["decision_id"]
            seen_prompt_decisions.append([decision_id])
            selected = decision["candidate_values"][0]
            response = {
                "schema_version": "prideqc-llm-refinement-decision-v1",
                "request_id": request["request_id"],
                "project_accession": "PXDTEST",
                "decisions": [
                    {
                        "decision_id": decision_id,
                        "decision": "accept",
                        "selected_value": selected,
                        "reason": "The supplied evidence supports this candidate.",
                    }
                ],
            }
            return {"choices": [{"message": {"content": json.dumps(response)}}]}

        progress: list[tuple[int, int, str]] = []
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            server = root / "llama-server"
            model = root / "model.gguf"
            server.touch()
            model.touch()
            adapter = LocalLlamaCppAdapter(
                server_path=server,
                model_path=model,
                progress=lambda index, total, decision_id: progress.append(
                    (index, total, decision_id)
                ),
            )
            with (
                patch("prideqc.refinement_model._LocalServer", FakeServer),
                patch("prideqc.refinement_model._json_request", side_effect=fake_json_request),
            ):
                result = adapter.adjudicate(request)

        self.assertEqual(len(result.decisions["decisions"]), 2)
        self.assertEqual(len(seen_prompt_decisions), 2)
        self.assertTrue(all(len(ids) == 1 for ids in seen_prompt_decisions))
        self.assertEqual(result.audit["inference_mode"], "one-decision-per-call")
        self.assertEqual(len(result.audit["decision_calls"]), 2)
        self.assertTrue(
            all("elapsed_seconds" in call for call in result.audit["decision_calls"])
        )
        self.assertEqual([item[:2] for item in progress], [(1, 2), (2, 2)])



if __name__ == "__main__":
    unittest.main()
