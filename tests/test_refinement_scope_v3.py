"""Regression coverage for v3 evidence scope versus canonical parameter scope."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from prideqc import __version__
from prideqc.cohort import synthesize_cohort
from prideqc.refinement_adjudication import build_llm_adjudication_request
from prideqc.refinement_packet import build_llm_refinement_packet
from prideqc.sdrf import SDRFDocument
from tests.test_cohort_refinement import _result


class RefinementScopeV3Tests(unittest.TestCase):
    def test_complete_accession_group_is_verified_scope(self) -> None:
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            results = [
                _result(
                    f"run-{index}.raw",
                    precursor=7.8 + index * 0.02,
                    fragment=10.0 + index * 0.02,
                    duration=7200,
                    mass_shift=14.0155 + index * 1e-5,
                    candidates=[("UniMod:34", "Methylation")],
                )
                for index in range(4)
            ]
            path = root / "PXDTEST.sdrf.tsv"
            path.write_text(
                "comment[data file]\tcomment[precursor mass tolerance]"
                "\tcomment[fragment mass tolerance]\tcomment[modification parameters]\n"
                + "".join(
                    f"run-{index}.raw\tnot available\tnot available"
                    "\tNT=Oxidation;AC=UniMod:35\n"
                    for index in range(4)
                ),
                encoding="utf-8",
            )
            document = SDRFDocument.read(path)
            synthesis = synthesize_cohort(results, project_accession="PXDTEST")
            self.assertEqual(len(synthesis.groups), 1)
            packet = build_llm_refinement_packet(
                document,
                results,
                synthesis,
                project_accession="PXDTEST",
                sdrf_path=path,
                prideqc_version=__version__,
            )

        for decision in packet["decision_candidates"]:
            self.assertEqual(decision["parameter_scope_status"], "accession-complete")
            self.assertEqual(decision["target_runs"], decision["parameter_scope_runs"])
            self.assertEqual(decision["target_rows"], decision["parameter_scope_rows"])
            self.assertTrue(
                set(decision["evidence_runs"]).issubset(
                    set(decision["parameter_scope_runs"])
                )
            )

        request = build_llm_adjudication_request(packet)
        for decision in request["decisions"]:
            self.assertEqual(decision["target_runs"], decision["parameter_scope_runs"])
            self.assertIn("evidence_runs", decision)
            self.assertIn("parameter_scope_basis", decision)

    def test_split_cohort_is_explicitly_unverified_for_search_parameter_writeback(self) -> None:
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            results = [
                *[
                    _result(
                        f"short-{index}.raw",
                        precursor=7.0 + index * 0.05,
                        fragment=10.0 + index * 0.05,
                        duration=7200 + index,
                        mass_shift=14.0155 + index * 1e-5,
                        candidates=[("UniMod:34", "Methylation")],
                    )
                    for index in range(4)
                ],
                *[
                    _result(
                        f"long-{index}.raw",
                        precursor=7.4 + index * 0.05,
                        fragment=12.0 + index * 0.05,
                        duration=15900 + index,
                        mass_shift=14.0155 + index * 1e-5,
                        candidates=[("UniMod:34", "Methylation")],
                    )
                    for index in range(6)
                ],
            ]
            path = root / "PXDTEST.sdrf.tsv"
            path.write_text(
                "comment[data file]\tcomment[precursor mass tolerance]"
                "\tcomment[fragment mass tolerance]\tcomment[modification parameters]\n"
                + "".join(
                    f"{result.input_path.name}\tnot available\tnot available"
                    "\tNT=Oxidation;AC=UniMod:35\n"
                    for result in results
                ),
                encoding="utf-8",
            )
            document = SDRFDocument.read(path)
            synthesis = synthesize_cohort(results, project_accession="PXDTEST")
            self.assertEqual(len(synthesis.groups), 2)
            packet = build_llm_refinement_packet(
                document,
                results,
                synthesis,
                project_accession="PXDTEST",
                sdrf_path=path,
                prideqc_version=__version__,
            )

        for decision in packet["decision_candidates"]:
            self.assertEqual(
                decision["parameter_scope_status"],
                "cohort-group-unverified",
            )
            self.assertLess(
                len(decision["parameter_scope_runs"]),
                len(results),
            )


if __name__ == "__main__":
    unittest.main()
