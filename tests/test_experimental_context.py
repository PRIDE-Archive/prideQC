"""Conservative acquisition and labelling annotation contracts."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from prideqc.annotations import TechnicalAnnotator, labeling_candidates
from prideqc.metrics import RunSummary
from prideqc.models import AnalysisResult, Annotation, CVTerm, EvidenceKind, RunMetadata
from prideqc.mzqc import MzQCWriter, annotation_definition
from prideqc.pipeline import Analyzer
from prideqc.readers import read_header
from tests.helpers import MemoryReader, spectrum


class ExperimentalContextTests(unittest.TestCase):
    def test_wide_isolation_inference_includes_pride_term_and_support(self):
        spectra = [spectrum(i, [10], 2, charge=2, width=20) for i in range(100)]
        result = Analyzer(MemoryReader(spectra)).analyze("dia.mzML")
        annotations = {a.field: a for a in result.annotations}
        acq = annotations["acquisition_method"]
        evidence = annotations["acquisition_method_evidence"]
        self.assertEqual(acq.kind, EvidenceKind.INFERRED)
        self.assertEqual(acq.sdrf_value, "NT=Data-independent acquisition;AC=PRIDE:0000450")
        self.assertEqual(evidence.value["cv_term"]["accession"], "PRIDE:0000450")
        self.assertEqual(evidence.value["supporting_scans"], 100)
        self.assertEqual(evidence.value["assessed_ms2_scans"], 100)
        self.assertIsNone(evidence.sdrf_value)
        with tempfile.TemporaryDirectory() as tmp:
            mzqc = MzQCWriter().build(result, Path(tmp) / "evidence.obo")
            named = {m["name"]: m for m in mzqc["mzQC"]["runQualities"][0]["qualityMetrics"]}
            self.assertEqual(
                named["acquisition method inference evidence"]["value"]["cv_term"]["accession"],
                "PRIDE:0000450",
            )

    def test_ambiguous_targeted_pattern_does_not_automatically_become_prm_or_srm(self):
        target_mzs = [400, 500, 600, 700]
        spectra = []
        rt = 0
        for _ in range(25):
            spectra.append(spectrum(rt, [10], 1))
            rt += 1
            for target in target_mzs:
                spectra.append(spectrum(rt, [10], 2, charge=2, width=2.0,
                                        precursor_mz=target))
                rt += 1
        result = Analyzer(MemoryReader(spectra)).analyze("targeted.mzML")
        named = {a.field: a for a in result.annotations}
        self.assertIsNone(named["acquisition_method"].value)
        self.assertNotIn("acquisition_method_evidence", named)
        self.assertIsNone(named["labeling_strategy"].value)
        self.assertIsNone(named["labeling_strategy"].sdrf_value)
        self.assertIsNone(named["chromatography_type"].value)

    def test_recorded_acquisition_metadata_survives_result_roundtrip(self):
        metadata = RunMetadata(acquisition_methods=[
            CVTerm("PRIDE:0000630", "Selected reaction monitoring")
        ])
        result = Analyzer(MemoryReader([spectrum(0, [10])], metadata)).analyze("srm.mzML")
        restored = AnalysisResult.from_dict(result.to_dict())
        self.assertEqual(restored.metadata.acquisition_methods,
                         [CVTerm("PRIDE:0000630", "Selected reaction monitoring")])
        self.assertEqual(next(a for a in restored.annotations
                              if a.field == "acquisition_method").kind, EvidenceKind.OBSERVED)

    def test_explicit_prm_is_recorded_only_with_exact_header_cv_term(self):
        metadata = RunMetadata(acquisition_methods=[
            CVTerm("PRIDE:0000629", "Parallel reaction monitoring")
        ])
        result = Analyzer(MemoryReader([], metadata)).analyze("prm.mzML")
        named = {a.field: a for a in result.annotations}
        self.assertEqual(named["acquisition_method"].kind, EvidenceKind.OBSERVED)
        self.assertEqual(named["acquisition_method"].sdrf_value,
                         "NT=Parallel reaction monitoring;AC=PRIDE:0000629")
        self.assertEqual(named["acquisition_method_evidence"].value["classification"],
                         "OBSERVED_MZML_FILECONTENT")

    def test_recorded_prm_conflicting_with_wide_window_dia_abstains(self):
        metadata = RunMetadata(acquisition_methods=[
            CVTerm("PRIDE:0000629", "Parallel reaction monitoring")
        ])
        spectra = [spectrum(i, [10], 2, charge=2, width=20) for i in range(100)]
        result = Analyzer(MemoryReader(spectra, metadata)).analyze("conflict.mzML")
        named = {a.field: a for a in result.annotations}
        self.assertIsNone(named["acquisition_method"].value)
        self.assertIsNone(named["acquisition_method"].sdrf_value)
        self.assertEqual(named["acquisition_method_conflict"].value["inferred"],
                         "Data-independent acquisition")
        self.assertNotIn("acquisition_method_evidence", named)

    def test_multiple_explicit_methods_are_preserved_without_auto_selection(self):
        metadata = RunMetadata(acquisition_methods=[
            CVTerm("PRIDE:0000629", "Parallel reaction monitoring"),
            CVTerm("PRIDE:0000630", "Selected reaction monitoring"),
        ])
        result = Analyzer(MemoryReader([], metadata)).analyze("mixed.mzML")
        named = {a.field: a for a in result.annotations}
        self.assertIsNone(named["acquisition_method"].value)
        self.assertIsNone(named["acquisition_method"].sdrf_value)
        self.assertEqual(len(named["acquisition_recorded_terms"].value), 2)
        self.assertNotIn("acquisition_method_evidence", named)

    def test_mzml_header_only_exact_acquisition_cv_term_is_read(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "header.mzML"
            path.write_text("""<?xml version="1.0" encoding="utf-8"?>
<mzML xmlns="http://psi.hupo.org/ms/mzml"><fileDescription><fileContent>
<cvParam accession="PRIDE:0000630" name="Selected reaction monitoring"/>
<cvParam accession="PRIDE:9999999" name="untrusted acquisition label"/>
</fileContent></fileDescription><run id="test"><spectrumList count="0"/>
</run></mzML>""", encoding="utf-8")
            metadata = read_header(path)
            self.assertEqual(metadata.acquisition_methods,
                             [CVTerm("PRIDE:0000630", "Selected reaction monitoring")])

    def test_ionization_is_preserved_from_supplied_cv_not_guessed(self):
        metadata = RunMetadata(ionization=[CVTerm("MS:1000073", "electrospray ionization")])
        run = RunSummary()
        named = {a.field: a for a in TechnicalAnnotator().annotate(metadata, run)}
        self.assertEqual(named["ionization"].kind, EvidenceKind.OBSERVED)
        self.assertEqual(named["ionization"].value[0]["accession"], "MS:1000073")
        self.assertIsNone(named["chromatography_type"].value)
        result = Analyzer(MemoryReader([spectrum(0, [10])], metadata)).analyze("esi.mzML")
        with tempfile.TemporaryDirectory() as tmp:
            run = MzQCWriter().build(result, Path(tmp) / "ion.obo")["mzQC"]["runQualities"][0]
            named_metrics = {m["name"]: m["value"] for m in run["qualityMetrics"]}
            self.assertEqual(named_metrics["ionization source CV evidence"][0]["accession"],
                             "MS:1000073")

    def test_tmt_candidates_do_not_assert_plex_or_sample_label(self):
        evidence = Annotation(
            "tmt_reporter_evidence",
            {"levels": {"MS2": {"methods": {
                "TMT6": {"reporter_positive_spectra": 7},
                "TMTpro16": {"reporter_positive_spectra": 6},
            }}}},
            EvidenceKind.INFERRED,
            "test reporter collector",
        )
        result = labeling_candidates([evidence])
        self.assertIsNotNone(result)
        self.assertEqual(result.value["candidate_terms"][0]["cv_term"]["accession"],
                         "PRIDE:0000314")
        self.assertIsNone(result.sdrf_column)
        self.assertIsNone(result.sdrf_value)
        self.assertEqual(result.value["scientific_status"],
                         "CANDIDATE_ONLY_NOT_VERIFIED_SAMPLE_LABEL")
        self.assertIsNotNone(annotation_definition(result))

    def test_negative_or_unrun_reporter_screen_does_not_prove_label_free(self):
        empty = Annotation("tmt_reporter_evidence", {"levels": {}},
                           EvidenceKind.OBSERVED, "test")
        self.assertIsNone(labeling_candidates([empty]))
        self.assertIsNone(labeling_candidates([]))

    def test_diagnostic_tmt_family_is_candidate_only(self):
        evidence = Annotation("signature_TMT_family", {"matching_spectra": 10},
                              EvidenceKind.INFERRED, "diagnostic ion collector",
                              support=10, total=100)
        result = labeling_candidates([evidence])
        self.assertIsNotNone(result)
        self.assertEqual(result.value["candidate_terms"][0]["cv_term"]["accession"],
                         "PRIDE:0000314")
        self.assertIsNone(result.sdrf_value)

    def test_itraq_candidate_has_separate_pride_term(self):
        evidence = Annotation("signature_iTRAQ_family", {"matching_spectra": 20},
                              EvidenceKind.INFERRED, "diagnostic ion collector",
                              support=20, total=100)
        result = labeling_candidates([evidence])
        self.assertIsNotNone(result)
        self.assertEqual(result.value["candidate_terms"][0]["cv_term"]["accession"],
                         "PRIDE:0000313")

    def test_reference_catalog_generated_from_exact_registry(self):
        from scripts.generate_qc_metrics_reference import DEST, REGISTRY, render

        self.assertTrue(REGISTRY.exists())
        self.assertEqual(DEST.read_text(encoding="utf-8"), render())
        self.assertIn("Registered core QC metric definitions: 104", render())
        self.assertIn("`MS2_PrecursorCharge_Fractions`", render())
        self.assertIn("`acquisition_method_evidence`", render())
        self.assertIn("PRIDE:0000315", render())


if __name__ == "__main__":
    unittest.main()
