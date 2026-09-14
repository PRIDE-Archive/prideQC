"""Evidence, SDRF, serialization and orchestration contracts."""

import importlib.util
import io
import json
import math
import subprocess
import sys
import tempfile
import unittest
from contextlib import redirect_stderr
from dataclasses import replace
from pathlib import Path
from unittest.mock import patch

import numpy as np

from prideqc.annotations import DiagnosticIonCollector
from prideqc.cli import main
from prideqc.conversion import ExternalConverter
from prideqc.io import atomic_text, json_safe
from prideqc.mass_error import RepeatSpectrumMassErrorCollector
from prideqc.models import Annotation, EvidenceKind, Metric, RunMetadata
from prideqc.mzqc import MzQCWriter, definition
from prideqc.pipeline import Analyzer, FileOutcome, Workflow, WorkflowOptions
from prideqc.readers import read_header
from prideqc.sdrf import SDRFDocument
from tests.helpers import FUSION, MemoryReader, spectrum


def analyze(spectra=None, name="run.mzML", metadata=None):
    reader = MemoryReader(spectra or [spectrum(0, [10]), spectrum(10, [20])], metadata)
    return Analyzer(reader).analyze(name)


class AnnotationTests(unittest.TestCase):
    def test_one_pass_shared_by_metrics_and_diagnostics(self):
        collector = DiagnosticIonCollector()
        reader = MemoryReader(
            [spectrum(0, [10, 10, 10], 2, mz=[126.127726, 127.131081, 128.134436])],
        )
        result = Analyzer(reader).analyze("a.mzML", collectors=[collector])
        self.assertEqual(reader.read_count, 1)
        evidence = next(a for a in result.annotations if a.field == "signature_TMT_family")
        self.assertEqual(evidence.support, 1)
        self.assertIsNone(evidence.sdrf_value)


    def test_repeat_spectrum_mass_error_estimator_reports_precision(self):
        collector = RepeatSpectrumMassErrorCollector(
            min_spectrum_pairs=10,
            min_precursor_clusters=1,
            min_fragment_pairs=80,
        )
        spectra = []
        base_fragments = np.arange(100.0, 110.0)
        intensities = [100.0 - i for i in range(10)]
        for i in range(30):
            shift_da = math.sin(i * 0.73) * 0.002
            precursor = 500.0 + math.sin(i * 0.61) * 0.001
            spectra.append(spectrum(
                i * 2.0,
                intensities,
                2,
                charge=2,
                precursor_mz=precursor,
                mz=(base_fragments + shift_da).tolist(),
                representation="centroid",
            ))
        result = Analyzer(MemoryReader(spectra)).analyze("repeat.mzML", collectors=[collector])
        precursor = next(
            a for a in result.annotations if a.field == "estimated_precursor_mass_error_ppm"
        )
        fragment = next(
            a for a in result.annotations if a.field == "estimated_fragment_mass_error_da"
        )
        self.assertEqual(precursor.kind, EvidenceKind.INFERRED)
        self.assertEqual(fragment.kind, EvidenceKind.INFERRED)
        self.assertGreater(precursor.value["single_measurement_sigma"], 0)
        self.assertGreater(fragment.value["single_measurement_sigma"], 0)
        self.assertIsNone(precursor.sdrf_value)
        self.assertIsNone(fragment.sdrf_value)
        self.assertGreaterEqual(precursor.support, 10)
        self.assertGreaterEqual(fragment.support, 80)

    def test_unknown_peak_type_estimated_centroid_supports_fragment_precision(self):
        collector = RepeatSpectrumMassErrorCollector(
            min_spectrum_pairs=10,
            min_precursor_clusters=1,
            min_fragment_pairs=80,
        )
        spectra = []
        base_fragments = np.arange(100.0, 110.0)
        intensities = [100.0 - i for i in range(10)]
        for i in range(30):
            item = spectrum(
                i * 2.0,
                intensities,
                2,
                charge=2,
                precursor_mz=500.0 + math.sin(i * 0.61) * 0.001,
                mz=(base_fragments + math.sin(i * 0.73) * 0.002).tolist(),
                representation="unknown",
            )
            spectra.append(replace(item, estimated_representation="centroid"))
        result = Analyzer(MemoryReader(spectra)).analyze(
            "estimated-centroid.mzML",
            collectors=[collector],
        )
        fragment = next(
            a for a in result.annotations if a.field == "estimated_fragment_mass_error_da"
        )
        self.assertEqual(fragment.kind, EvidenceKind.INFERRED)
        self.assertGreater(fragment.value["single_measurement_sigma"], 0)

    def test_profile_fragment_precision_uses_ephemeral_peak_centers(self):
        collector = RepeatSpectrumMassErrorCollector(
            min_spectrum_pairs=10,
            min_precursor_clusters=1,
            min_fragment_pairs=80,
        )
        spectra = []
        base_fragments = np.arange(100.0, 110.0)
        offsets = np.linspace(-0.03, 0.03, 7)
        for i in range(30):
            shift_da = math.sin(i * 0.73) * 0.002
            precursor = 500.0 + math.sin(i * 0.61) * 0.001
            mz = []
            intensities = []
            for peak_index, center in enumerate(base_fragments):
                mz.extend(center + shift_da + offsets)
                intensities.extend(
                    (100.0 - peak_index) * np.exp(-0.5 * (offsets / 0.012) ** 2)
                )
            spectra.append(spectrum(
                i * 2.0,
                list(intensities),
                2,
                charge=2,
                precursor_mz=precursor,
                mz=list(mz),
                representation="profile",
            ))
        result = Analyzer(MemoryReader(spectra)).analyze(
            "profile-peaks.mzML",
            collectors=[collector],
        )
        fragment = next(
            a for a in result.annotations if a.field == "estimated_fragment_mass_error_da"
        )
        diagnostic = next(
            a for a in result.annotations if a.field == "mass_error_estimator_diagnostics"
        )
        self.assertEqual(fragment.kind, EvidenceKind.INFERRED)
        self.assertGreater(fragment.value["single_measurement_sigma"], 0)
        self.assertGreaterEqual(fragment.support, 80)
        self.assertEqual(diagnostic.value["fragment_profile_spectra"], 30)
        self.assertEqual(diagnostic.value["fragment_centroid_spectra"], 0)
        self.assertGreaterEqual(diagnostic.value["fragment_profile_centroids"], 300)
        self.assertEqual(diagnostic.value["excluded_profile_peak_pick_failure"], 0)

    def test_profile_fragment_precision_abstains_without_resolvable_peaks(self):
        collector = RepeatSpectrumMassErrorCollector(
            min_spectrum_pairs=10,
            min_fragment_pairs=10,
        )
        spectra = []
        for cycle in range(4):
            for target_index in range(20):
                base = 500.0 + target_index * 2.0
                precursor = base + math.sin(cycle * 0.71 + target_index * 0.13) * 0.001
                spectra.append(spectrum(
                    cycle * 30 + target_index,
                    [10.0] * 10,
                    2,
                    charge=2,
                    precursor_mz=precursor,
                    mz=[100.0 + j for j in range(10)],
                    representation="profile",
                ))
        result = Analyzer(MemoryReader(spectra)).analyze("profile-flat.mzML", collectors=[collector])
        precursor = next(
            a for a in result.annotations if a.field == "estimated_precursor_mass_error_ppm"
        )
        fragment = next(
            a for a in result.annotations if a.field == "estimated_fragment_mass_error_da"
        )
        diagnostic = next(
            a for a in result.annotations if a.field == "mass_error_estimator_diagnostics"
        )
        self.assertEqual(precursor.kind, EvidenceKind.INFERRED)
        self.assertGreater(precursor.value["single_measurement_sigma"], 0)
        self.assertEqual(fragment.kind, EvidenceKind.UNAVAILABLE)
        self.assertIsNone(fragment.value)
        self.assertEqual(diagnostic.value["precursor_eligible_ms2"], 80)
        self.assertGreaterEqual(diagnostic.value["precursor_paired_spectra"], 10)
        self.assertEqual(diagnostic.value["fragment_eligible_ms2"], 0)
        self.assertEqual(diagnostic.value["fragment_profile_spectra"], 80)
        self.assertEqual(diagnostic.value["excluded_profile_peak_pick_failure"], 80)

    def test_precursor_search_tolerance_suggestion_requires_strong_diverse_support(self):
        collector = RepeatSpectrumMassErrorCollector(
            min_spectrum_pairs=10,
            min_precursor_clusters=2,
            min_tolerance_pairs=20,
            min_tolerance_clusters=5,
            tolerance_sigma_multiplier=6.0,
        )
        spectra = []
        for cycle in range(8):
            for target_index in range(12):
                base = 450.0 + target_index * 5.0
                precursor = base + math.sin(cycle * 0.83 + target_index * 0.19) * 0.001
                spectra.append(spectrum(
                    cycle * 30 + target_index,
                    [10.0] * 10,
                    2,
                    charge=2,
                    precursor_mz=precursor,
                    mz=[100.0 + j for j in range(10)],
                    representation="profile",
                ))
        result = Analyzer(MemoryReader(spectra)).analyze("tolerance.mzML", collectors=[collector])
        suggestion = next(
            a for a in result.annotations
            if a.field == "suggested_precursor_search_tolerance_ppm"
        )
        precision = next(
            a for a in result.annotations if a.field == "estimated_precursor_mass_error_ppm"
        )
        self.assertEqual(suggestion.kind, EvidenceKind.INFERRED)
        self.assertIsNone(suggestion.sdrf_value)
        self.assertAlmostEqual(
            suggestion.value["suggested_tolerance"],
            6.0 * precision.value["single_measurement_sigma"],
        )
        self.assertIn(suggestion.value["confidence"], {"moderate", "high"})

    def test_fragment_search_tolerance_suggestion_uses_ppm_for_high_resolution(self):
        collector = RepeatSpectrumMassErrorCollector(
            min_fragment_pairs=10,
            min_fragment_tolerance_pairs=20,
            min_fragment_tolerance_spectra=5,
        )
        collector.fragment_errors_da = [
            math.sin(i * 0.47) * 0.001
            for i in range(200)
        ]
        collector.fragment_errors_ppm = [
            math.sin(i * 0.47) * 2.0
            for i in range(200)
        ]
        collector.fragment_paired_spectra = 25
        collector.fragment_eligible_ms2 = 200
        annotations = collector.annotations()
        ppm = next(
            a for a in annotations
            if a.field == "suggested_fragment_search_tolerance_ppm"
        )
        da = next(
            a for a in annotations
            if a.field == "suggested_fragment_search_tolerance_da"
        )
        precision = next(
            a for a in annotations
            if a.field == "estimated_fragment_mass_error_ppm"
        )
        self.assertEqual(ppm.kind, EvidenceKind.INFERRED)
        self.assertEqual(ppm.value["resolution_regime"], "high-resolution")
        self.assertAlmostEqual(
            ppm.value["suggested_tolerance"],
            6.0 * precision.value["single_measurement_sigma"],
        )
        self.assertEqual(da.kind, EvidenceKind.UNAVAILABLE)
        self.assertIsNone(da.value)

    def test_fragment_tolerance_reports_robust_inlier_diagnostics(self):
        collector = RepeatSpectrumMassErrorCollector(
            min_fragment_pairs=10,
            min_fragment_tolerance_pairs=20,
            min_fragment_tolerance_spectra=5,
        )
        core_da = [math.sin(i * 0.47) * 0.001 for i in range(200)]
        core_ppm = [math.sin(i * 0.47) * 2.0 for i in range(200)]
        collector.fragment_errors_da = core_da + [0.05, -0.05] * 10
        collector.fragment_errors_ppm = core_ppm + [100.0, -100.0] * 10
        collector.fragment_paired_spectra = 25
        collector.fragment_eligible_ms2 = 220

        annotations = collector.annotations()
        precision = next(
            a for a in annotations
            if a.field == "estimated_fragment_mass_error_ppm"
        )
        suggestion = next(
            a for a in annotations
            if a.field == "suggested_fragment_search_tolerance_ppm"
        )

        self.assertEqual(suggestion.kind, EvidenceKind.INFERRED)
        self.assertGreater(precision.value["robust_inlier_fraction"], 0.8)
        self.assertLess(precision.value["robust_inlier_fraction"], 1.0)
        self.assertEqual(
            precision.value["robust_inlier_count"]
            + precision.value["robust_outlier_count"],
            220,
        )
        self.assertAlmostEqual(
            suggestion.value["robust_inlier_fraction"],
            precision.value["robust_inlier_fraction"],
        )
        self.assertEqual(
            suggestion.value["robust_inlier_count"],
            precision.value["robust_inlier_count"],
        )
        self.assertEqual(
            suggestion.value["robust_outlier_count"],
            precision.value["robust_outlier_count"],
        )

    def test_fragment_search_tolerance_suggestion_uses_da_for_low_resolution(self):
        collector = RepeatSpectrumMassErrorCollector(
            min_fragment_pairs=10,
            min_fragment_tolerance_pairs=20,
            min_fragment_tolerance_spectra=5,
        )
        collector.fragment_errors_da = [
            math.sin(i * 0.47) * 0.08
            for i in range(200)
        ]
        collector.fragment_errors_ppm = [
            math.sin(i * 0.47) * 120.0
            for i in range(200)
        ]
        collector.fragment_paired_spectra = 25
        collector.fragment_eligible_ms2 = 200
        annotations = collector.annotations()
        ppm = next(
            a for a in annotations
            if a.field == "suggested_fragment_search_tolerance_ppm"
        )
        da = next(
            a for a in annotations
            if a.field == "suggested_fragment_search_tolerance_da"
        )
        precision = next(
            a for a in annotations
            if a.field == "estimated_fragment_mass_error_da"
        )
        self.assertEqual(da.kind, EvidenceKind.INFERRED)
        self.assertEqual(da.value["resolution_regime"], "low-resolution")
        self.assertAlmostEqual(
            da.value["suggested_tolerance"],
            6.0 * precision.value["single_measurement_sigma"],
        )
        self.assertEqual(ppm.kind, EvidenceKind.UNAVAILABLE)
        self.assertIsNone(ppm.value)

    def test_fragment_search_tolerance_suggestion_abstains_on_ambiguous_regime(self):
        collector = RepeatSpectrumMassErrorCollector(
            min_fragment_pairs=10,
            min_fragment_tolerance_pairs=20,
            min_fragment_tolerance_spectra=5,
        )
        collector.fragment_errors_da = [
            math.sin(i * 0.47) * 0.02
            for i in range(200)
        ]
        collector.fragment_errors_ppm = [
            math.sin(i * 0.47) * 12.0
            for i in range(200)
        ]
        collector.fragment_paired_spectra = 25
        collector.fragment_eligible_ms2 = 200
        annotations = collector.annotations()
        for field_name in (
            "suggested_fragment_search_tolerance_ppm",
            "suggested_fragment_search_tolerance_da",
        ):
            suggestion = next(a for a in annotations if a.field == field_name)
            self.assertEqual(suggestion.kind, EvidenceKind.UNAVAILABLE)
            self.assertIsNone(suggestion.value)

    def test_precursor_search_tolerance_suggestion_abstains_on_small_target_grid(self):
        collector = RepeatSpectrumMassErrorCollector(
            min_spectrum_pairs=10,
            min_precursor_clusters=2,
            min_tolerance_pairs=20,
            min_tolerance_clusters=10,
        )
        spectra = []
        for cycle in range(20):
            for target_index in range(4):
                base = 500.0 + target_index * 10.0
                precursor = base + math.sin(cycle * 0.67 + target_index) * 0.001
                spectra.append(spectrum(
                    cycle * 10 + target_index,
                    [10.0] * 10,
                    2,
                    charge=2,
                    precursor_mz=precursor,
                    mz=[100.0 + j for j in range(10)],
                    representation="profile",
                ))
        result = Analyzer(MemoryReader(spectra)).analyze("fixed-grid.mzML", collectors=[collector])
        precision = next(
            a for a in result.annotations if a.field == "estimated_precursor_mass_error_ppm"
        )
        suggestion = next(
            a for a in result.annotations
            if a.field == "suggested_precursor_search_tolerance_ppm"
        )
        self.assertEqual(precision.kind, EvidenceKind.INFERRED)
        self.assertEqual(suggestion.kind, EvidenceKind.UNAVAILABLE)
        self.assertIsNone(suggestion.value)

    def test_profile_and_unknown_spectra_do_not_become_reporter_evidence(self):
        collector = DiagnosticIonCollector()
        collector.consume_spectrum(
            spectrum(
                0,
                [10, 10, 10],
                2,
                mz=[126.127726, 127.131081, 128.134436],
                representation="profile",
            ),
        )
        self.assertEqual(collector.total, 0)
        self.assertEqual(collector.skipped_profile, 1)

    def test_acquisition_heuristic_is_not_observed_fact(self):
        result = analyze([spectrum(i, [1], 2, charge=2, width=20) for i in range(100)])
        evidence = next(a for a in result.annotations if a.field == "acquisition_method")
        self.assertEqual(evidence.kind, EvidenceKind.INFERRED)
        self.assertIn("PRIDE:0000450", evidence.sdrf_value)

    def test_small_fixed_target_narrow_run_abstains_for_prm_or_narrow_dia_ambiguity(self):
        targets = [400.0, 500.0, 600.0, 700.0]
        spectra = []
        rt = 0
        for _ in range(25):
            spectra.append(spectrum(rt, [1], 1))
            rt += 1
            for target in targets:
                spectra.append(spectrum(rt, [1], 2, charge=2, width=2.0, precursor_mz=target))
                rt += 1
        result = analyze(spectra)
        evidence = next(a for a in result.annotations if a.field == "acquisition_method")
        self.assertIsNone(evidence.value)

    def test_narrow_dynamic_many_target_run_is_inferred_dda(self):
        spectra = []
        rt = 0
        precursor = 400.0
        for _ in range(60):
            spectra.append(spectrum(rt, [1], 1))
            rt += 1
            for _ in range(10):
                spectra.append(spectrum(rt, [1], 2, charge=2, width=2.0, precursor_mz=precursor))
                precursor += 0.2
                rt += 1
        result = analyze(spectra)
        evidence = next(a for a in result.annotations if a.field == "acquisition_method")
        self.assertEqual(evidence.value, "Data-dependent acquisition")
        self.assertIn("PRIDE:0000627", evidence.sdrf_value)

    def test_stable_intermediate_width_cycles_are_inferred_dia(self):
        targets = [400.0 + 10 * i for i in range(8)]
        spectra = []
        rt = 0
        for _ in range(25):
            spectra.append(spectrum(rt, [1], 1))
            rt += 1
            for target in targets:
                spectra.append(spectrum(rt, [1], 2, charge=2, width=10.0, precursor_mz=target))
                rt += 1
        result = analyze(spectra)
        evidence = next(a for a in result.annotations if a.field == "acquisition_method")
        self.assertEqual(evidence.value, "Data-independent acquisition")
        self.assertEqual(evidence.support, 200)
        self.assertEqual(evidence.total, 200)
        self.assertIn("PRIDE:0000450", evidence.sdrf_value)

    def test_changing_intermediate_target_grid_abstains(self):
        spectra = []
        rt = 0
        target = 400.0
        for _ in range(25):
            spectra.append(spectrum(rt, [1], 1))
            rt += 1
            for _ in range(8):
                spectra.append(spectrum(rt, [1], 2, charge=2, width=10.0, precursor_mz=target))
                target += 1.0
                rt += 1
        result = analyze(spectra)
        evidence = next(a for a in result.annotations if a.field == "acquisition_method")
        self.assertIsNone(evidence.value)

    def test_sparse_or_mixed_widths_abstain(self):
        for widths in ([20] * 10, [2] * 50 + [20] * 50):
            result = analyze([spectrum(i, [1], 2, charge=2, width=w) for i, w in enumerate(widths)])
            evidence = next(a for a in result.annotations if a.field == "acquisition_method")
            self.assertIsNone(evidence.value)

    def test_tolerance_is_not_guessed_from_instrument_or_width(self):
        result = analyze()
        evidence = next(a for a in result.annotations if a.field == "precursor_mass_tolerance")
        self.assertEqual(evidence.kind, EvidenceKind.UNAVAILABLE)
        self.assertIsNone(evidence.sdrf_value)


class SDRFTests(unittest.TestCase):
    def test_not_applicable_is_preserved_as_an_existing_assertion(self):
        document = SDRFDocument(
            ["comment[data file]", "comment[instrument]"],
            [["run.mzML", "not applicable"]],
        )
        self.assertEqual(document.annotate([analyze()])[0].status, "conflict")
        self.assertEqual(document.rows[0][1], "not applicable")

    def test_duplicate_modification_columns_and_repeated_sample_rows_survive(self):
        columns = ["source name", "comment[data file]", "comment[modification parameters]", "comment[modification parameters]"]
        document = SDRFDocument(columns.copy(), [["s1", "run.raw", "Oxidation", "Carbamidomethyl"],
                                                 ["s2", "run.raw", "", ""]], ["# meta\n"])
        result = analyze(metadata=RunMetadata(instruments=[FUSION], source_files=["run.raw"]))
        changes = document.annotate([result])
        self.assertEqual(document.columns[:4], columns)
        self.assertEqual(document.rows[0][:4], ["s1", "run.raw", "Oxidation", "Carbamidomethyl"])
        self.assertEqual(sum(c.status == "filled" for c in changes), 2)
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "out.tsv"
            document.write(path)
            again = SDRFDocument.read(path)
            self.assertEqual(again.columns, document.columns)
            self.assertEqual(again.rows, document.rows)
            self.assertEqual(again.preamble, ["# meta\n"])

    def test_conflict_preserved_until_explicit_overwrite(self):
        document = SDRFDocument(
            ["comment[data file]", "comment[instrument]"],
            [["run.mzML", "existing model"]],
        )
        changes = document.annotate([analyze()])
        self.assertEqual(changes[0].status, "conflict")
        self.assertEqual(document.rows[0][1], "existing model")
        document.annotate([analyze()], overwrite=True)
        self.assertIn("MS:1002416", document.rows[0][1])

    def test_cv_accession_match_preserves_original_label(self):
        document = SDRFDocument(
            ["comment[data file]", "comment[instrument]"],
            [["run.mzML", "AC=MS:1002416;NT=old spelling"]],
        )
        self.assertEqual(document.annotate([analyze()])[0].status, "match")
        self.assertEqual(document.rows[0][1], "AC=MS:1002416;NT=old spelling")

    def test_no_implicit_stem_matching(self):
        document = SDRFDocument(["comment[data file]"], [["run.raw"]])
        self.assertEqual(document.annotate([analyze()])[0].status, "unmatched")
        changes = document.annotate([analyze()], aliases={"run.raw": "run.mzML"})
        self.assertEqual(changes[0].status, "filled")

    def test_ambiguous_source_alias_is_rejected_before_mutation(self):
        metadata = RunMetadata(instruments=[FUSION], source_files=["run.raw"])
        document = SDRFDocument(["comment[data file]"], [["run.raw"]])
        with self.assertRaisesRegex(ValueError, "Ambiguous"):
            document.annotate(
                [analyze(
                    name="a.mzML",
                    metadata=metadata,
                ), analyze(
                    name="b.mzML",
                    metadata=metadata,
                )],
            )
        self.assertEqual(len(document.columns), 1)

    def test_inferred_annotation_requires_explicit_option(self):
        result = analyze()
        result.annotations = [Annotation(
            "mode",
            "DIA",
            EvidenceKind.INFERRED,
            "test",
            sdrf_column="comment[proteomics data acquisition method]",
            sdrf_value="NT=Data-independent acquisition;AC=PRIDE:0000450",
        )]
        document = SDRFDocument(["comment[data file]"], [["run.mzML"]])
        self.assertEqual(document.annotate([result])[0].status, "suggestion")
        self.assertEqual(len(document.columns), 1)
        self.assertEqual(document.annotate([result], include_inferred=True)[0].status, "filled")

    def test_duplicate_target_columns_are_reported_without_change(self):
        document = SDRFDocument(
            ["comment[data file]", "comment[instrument]", "comment[instrument]"],
            [["run.mzML", "", ""]],
        )
        self.assertEqual(document.annotate([analyze()])[0].status, "ambiguous_column")
        self.assertEqual(document.rows[0][1:], ["", ""])

    def test_invalid_tsv_width_or_duplicate_data_file_column_rejected(self):
        with self.assertRaises(ValueError):
            SDRFDocument(["comment[data file]"], [["a", "extra"]])
        with self.assertRaises(ValueError):
            SDRFDocument(["comment[data file]", "comment[data file]"], [])


class SerializationTests(unittest.TestCase):
    def test_json_nonfinite_values_are_null_recursively(self):
        value = json_safe({"a": np.array([1, math.nan]), "b": {"x": math.inf}})
        self.assertEqual(value, {"a": [1.0, None], "b": {"x": None}})
        json.dumps(value, allow_nan=False)

    def test_mzqc_one_run_has_complete_provenance_and_no_duplicate_accessions(self):
        result = analyze()
        data = MzQCWriter().build(result, Path("run.obo"))["mzQC"]
        self.assertEqual(len(data["runQualities"]), 1)
        run = data["runQualities"][0]
        self.assertTrue(run["metadata"]["inputFiles"][0]["location"].startswith("file://"))
        self.assertEqual(run["metadata"]["analysisSoftware"][0]["name"], "prideqc")
        accessions = [m["accession"] for m in run["qualityMetrics"]]
        self.assertEqual(len(accessions), len(set(accessions)))
        self.assertNotIn("setQualities", data)
        self.assertTrue(data["creationDate"].endswith("+00:00"))

    @unittest.skipUnless(importlib.util.find_spec("jsonschema"), "jsonschema not installed")
    def test_official_mzqc_schema_with_format_checks(self):
        import jsonschema

        schema = json.loads((Path(__file__).parent / "data" / "mzqc_schema.json").read_text())
        for items in ([], [spectrum(0, [1])], [spectrum(0, [1], 2, charge=2)]):
            result = Analyzer(MemoryReader(items)).analyze("run.mzML")
            jsonschema.Draft7Validator(schema, format_checker=jsonschema.FormatChecker()).validate(
                MzQCWriter().build(result, Path("local.obo")))

    def test_cv_shape_stable_when_unavailable(self):
        self.assertEqual(definition(Metric("MzRange_MS2", None)).shape, "tuple")
        self.assertEqual(definition(Metric("MzRange_MS2", [400, 600])).shape, "tuple")

    def test_unknown_metric_cannot_be_mislabeled(self):
        with self.assertRaisesRegex(ValueError, "Unregistered"):
            definition(Metric("made_up", 1))

    def test_atomic_write_failure_leaves_previous_file(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "result.json"
            path.write_text("old")
            with self.assertRaises(RuntimeError), atomic_text(path) as handle:
                handle.write("partial")
                raise RuntimeError("failure")
            self.assertEqual(path.read_text(), "old")
            self.assertEqual(len(list(Path(folder).iterdir())), 1)


class ReaderHeaderTests(unittest.TestCase):
    def test_preserves_cv_and_reference_groups_without_reading_payload(self):
        xml = '''<mzML xmlns="http://psi.hupo.org/ms/mzml">
        <referenceableParamGroupList><referenceableParamGroup id="inst">
        <cvParam accession="MS:1002416" name="Orbitrap Fusion" value=""/>
        </referenceableParamGroup></referenceableParamGroupList>
        <fileDescription><sourceFileList><sourceFile id="x" name="run.raw" location="file:///data"/></sourceFileList></fileDescription>
        <instrumentConfigurationList><instrumentConfiguration id="IC1">
        <referenceableParamGroupRef ref="inst"/>
        <cvParam accession="MS:1000529" name="instrument serial number" value="ABC"/>
        <componentList><analyzer><cvParam accession="MS:1000484" name="orbitrap"/></analyzer></componentList>
        </instrumentConfiguration></instrumentConfigurationList>
        <run startTimeStamp="2026-01-01T00:00:00Z"><spectrumList count="1">'''
        # Malformed content is far beyond the header parser's read-ahead buffer.
        xml += " " * 65536 + "<invalid payload"
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "header.mzML"
            path.write_text(xml)
            metadata = read_header(path)
        self.assertEqual(metadata.instruments, [FUSION])
        self.assertEqual(metadata.serial_numbers, ["ABC"])
        self.assertEqual(metadata.source_files, ["run.raw"])
        self.assertEqual(metadata.analyzers[0].accession, "MS:1000484")


class WorkflowTests(unittest.TestCase):
    def test_progress_reports_completed_files_and_can_be_disabled(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            source = root / "run.mzML"
            source.touch()
            result = analyze(name=str(source))
            with patch("prideqc.pipeline._analyze_file", return_value=FileOutcome(source, result)):
                stream = io.StringIO()
                with redirect_stderr(stream):
                    Workflow(WorkflowOptions(progress=True)).run([source], root / "progress")
                self.assertIn("Analyzing files: 1/1", stream.getvalue())
                self.assertIn("Analyzing files: 0/1", stream.getvalue())
                self.assertNotIn("\r", stream.getvalue())
            with patch("prideqc.pipeline._analyze_file", return_value=FileOutcome(source, result)):
                stream = io.StringIO()
                with redirect_stderr(stream):
                    Workflow(WorkflowOptions(progress=False)).run([source], root / "quiet")
                self.assertEqual(stream.getvalue(), "")

    def test_partial_failure_writes_manifest_and_successful_mzqc(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            good, bad = root / "good.mzML", root / "bad.mzML"
            good.touch()
            bad.touch()
            result = analyze(name=str(good))
            with patch(
                "prideqc.pipeline._analyze_file",
                side_effect=[FileOutcome(good, result), FileOutcome(bad, error="broken file")],
            ):
                manifest = Workflow(WorkflowOptions(continue_on_error=True)).run(
                    [good, bad],
                    root / "out",
                )
            self.assertFalse(manifest["success"])
            self.assertTrue((root / "out/good.mzML.mzQC").exists())
            self.assertFalse((root / "out/bad.mzML.mzQC").exists())
            self.assertEqual(
                json.loads((root / "out/manifest.json").read_text())["files"][1]["status"],
                "failed",
            )

    def test_output_collision_and_existing_results_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            source = root / "run.mzML"
            source.touch()
            with self.assertRaises(ValueError):
                Workflow().run([source, source], root / "out")
            (root / "out").mkdir()
            (root / "out" / "existing").touch()
            with self.assertRaises(FileExistsError):
                Workflow().run([source], root / "out")

    def test_overwrite_clears_existing_results_directory(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            source = root / "run.mzML"
            source.touch()
            output = root / "out"
            output.mkdir()
            (output / "old-result.txt").write_text("stale")
            result = analyze(name=str(source))
            with patch("prideqc.pipeline._analyze_file", return_value=FileOutcome(source, result)):
                manifest = Workflow(WorkflowOptions(overwrite=True)).run([source], output)
            self.assertTrue(manifest["success"])
            self.assertFalse((output / "old-result.txt").exists())
            self.assertTrue((output / "run.mzML.mzQC").exists())

    def test_fail_fast_records_unprocessed_files(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            first, second = root / "a.mzML", root / "b.mzML"
            first.touch()
            second.touch()
            with patch(
                "prideqc.pipeline._analyze_file",
                return_value=FileOutcome(first, error="bad"),
            ) as worker:
                manifest = Workflow().run([first, second], root / "out")
            self.assertEqual(worker.call_count, 1)
            self.assertEqual(manifest["unprocessed"], [str(second)])

    def test_cli_failed_run_returns_nonzero(self):
        with patch(
            "prideqc.pipeline.Workflow.run",
            return_value={"files": [], "errors": [], "success": False},
        ):
            self.assertEqual(main(["analyze", "a.mzML", "-o", "out"]), 1)

    def test_help_does_not_import_runtime_dependencies(self):
        code = (
            "from prideqc.cli import parser; import sys; parser(); "
            "assert not {'pyopenms', 'numpy', 'pridepy', 'sdrf_pipelines'} & sys.modules.keys()"
        )
        subprocess.run([sys.executable, "-c", code], check=True)

    def test_converter_arguments_preserve_paths_and_never_use_shell(self):
        converter = ExternalConverter("msconvert")
        source = Path("/data/run; echo secret.wiff")
        command = converter.command(source, Path("/data/output folder"))
        self.assertEqual(command[1], str(source))
        self.assertEqual(command[-1], "/data/output folder")

    def test_input_validation(self):
        with self.assertRaises(ValueError):
            WorkflowOptions(workers=0)
        with self.assertRaises(ValueError):
            DiagnosticIonCollector(ppm=math.nan)
