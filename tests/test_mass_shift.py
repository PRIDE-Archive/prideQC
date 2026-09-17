"""Identification-free recurrent mass-shift scout contracts."""

from __future__ import annotations

import json
import math
import tempfile
import unittest
from dataclasses import replace
from pathlib import Path
from unittest.mock import patch

import numpy as np

from prideqc.cli import parser
from prideqc.mass_error import RepeatSpectrumMassErrorCollector
from prideqc.mass_shift import MassShiftCollector, ModificationRecord, load_openms_modifications
from prideqc.models import AnalysisResult, EvidenceKind, RunMetadata
from prideqc.mzqc import MzQCWriter
from prideqc.pipeline import FileOutcome, Workflow, WorkflowOptions
from tests.helpers import spectrum

PHOSPHO = ModificationRecord(
    "UNIMOD:21",
    "Phospho",
    79.966331,
    ("S", "T", "Y"),
    ("Anywhere",),
    "Post-translational",
)


def related_pair(
    family: int,
    delta_da: float,
    *,
    rt: float | None = None,
    representation: str = "centroid",
):
    """Return a synthetic related MS2 pair with unchanged and shifted fragments."""

    # Each family gets a different coarse fingerprint so unrelated families are
    # not candidate neighbours. Fifteen peaks stay unchanged and fifteen move by
    # the precursor delta, which is enough for the v22 relatedness screen.
    base = 100.0 + family * 37.0 + np.arange(30, dtype=float) * 3.11
    shifted = base.copy()
    shifted[15:] += delta_da
    intensities = (200.0 - np.arange(30, dtype=float)).tolist()
    precursor_mz = 500.0 + family * 12.5
    start_rt = float(family * 20 if rt is None else rt)
    left = spectrum(
        start_rt,
        intensities,
        2,
        charge=2,
        precursor_mz=precursor_mz,
        mz=base.tolist(),
        representation=representation,
    )
    right = spectrum(
        start_rt + 2.0,
        intensities,
        2,
        charge=2,
        precursor_mz=precursor_mz + delta_da / 2.0,
        mz=shifted.tolist(),
        representation=representation,
    )
    return left, right


class CliMassShiftTests(unittest.TestCase):
    def test_mass_shift_scout_is_explicitly_opt_in(self):
        arguments = parser().parse_args(
            ["analyze", "run.mzML", "-o", "out", "--estimate-mass-shifts"]
        )
        self.assertTrue(arguments.estimate_mass_shifts)
        self.assertFalse(arguments.estimate_mass_error)


class MassShiftCollectorTests(unittest.TestCase):
    def test_recurrent_related_spectra_find_putative_unimod_candidate(self):
        collector = MassShiftCollector(
            modifications=(PHOSPHO,),
            minimum_cluster_pairs=6,
            minimum_cluster_unique_spectra=10,
        )
        for family in range(10):
            delta = 79.966331 + math.sin(family * 0.73) * 0.0015
            for item in related_pair(family, delta):
                collector.consume_spectrum(item)

        annotations = collector.annotations()
        shifts = next(
            item for item in annotations if item.field == "putative_modification_mass_shifts"
        )
        diagnostics = next(
            item for item in annotations if item.field == "mass_shift_scout_diagnostics"
        )

        self.assertEqual(shifts.kind, EvidenceKind.INFERRED)
        self.assertEqual(len(shifts.value), 1)
        cluster = shifts.value[0]
        self.assertAlmostEqual(cluster["delta_mass_da"], 79.966331, places=2)
        self.assertGreaterEqual(cluster["pair_support"], 6)
        self.assertGreaterEqual(cluster["unique_spectrum_support"], 10)
        self.assertEqual(cluster["classification"], "putative-ptm")
        self.assertEqual(cluster["unimod_candidates"][0]["unimod_accession"], "UNIMOD:21")
        self.assertEqual(cluster["unimod_candidates"][0]["name"], "Phospho")
        self.assertIn("localization", shifts.detail)
        self.assertGreater(diagnostics.value["scored_candidate_pairs"], 0)
        self.assertGreater(diagnostics.value["accepted_related_pairs"], 0)

    def test_isotope_shift_is_classified_before_unimod(self):
        isotope = 1.00335483507
        fake = ModificationRecord(
            "UNIMOD:99999",
            "Fake isotope-mass modification",
            isotope,
            (),
            (),
            "Post-translational",
        )
        collector = MassShiftCollector(
            modifications=(fake,),
            minimum_cluster_pairs=4,
            minimum_cluster_unique_spectra=6,
        )
        for family in range(8):
            for item in related_pair(family, isotope):
                collector.consume_spectrum(item)
        cluster = collector.annotations()[0].value[0]
        self.assertEqual(cluster["classification"], "isotope-like")
        self.assertIn("isotope-selection", cluster["artifact_candidate"]["name"])
        self.assertEqual(cluster["unimod_candidates"], [])

    def test_unknown_recurrent_shift_remains_visible(self):
        collector = MassShiftCollector(
            modifications=(),
            minimum_cluster_pairs=4,
            minimum_cluster_unique_spectra=6,
        )
        for family in range(8):
            for item in related_pair(family, 23.4567):
                collector.consume_spectrum(item)
        cluster = collector.annotations()[0].value[0]
        self.assertEqual(cluster["classification"], "unknown")
        self.assertEqual(cluster["unimod_candidates"], [])

    def test_unrelated_spectra_with_same_precursor_delta_do_not_create_cluster(self):
        collector = MassShiftCollector(
            modifications=(PHOSPHO,),
            minimum_cluster_pairs=3,
            minimum_cluster_unique_spectra=4,
        )
        intensities = (200.0 - np.arange(30, dtype=float)).tolist()
        for family in range(8):
            precursor_mz = 500.0 + family * 10.0
            left_mz = (100.0 + family * 50.0 + np.arange(30) * 2.37).tolist()
            right_mz = (1700.0 + family * 50.0 + np.arange(30) * 2.61).tolist()
            collector.consume_spectrum(
                spectrum(
                    family * 20.0, intensities, 2, charge=2,
                    precursor_mz=precursor_mz, mz=left_mz, representation="centroid",
                )
            )
            collector.consume_spectrum(
                spectrum(
                    family * 20.0 + 2.0, intensities, 2, charge=2,
                    precursor_mz=precursor_mz + 79.966331 / 2.0,
                    mz=right_mz, representation="centroid",
                )
            )
        shifts = collector.annotations()[0]
        self.assertIsNone(shifts.value)
        self.assertEqual(shifts.kind, EvidenceKind.UNAVAILABLE)

    def test_cluster_result_is_stable_when_pair_order_is_shuffled(self):
        items = []
        for family in range(10):
            items.extend(related_pair(family, 79.966331 + math.sin(family) * 0.001))

        first = MassShiftCollector(
            modifications=(PHOSPHO,),
            minimum_cluster_pairs=6,
            minimum_cluster_unique_spectra=10,
        )
        second = MassShiftCollector(
            modifications=(PHOSPHO,),
            minimum_cluster_pairs=6,
            minimum_cluster_unique_spectra=10,
        )
        for item in items:
            first.consume_spectrum(item)
        order = np.random.default_rng(17).permutation(len(items))
        for index in order:
            second.consume_spectrum(items[int(index)])

        left = first.annotations()[0].value[0]
        right = second.annotations()[0].value[0]
        self.assertAlmostEqual(left["delta_mass_da"], right["delta_mass_da"], places=9)
        self.assertEqual(left["pair_support"], right["pair_support"])
        self.assertEqual(left["unique_spectrum_support"], right["unique_spectrum_support"])
        self.assertEqual(left["classification"], right["classification"])

    def test_missing_or_nonpositive_precursor_charge_abstains(self):
        collector = MassShiftCollector(modifications=())
        intensities = [100.0] * 20
        collector.consume_spectrum(
            spectrum(0.0, intensities, 2, charge=None, mz=list(np.arange(100.0, 120.0)))
        )
        collector.consume_spectrum(
            spectrum(2.0, intensities, 2, charge=0, mz=list(np.arange(100.0, 120.0)))
        )
        diagnostics = collector.annotations()[1].value
        self.assertEqual(diagnostics["eligible_centroid_ms2"], 0)
        self.assertEqual(diagnostics["excluded_missing_precursor"], 2)

    def test_profile_ms2_uses_ephemeral_openms_peak_picker(self):
        collector = MassShiftCollector(
            modifications=(PHOSPHO,),
            oms=FakeOpenMS,
            minimum_cluster_pairs=1,
            minimum_cluster_unique_spectra=2,
        )
        left, right = related_pair(0, 79.966331, representation="profile")
        collector.consume_spectrum(left)
        collector.consume_spectrum(right)
        shifts, diagnostics = collector.annotations()
        self.assertEqual(shifts.kind, EvidenceKind.INFERRED)
        self.assertEqual(shifts.value[0]["classification"], "putative-ptm")
        self.assertEqual(diagnostics.value["explicit_profile_ms2"], 2)
        self.assertEqual(diagnostics.value["profile_ms2_centroided"], 2)
        self.assertEqual(diagnostics.value["profile_peak_pick_failures"], 0)
        self.assertEqual(diagnostics.value["excluded_profile_or_unknown"], 0)

    def test_unknown_estimated_profile_uses_peak_picker(self):
        collector = MassShiftCollector(
            modifications=(PHOSPHO,),
            oms=FakeOpenMS,
            minimum_cluster_pairs=1,
            minimum_cluster_unique_spectra=2,
        )
        left, right = related_pair(0, 79.966331, representation="unknown")
        collector.consume_spectrum(replace(left, estimated_representation="profile"))
        collector.consume_spectrum(replace(right, estimated_representation="profile"))
        shifts, diagnostics = collector.annotations()
        self.assertEqual(shifts.kind, EvidenceKind.INFERRED)
        self.assertEqual(diagnostics.value["unknown_ms2"], 2)
        self.assertEqual(diagnostics.value["unknown_ms2_inferred_profile"], 2)
        self.assertEqual(diagnostics.value["profile_ms2_centroided"], 2)

    def test_unresolved_unknown_representation_abstains(self):
        collector = MassShiftCollector(modifications=())
        left, right = related_pair(0, 79.966331, representation="unknown")
        collector.consume_spectrum(left)
        collector.consume_spectrum(right)
        shifts, diagnostics = collector.annotations()
        self.assertIsNone(shifts.value)
        self.assertEqual(diagnostics.value["unknown_ms2"], 2)
        self.assertEqual(diagnostics.value["unresolved_spectrum_type"], 2)
        self.assertEqual(diagnostics.value["excluded_profile_or_unknown"], 2)

    def test_estimated_centroid_representation_is_eligible(self):
        collector = MassShiftCollector(
            modifications=(PHOSPHO,),
            minimum_cluster_pairs=1,
            minimum_cluster_unique_spectra=2,
        )
        left, right = related_pair(0, 79.966331, representation="unknown")
        collector.consume_spectrum(replace(left, estimated_representation="centroid"))
        collector.consume_spectrum(replace(right, estimated_representation="centroid"))
        shifts = collector.annotations()[0]
        self.assertEqual(shifts.kind, EvidenceKind.INFERRED)
        self.assertEqual(shifts.value[0]["unimod_candidates"][0]["unimod_accession"], "UNIMOD:21")

    def test_inverted_index_is_bounded(self):
        collector = MassShiftCollector(
            modifications=(),
            maximum_indexed_spectra=3,
            minimum_cluster_pairs=1,
            minimum_cluster_unique_spectra=2,
        )
        for family in range(10):
            left, _ = related_pair(family, 40.0)
            collector.consume_spectrum(left)
        diagnostics = collector.annotations()[1].value
        self.assertEqual(diagnostics["indexed_spectra_current"], 3)
        self.assertEqual(diagnostics["index_evictions"], 7)

    def test_v19_precursor_precision_is_consumed_read_only(self):
        precision = RepeatSpectrumMassErrorCollector()
        precision.precursor_errors_ppm.extend(
            [2.0 * math.sin(index * 0.31) for index in range(240)]
        )
        precision.precursor_paired_spectra = precision.min_tolerance_pairs
        precision.precursor_clusters_used = precision.min_tolerance_clusters
        before = list(precision.precursor_errors_ppm)

        collector = MassShiftCollector(
            precision_source=precision,
            modifications=(PHOSPHO,),
            minimum_cluster_pairs=4,
            minimum_cluster_unique_spectra=6,
        )
        for family in range(8):
            for item in related_pair(family, 79.966331):
                collector.consume_spectrum(item)
        diagnostics = collector.annotations()[1].value

        self.assertEqual(
            diagnostics["calibration"]["source"],
            "v19-repeat-precursor-precision-read-only",
        )
        self.assertGreater(
            diagnostics["calibration"]["unimod_tolerance_da"],
            0.0,
        )
        self.assertEqual(precision.precursor_errors_ppm, before)


    def test_unimod_window_does_not_inherit_widened_cluster_tolerance(self):
        precision = RepeatSpectrumMassErrorCollector()
        precision.precursor_errors_ppm.extend(
            [80.0 * math.sin(index * 0.31) for index in range(240)]
        )
        precision.precursor_paired_spectra = precision.min_tolerance_pairs
        precision.precursor_clusters_used = precision.min_tolerance_clusters
        near_but_outside = ModificationRecord(
            "UNIMOD:90001",
            "Near but outside annotation window",
            79.996331,
            (),
            (),
            "Post-translational",
        )
        collector = MassShiftCollector(
            precision_source=precision,
            modifications=(near_but_outside,),
            minimum_cluster_pairs=1,
            minimum_cluster_unique_spectra=2,
        )
        for item in related_pair(0, 79.966331):
            collector.consume_spectrum(item)
        shifts, diagnostics = collector.annotations()
        calibration = diagnostics.value["calibration"]
        self.assertGreater(calibration["cluster_tolerance_da"], 0.02)
        self.assertEqual(calibration["unimod_tolerance_da"], 0.02)
        self.assertEqual(calibration["unimod_tolerance_policy"], "independent-fixed-window")
        self.assertEqual(shifts.value[0]["unimod_candidates"], [])

    def test_decoys_and_amino_acid_substitutions_are_hidden_from_default_candidates(self):
        modifications = (
            PHOSPHO,
            ModificationRecord(
                "UNIMOD:99913",
                "Phosphorylation Decoy",
                79.966331,
                (),
                (),
                "Post-translational",
            ),
            ModificationRecord(
                "UNIMOD:99914",
                "Ser->Tyr substitution",
                79.966331,
                (),
                (),
                "AA substitution",
            ),
        )
        collector = MassShiftCollector(
            modifications=modifications,
            minimum_cluster_pairs=1,
            minimum_cluster_unique_spectra=2,
        )
        for item in related_pair(0, 79.966331):
            collector.consume_spectrum(item)
        cluster = collector.annotations()[0].value[0]
        self.assertEqual(
            [item["unimod_accession"] for item in cluster["unimod_candidates"]],
            ["UNIMOD:21"],
        )
        self.assertEqual(cluster["diagnostic_unimod_candidate_count"], 3)
        self.assertEqual(cluster["diagnostic_suppressed_unimod_candidates"], 2)

    def test_isotope_adjacent_satellite_is_suppressed_but_counted(self):
        collector = MassShiftCollector(
            modifications=(),
            minimum_cluster_pairs=1,
            minimum_cluster_unique_spectra=2,
        )
        for item in related_pair(0, 0.92):
            collector.consume_spectrum(item)
        shifts, diagnostics = collector.annotations()
        self.assertIsNone(shifts.value)
        self.assertEqual(diagnostics.value["raw_recurrent_clusters"], 1)
        self.assertEqual(diagnostics.value["suppressed_isotope_adjacent_clusters"], 1)
        self.assertEqual(diagnostics.value["reported_clusters"], 0)

    def test_unknown_reporting_is_bounded_without_discarding_raw_cluster_count(self):
        collector = MassShiftCollector(
            modifications=(),
            minimum_cluster_pairs=1,
            minimum_cluster_unique_spectra=2,
            maximum_reported_unknown_clusters=3,
        )
        for family in range(8):
            for item in related_pair(family, 20.0 + family * 7.0):
                collector.consume_spectrum(item)
        shifts, diagnostics = collector.annotations()
        self.assertEqual(len(shifts.value), 3)
        self.assertEqual(diagnostics.value["raw_recurrent_clusters"], 8)
        self.assertEqual(diagnostics.value["reported_clusters"], 3)
        self.assertEqual(diagnostics.value["suppressed_unknown_clusters"], 5)



class FakeModification:
    def __init__(self, accession, name, delta, origin, term, classification):
        self.accession = accession
        self.name = name
        self.delta = delta
        self.origin = origin
        self.term = term
        self.classification = classification

    def getUniModAccession(self):
        return self.accession

    def getDiffMonoMass(self):
        return self.delta

    def getFullName(self):
        return self.name

    def getId(self):
        return self.name

    def getOrigin(self):
        return self.origin

    def getTermSpecificity(self):
        return self.term

    def getTermSpecificityName(self, _value):
        return self.term

    def getSourceClassification(self):
        return self.classification

    def getSourceClassificationName(self, _value):
        return self.classification


class FakeModificationDB:
    def __init__(self):
        self.items = [
            FakeModification("UNIMOD:35", "Oxidation", 15.994915, "M", "Anywhere", "Post-translational"),
            FakeModification("UNIMOD:35", "Oxidation", 15.994915, "W", "Anywhere", "Post-translational"),
            FakeModification("MOD:00000", "Not UniMod", 12.3, "X", "Anywhere", "Other"),
        ]

    def getNumberOfModifications(self):
        return len(self.items)

    def getModification(self, index):
        return self.items[index]


class FakeMSSpectrum:
    def __init__(self):
        self._peaks = (np.array([], dtype=float), np.array([], dtype=float))

    def set_peaks(self, peaks):
        self._peaks = tuple(np.asarray(item, dtype=float).copy() for item in peaks)

    def get_peaks(self):
        return self._peaks


class FakePeakPickerHiRes:
    def pick(self, source, target):
        target.set_peaks(source.get_peaks())


class FakeOpenMS:
    class Constants:
        PROTON_MASS_U = 1.007276466621
        C13C12_MASSDIFF_U = 1.00335483507

    MSSpectrum = FakeMSSpectrum
    PeakPickerHiRes = FakePeakPickerHiRes

    @staticmethod
    def ModificationsDB():
        return FakeModificationDB()


class OpenMSCatalogTests(unittest.TestCase):
    def test_openms_unimod_records_are_deduplicated_and_keep_specificities(self):
        records = load_openms_modifications(FakeOpenMS)
        self.assertEqual(len(records), 1)
        record = records[0]
        self.assertEqual(record.accession, "UNIMOD:35")
        self.assertEqual(record.name, "Oxidation")
        self.assertEqual(record.origins, ("M", "W"))
        self.assertEqual(record.term_specificities, ("Anywhere",))
        self.assertEqual(record.source_classification, "Post-translational")


class MzQCMassShiftTests(unittest.TestCase):
    def test_structured_mass_shift_annotation_is_serialized_into_mzqc(self):
        collector = MassShiftCollector(
            modifications=(PHOSPHO,),
            minimum_cluster_pairs=1,
            minimum_cluster_unique_spectra=2,
        )
        for item in related_pair(0, 79.966331):
            collector.consume_spectrum(item)
        annotations = collector.annotations()
        result = AnalysisResult(
            Path("run.mzML"),
            RunMetadata(),
            [],
            annotations,
            [],
            "test-reader",
            0.1,
        )
        with tempfile.TemporaryDirectory() as tmp:
            cv = Path(tmp) / "prideqc.obo"
            MzQCWriter().write_vocabulary(result.metrics, cv, result.annotations)
            document = MzQCWriter().build(result, cv)
            import jsonschema

            schema = json.loads(
                (Path(__file__).parent / "data" / "mzqc_schema.json").read_text(encoding="utf-8")
            )
            jsonschema.Draft7Validator(
                schema, format_checker=jsonschema.FormatChecker()
            ).validate(document)
            metrics = document["mzQC"]["runQualities"][0]["qualityMetrics"]
            candidate = next(
                item for item in metrics if item["name"] == "putative modification mass shifts"
            )
            self.assertEqual(candidate["value"][0]["unimod_candidates"][0]["unimod_accession"], "UNIMOD:21")
            self.assertIn("hypotheses", candidate["description"].lower())
            self.assertIn("putative modification mass shifts", cv.read_text())

    def test_workflow_persists_per_file_and_aggregate_mass_shift_tables(self):
        collector = MassShiftCollector(
            modifications=(PHOSPHO,),
            minimum_cluster_pairs=1,
            minimum_cluster_unique_spectra=2,
        )
        for item in related_pair(0, 79.966331):
            collector.consume_spectrum(item)

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "run.raw"
            source.touch()
            result = AnalysisResult(
                source, RunMetadata(), [], collector.annotations(), [], "test-reader", 0.1
            )
            output = root / "out"
            with patch(
                "prideqc.pipeline._analyze_file",
                return_value=FileOutcome(source, result),
            ):
                manifest = Workflow(
                    WorkflowOptions(estimate_mass_shifts=True)
                ).run([source], output)

            self.assertEqual(
                manifest["files"][0]["mass_shifts"], "run.raw.mass-shifts.tsv"
            )
            self.assertTrue((output / "run.raw.mass-shifts.tsv").is_file())
            self.assertTrue((output / "mass-shifts.tsv").is_file())
            self.assertIn("UNIMOD:21", (output / "mass-shifts.tsv").read_text())

    def test_mass_shift_tsv_preserves_all_candidate_rows(self):
        annotation = MassShiftCollector(
            modifications=(
                PHOSPHO,
                ModificationRecord(
                    "UNIMOD:99998",
                    "Mass-isobaric candidate",
                    79.9664,
                    ("S",),
                    ("Anywhere",),
                    "Other",
                ),
            ),
            minimum_cluster_pairs=1,
            minimum_cluster_unique_spectra=2,
        )
        for item in related_pair(0, 79.966331):
            annotation.consume_spectrum(item)
        result = AnalysisResult(
            Path("run.raw"),
            RunMetadata(),
            [],
            annotation.annotations(),
            [],
            "test-reader",
            0.1,
        )
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "mass-shifts.tsv"
            Workflow(WorkflowOptions(estimate_mass_shifts=True))._write_mass_shift_table([result], output)
            lines = output.read_text().splitlines()
            self.assertEqual(len(lines), 3)
            self.assertIn("UNIMOD:21", lines[1] + lines[2])
            self.assertIn("UNIMOD:99998", lines[1] + lines[2])


if __name__ == "__main__":
    unittest.main()
