"""The single analysis pipeline shared by the Python API and command line."""

from __future__ import annotations

import csv
import multiprocessing
import re
import shutil
import sys
import time
from collections.abc import Iterable
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any

from prideqc import __version__
from prideqc.annotations import DiagnosticIonCollector, TechnicalAnnotator
from prideqc.cohort import (
    COHORT_OVERWRITE_FIELDS,
    COHORT_SDRF_FIELDS,
    PUTATIVE_MODIFICATION_COLUMN,
    SemanticEvidence,
    read_semantic_evidence,
    synthesize_cohort,
)
from prideqc.conversion import ExternalConverter
from prideqc.io import atomic_text, json_safe, write_json
from prideqc.mass_error import RepeatSpectrumMassErrorCollector
from prideqc.mass_shift import MassShiftCollector
from prideqc.metrics import QCMetricCalculator, RunSummary
from prideqc.models import AnalysisResult, EvidenceCollector, FloatArray, Spectrum, SpectrumReader
from prideqc.mzqc import MzQCWriter
from prideqc.refinement_adjudication import (
    build_llm_adjudication_request,
    validate_llm_adjudication_request,
)
from prideqc.refinement_packet import (
    build_llm_refinement_packet,
    validate_llm_refinement_packet,
)
from prideqc.sdrf import SDRFChange, SDRFDocument
from prideqc.validation import SDRFPipelinesValidator, SDRFValidator

if TYPE_CHECKING:
    from prideqc.pride import DownloadOptions, PrideRepository


class _AnalysisSink:
    def __init__(self, summary: RunSummary, collectors: Iterable[EvidenceCollector]) -> None:
        self.summary = summary
        self.collectors = tuple(collectors)

    def consume_spectrum(self, spectrum: Spectrum) -> None:
        self.summary.consume_spectrum(spectrum)
        for collector in self.collectors:
            collector.consume_spectrum(spectrum)

    def consume_chromatogram(self, rt: FloatArray, kind: int) -> None:
        self.summary.consume_chromatogram(rt, kind)


class Analyzer:
    """Inject a reader or evidence collector without coupling it to mzQC/SDRF."""

    def __init__(
        self,
        reader: SpectrumReader | None = None,
        calculator: QCMetricCalculator | None = None,
        annotator: TechnicalAnnotator | None = None,
        *,
        estimate_peak_type: bool = False,
    ) -> None:
        if reader is None:
            from prideqc.readers import PyOpenMSReader

            reader = PyOpenMSReader(estimate_peak_type=estimate_peak_type)
        self.reader = reader
        self.calculator = calculator or QCMetricCalculator()
        self.annotator = annotator or TechnicalAnnotator()

    def analyze(
        self,
        path: Path | str,
        *,
        collectors: Iterable[EvidenceCollector] = (),
    ) -> AnalysisResult:
        started = time.perf_counter()
        path = Path(path).resolve()
        collectors = tuple(collectors)
        summary = RunSummary()
        metadata = self.reader.read(path, _AnalysisSink(summary, collectors))
        annotations = self.annotator.annotate(metadata, summary)
        for collector in collectors:
            annotations.extend(collector.annotations())
        return AnalysisResult(
            path,
            metadata,
            self.calculator.calculate(summary),
            annotations,
            summary.warnings(),
            self.reader.engine_version,
            time.perf_counter() - started,
        )


@dataclass(frozen=True, slots=True)
class WorkflowOptions:
    workers: int = 1
    diagnostics: bool = False
    estimate_mass_error: bool = False
    estimate_mass_shifts: bool = False
    estimate_peak_type: bool = False
    converter: str | None = None
    converter_executable: str | None = None
    conversion_timeout: float = 3600
    continue_on_error: bool = False
    include_inferred: bool = False
    overwrite_sdrf_values: bool = False
    refine_sdrf_qc: bool = False
    ptm_study_evidence: str | None = None
    project_accession: str | None = None
    sdrf_template: str = "ms-proteomics"
    validate_ontology: bool = False
    progress: bool = False
    overwrite: bool = False

    def __post_init__(self) -> None:
        if self.workers < 1:
            raise ValueError("workers must be >=1")
        if self.conversion_timeout <= 0:
            raise ValueError("conversion_timeout must be positive")
        if not self.sdrf_template.strip():
            raise ValueError("sdrf_template must be nonempty")


@dataclass(slots=True)
class FileOutcome:
    source: Path
    result: AnalysisResult | None = None
    error: str | None = None


class _Progress:
    """Dependency-free progress display updated as workers finish files."""

    def __init__(self, total: int, enabled: bool) -> None:
        self.total = total
        self.done = 0
        self.enabled = enabled

    def start(self) -> None:
        if not self.enabled:
            return
        self._write(self._line("starting", ""))

    def _line(self, status: str, filename: str) -> str:
        width = 20
        filled = round(width * self.done / self.total) if self.total else width
        bar = "#" * filled + "-" * (width - filled)
        suffix = f" ({status}: {filename})" if filename else ""
        return f"Analyzing files: {self.done}/{self.total} [{bar}]{suffix}"

    def _write(self, line: str) -> None:
        # A newline per event is deliberate: carriage-return repainting is
        # commonly hidden or interleaved by CI/log collectors and native
        # OpenMS messages. Flushing makes each completed file visible while
        # the remaining files are still being decoded.
        sys.stderr.write(line + "\n")
        sys.stderr.flush()

    def update(self, outcome: FileOutcome) -> None:
        if not self.enabled:
            return
        self.done += 1
        status = "ok" if outcome.error is None else "failed"
        self._write(self._line(status, outcome.source.name))

    def close(self, *, unprocessed: int = 0) -> None:
        if self.enabled and unprocessed:
            self._write(
                f"Analysis stopped: {unprocessed} file(s) were not processed after the first failure."
            )


def _analyze_file(task: tuple[Path, Path, WorkflowOptions]) -> FileOutcome:
    source, output, options = task
    try:
        path = source
        collectors: list[EvidenceCollector] = []
        if options.diagnostics:
            collectors.append(DiagnosticIonCollector())
        mass_error_collector = None
        if options.estimate_mass_error:
            mass_error_collector = RepeatSpectrumMassErrorCollector()
            collectors.append(mass_error_collector)
        if options.estimate_mass_shifts:
            collectors.append(MassShiftCollector(precision_source=mass_error_collector))
        from prideqc.readers import PyOpenMSReader, vendor_format

        reader = PyOpenMSReader(
            estimate_peak_type=(
                options.estimate_peak_type
                or options.estimate_mass_error
                or options.estimate_mass_shifts
            ),
        )
        is_mzml = source.name.casefold().endswith((".mzml", ".mzml.gz"))
        detected_vendor = vendor_format(source)
        if not is_mzml and (detected_vendor is None or not reader.supports_direct(source)):
            if options.converter is None:
                # Keep this as a reader-level capability error so users get the
                # installed version and upgrade/converter guidance.
                if detected_vendor is not None:
                    raise reader.unavailable_error(source)
                raise ValueError("Vendor data requires --converter or prior conversion to mzML.")
            path = ExternalConverter(
                options.converter,
                options.converter_executable,
                options.conversion_timeout,
            ).convert(source, output / "converted" / source.name)
        result = Analyzer(reader=reader).analyze(path, collectors=collectors)
        if path != source:
            result.source_path = source
        return FileOutcome(source, result)
    except Exception as exc:
        # This is the process/file boundary. A failed file never gets a success result.
        return FileOutcome(source, error=f"{type(exc).__name__}: {exc}")


def _path_project_accession(paths: Iterable[Path]) -> str | None:
    accessions: set[str] = set()
    for path in paths:
        accessions.update(match.upper() for match in re.findall(r"PXD\d+", str(path), re.I))
    if len(accessions) > 1:
        raise ValueError(
            "Multiple ProteomeXchange accessions found in refinement paths: "
            + ", ".join(sorted(accessions))
        )
    return next(iter(accessions), None)


class Workflow:
    """Produce per-file QC and batch summaries with explicit failure accounting."""

    def __init__(
        self,
        options: WorkflowOptions | None = None,
        *,
        validator: SDRFValidator | None = None,
    ) -> None:
        self.options = options or WorkflowOptions()
        self.validator = validator or SDRFPipelinesValidator(
            self.options.sdrf_template,
            ontology=self.options.validate_ontology,
        )

    def _ptm_semantic_evidence(
        self,
    ) -> dict[tuple[str, str], SemanticEvidence] | None:
        if not self.options.ptm_study_evidence:
            return None
        path = Path(self.options.ptm_study_evidence).resolve(strict=True)
        return read_semantic_evidence(path)

    def _cohort_project_accession(
        self,
        results: Iterable[AnalysisResult],
        *,
        extra_paths: Iterable[Path] = (),
    ) -> str | None:
        result_accessions = {
            str(result.project_accession).strip().upper()
            for result in results
            if result.project_accession
        }
        if len(result_accessions) > 1:
            raise ValueError(
                "Results contain multiple ProteomeXchange accessions: "
                + ", ".join(sorted(result_accessions))
            )
        observed = next(iter(result_accessions), None)
        explicit = (self.options.project_accession or "").strip().upper() or None
        inferred = _path_project_accession(extra_paths)
        candidates = {value for value in (observed, explicit, inferred) if value}
        if len(candidates) > 1:
            raise ValueError(
                "Conflicting ProteomeXchange accessions for cohort refinement: "
                + ", ".join(sorted(candidates))
            )
        return next(iter(candidates), None)

    def run_project(
        self,
        accession: str,
        output_directory: Path | str,
        *,
        download_directory: Path | str,
        filenames: Iterable[str] | None = None,
        sdrf: Path | None = None,
        aliases: dict[str, str] | None = None,
        download_options: DownloadOptions | None = None,
        repository: PrideRepository | None = None,
        username: str | None = None,
        password: str | None = None,
    ) -> dict[str, Any]:
        """Download an explicit PRIDE selection, then use the same local workflow.

        Selection defaults to the SDRF data-file cells. Repeated sample rows
        download once. Download failures stop before analysis and leave the
        acquisition manifest plus any successfully downloaded files.
        """
        from prideqc.pride import PrideRepository, selected_files

        output = Path(output_directory).resolve()
        download = Path(download_directory).resolve()
        if output.is_relative_to(download) or download.is_relative_to(output):
            raise ValueError("QC and download directories must be separate, non-nested directories.")
        if output.exists() and not output.is_dir():
            raise FileExistsError(f"Output path is not a directory: {output}")
        if output.exists() and any(output.iterdir()) and not self.options.overwrite:
            raise FileExistsError(f"Output directory is not empty: {output}")
        document = SDRFDocument.read(sdrf) if sdrf else None
        if filenames is None:
            filenames = [
                (aliases or {}).get(name, name) for name in document.data_files()
            ] if document else []
        names = selected_files(filenames)
        if self.options.converter is None:
            vendor_names = [
                name for name in names if not name.lower().endswith((".mzml", ".mzml.gz"))
            ]
            if vendor_names:
                from prideqc.readers import PyOpenMSReader

                try:
                    reader = PyOpenMSReader()
                    unsupported = [
                        name for name in vendor_names if not reader.supports_direct(Path(name))
                    ]
                except RuntimeError:
                    unsupported = vendor_names
                if unsupported:
                    raise ValueError(
                        "Selected vendor files require a native pyOpenMS reader or "
                        "--converter: " + ", ".join(unsupported),
                    )
        if sdrf:
            # Resolve validator/template/ontology setup before potentially large
            # downloads. run() validates the actual input again before analysis.
            self.validator.validate(sdrf)
        paths = (repository or PrideRepository()).download_files(
            accession,
            names,
            download,
            options=download_options,
            username=username,
            password=password,
        )
        normalized_accession = accession.strip().upper()
        manifest = self.run(
            paths,
            output,
            sdrf=sdrf,
            aliases=aliases,
            project_accession=normalized_accession,
        )
        manifest["acquisition"] = {
            "accession": normalized_accession,
            "manifest": str(download / "download-manifest.json"),
        }
        write_json(output / "manifest.json", manifest)
        return manifest

    def run(
        self,
        paths: Iterable[Path | str],
        output_directory: Path | str,
        *,
        sdrf: Path | None = None,
        aliases: dict[str, str] | None = None,
        project_accession: str | None = None,
    ) -> dict[str, Any]:
        inputs = [Path(path).resolve(strict=True) for path in paths]
        if not inputs:
            raise ValueError("At least one input file is required.")
        if len({p.name.casefold() for p in inputs}) != len(inputs):
            raise ValueError(
                "Input basenames must be unique (case-insensitive) to prevent output collisions."
            )
        if len(set(inputs)) != len(inputs):
            raise ValueError("Duplicate input paths are not allowed.")
        if self.options.refine_sdrf_qc and sdrf is None:
            raise ValueError("--refine-sdrf-qc requires an input SDRF.")
        document = SDRFDocument.read(sdrf) if sdrf else None
        output = Path(output_directory).resolve()
        if any(
            output == source or output.is_relative_to(source) or source.is_relative_to(output)
            for source in inputs
        ):
            raise ValueError("Output directory must be separate from all input files/directories.")
        if output.exists() and not output.is_dir():
            raise FileExistsError(f"Output path is not a directory: {output}")
        if output.exists() and any(output.iterdir()):
            if not self.options.overwrite:
                raise FileExistsError(
                    f"Output directory is not empty: {output}. Choose a new directory "
                    "or pass --overwrite."
                )
            # The user opted in explicitly. Remove only children of the
            # requested results directory; never remove the directory itself
            # or anything outside it.
            for child in output.iterdir():
                if child.is_dir() and not child.is_symlink():
                    shutil.rmtree(child)
                else:
                    child.unlink()
        input_validation = self.validator.validate(sdrf) if sdrf else None
        output.mkdir(parents=True, exist_ok=True)
        if self.options.refine_sdrf_qc and sdrf is not None:
            shutil.copyfile(sdrf, output / "original.sdrf.tsv")
        if input_validation:
            write_json(
                output / "sdrf-validation.json",
                {"input": input_validation.to_dict()},
            )
        tasks = [(source, output, self.options) for source in inputs]
        outcomes: list[FileOutcome] = []
        progress = _Progress(len(inputs), self.options.progress)
        progress.start()
        executor = None
        stream: Iterable[FileOutcome]
        if self.options.workers > 1:
            executor = ProcessPoolExecutor(
                max_workers=self.options.workers,
                mp_context=multiprocessing.get_context("spawn"),
            )
            # Submit individually so progress reflects completion order rather
            # than waiting behind a slow first file in executor.map().
            stream = (
                future.result()
                for future in as_completed(
                    [executor.submit(_analyze_file, task) for task in tasks]
                )
            )
        else:
            stream = map(_analyze_file, tasks)
        try:
            for outcome in stream:
                if outcome.result is not None:
                    try:
                        if project_accession:
                            outcome.result.project_accession = project_accession.strip().upper()
                        prefix = output / outcome.source.name
                        # Appending, not replacing, avoids .mzML/.gz stem collisions.
                        MzQCWriter().write(outcome.result, Path(f"{prefix}.mzQC"))
                        write_json(
                            Path(f"{prefix}.summary.json"),
                            outcome.result.to_dict(),
                        )
                        if self.options.estimate_mass_shifts:
                            self._write_mass_shift_table(
                                [outcome.result],
                                Path(f"{prefix}.mass-shifts.tsv"),
                            )
                    except Exception as exc:
                        outcome.error = f"OutputError: {exc}"
                        outcome.result = None
                        for suffix in (".mzQC", ".obo", ".summary.json", ".mass-shifts.tsv"):
                            Path(f"{prefix}{suffix}").unlink(missing_ok=True)
                outcomes.append(outcome)
                progress.update(outcome)
                if outcome.error and not self.options.continue_on_error:
                    break
        except Exception as exc:
            # A native worker crash must also leave a reviewable failure manifest.
            if len(outcomes) < len(inputs):
                outcomes.append(
                    FileOutcome(
                        inputs[len(outcomes)],
                        error=f"WorkerError: {type(exc).__name__}: {exc}",
                    )
                )
        finally:
            if executor:
                executor.shutdown(wait=True, cancel_futures=True)
            progress.close(unprocessed=len(inputs) - len(outcomes))
        # Completion order drives progress, while reports remain deterministic
        # and unprocessed inputs are calculated from the actual completed set.
        positions = {path: index for index, path in enumerate(inputs)}
        outcomes.sort(key=lambda outcome: positions[outcome.source])
        processed = {outcome.source for outcome in outcomes}
        results = [outcome.result for outcome in outcomes if outcome.result is not None]
        cohort_synthesis = None
        manifest: dict[str, Any] = {
            "prideqc_version": __version__,
            "options": asdict(self.options),
            "files": [
                {
                    "input": str(o.source),
                    "status": "failed" if o.error else "success",
                    "error": o.error,
                    "mzqc": f"{o.source.name}.mzQC" if o.result else None,
                    "mass_shifts": (
                        f"{o.source.name}.mass-shifts.tsv"
                        if o.result and self.options.estimate_mass_shifts
                        else None
                    ),
                    "elapsed_seconds": o.result.elapsed_seconds if o.result else None,
                }
                for o in outcomes
            ],
            "unprocessed": [str(p) for p in inputs if p not in processed],
            "errors": [],
            "sdrf_validation": {"input": input_validation.to_dict()}
            if input_validation
            else None,
            "cohort_refinement": {"enabled": self.options.refine_sdrf_qc},
        }
        try:
            if self.options.refine_sdrf_qc and results:
                if document is None or sdrf is None:
                    raise ValueError("Cohort SDRF refinement requires an input SDRF.")
                cohort_results, missing = document.resolve_results(results, aliases)
                if missing:
                    raise ValueError(
                        "Cohort SDRF refinement requires complete analyzed coverage; missing: "
                        + ", ".join(missing)
                    )
                project_accession = self._cohort_project_accession(
                    cohort_results, extra_paths=(sdrf,) if sdrf is not None else ()
                )
                semantic_evidence = self._ptm_semantic_evidence()
                cohort_synthesis = synthesize_cohort(
                    cohort_results,
                    semantic_evidence=semantic_evidence,
                    project_accession=project_accession,
                )
                write_json(output / "cohort-refinement.json", cohort_synthesis.to_dict())
                packet_artifacts: dict[str, dict[str, str]] = {}
                for outcome in outcomes:
                    if outcome.result is None:
                        continue
                    artifact_prefix = outcome.source.name
                    refs = {
                        "summary_json": f"{artifact_prefix}.summary.json",
                        "mzqc": f"{artifact_prefix}.mzQC",
                    }
                    if self.options.estimate_mass_shifts:
                        refs["mass_shifts_tsv"] = f"{artifact_prefix}.mass-shifts.tsv"
                    packet_artifacts[outcome.result.input_path.name] = refs
                packet = build_llm_refinement_packet(
                    document,
                    cohort_results,
                    cohort_synthesis,
                    project_accession=project_accession,
                    sdrf_path=sdrf.resolve(strict=True),
                    prideqc_version=__version__,
                    aliases=aliases,
                    source_artifacts=packet_artifacts,
                    semantic_evidence=semantic_evidence,
                )
                validate_llm_refinement_packet(packet)
                write_json(output / "llm-refinement-packet.json", packet)
                adjudication_request = build_llm_adjudication_request(packet)
                validate_llm_adjudication_request(adjudication_request)
                write_json(output / "llm-adjudication-request.json", adjudication_request)
                manifest["cohort_refinement"] = {
                    "enabled": True,
                    "artifact": "cohort-refinement.json",
                    "llm_refinement_packet": "llm-refinement-packet.json",
                    "llm_adjudication_request": "llm-adjudication-request.json",
                    "original_sdrf": "original.sdrf.tsv",
                    "experiment_groups": len(cohort_synthesis.groups),
                    "sdrf_eligible_ptm_families": len(cohort_synthesis.ptm_families),
                    "ptm_review_families": len(cohort_synthesis.ptm_review_families),
                }
                # Per-file mzQC/summary files were written as soon as each analysis
                # completed. Rewrite successful outputs once so cohort annotations
                # become part of the same provenance consumed by downstream tools.
                for outcome in outcomes:
                    if outcome.result is None:
                        continue
                    prefix = output / outcome.source.name
                    MzQCWriter().write(outcome.result, Path(f"{prefix}.mzQC"))
                    write_json(Path(f"{prefix}.summary.json"), outcome.result.to_dict())
            self._write_tables(results, output)
            if document:
                changes = document.annotate(
                    results,
                    aliases,
                    self.options.include_inferred,
                    self.options.overwrite_sdrf_values,
                    include_inferred_fields=(
                        COHORT_SDRF_FIELDS if self.options.refine_sdrf_qc else frozenset()
                    ),
                    overwrite_fields=(
                        COHORT_OVERWRITE_FIELDS
                        if self.options.refine_sdrf_qc
                        else frozenset()
                    ),
                    append_columns=(
                        frozenset(
                            {
                                "comment[modification parameters]",
                                PUTATIVE_MODIFICATION_COLUMN,
                            }
                        )
                        if self.options.refine_sdrf_qc
                        else frozenset()
                    ),
                )
                document.write(output / "refined.sdrf.tsv")
                self._write_sdrf_change_outputs(
                    changes,
                    output,
                    cohort_synthesis.to_dict() if cohort_synthesis is not None else None,
                )
                try:
                    refined_validation = self.validator.validate(
                        output / "refined.sdrf.tsv"
                    )
                except Exception as exc:
                    manifest["sdrf_validation"]["refined_error"] = (
                        f"{type(exc).__name__}: {exc}"
                    )
                    write_json(
                        output / "sdrf-validation.json",
                        manifest["sdrf_validation"],
                    )
                    raise
                manifest["sdrf_validation"]["refined"] = refined_validation.to_dict()
                write_json(
                    output / "sdrf-validation.json",
                    manifest["sdrf_validation"],
                )
                if not refined_validation.valid:
                    manifest["errors"].append(
                        "Refined SDRF has sdrf-pipelines validation findings; "
                        "see sdrf-validation.json."
                    )
        except Exception as exc:
            manifest["errors"].append(f"{type(exc).__name__}: {exc}")
        manifest["success"] = not (
            any(o.error for o in outcomes)
            or manifest["unprocessed"]
            or manifest["errors"]
        )
        write_json(output / "manifest.json", manifest)
        return manifest

    def refine_existing_sdrf(
        self,
        results_root: Path | str,
        output_directory: Path | str,
        *,
        sdrf: Path,
        aliases: dict[str, str] | None = None,
    ) -> dict[str, Any]:
        """Refine one SDRF from previously generated per-file summary artifacts.

        This is the accession-level companion to the Slurm one-file-per-task
        workflow: all successful ``*.summary.json`` files are loaded first, then
        experiment groups, shared tolerances and strict recurrent PTM families
        are synthesized once across the whole cohort.
        """
        import json

        root = Path(results_root).resolve()
        output = Path(output_directory).resolve()
        sdrf = Path(sdrf).resolve(strict=True)
        if not root.is_dir():
            raise ValueError(f"Results root is not a directory: {root}")
        if output == root or root.is_relative_to(output):
            raise ValueError(
                "Refinement output cannot equal or contain the results root. "
                "A subdirectory below the results root is allowed."
            )
        if output.exists() and not output.is_dir():
            raise FileExistsError(f"Output path is not a directory: {output}")
        if output.exists() and any(output.iterdir()):
            if not self.options.overwrite:
                raise FileExistsError(
                    f"Output directory is not empty: {output}. Choose a new directory "
                    "or pass --overwrite."
                )
            for child in output.iterdir():
                if child.is_dir() and not child.is_symlink():
                    shutil.rmtree(child)
                else:
                    child.unlink()

        summary_paths = sorted(root.rglob("*.summary.json"))
        if not summary_paths:
            raise ValueError(f"No *.summary.json files found below {root}")
        results: list[AnalysisResult] = []
        packet_artifacts: dict[str, dict[str, str]] = {}
        for summary_path in summary_paths:
            payload = json.loads(summary_path.read_text(encoding="utf-8"))
            if not isinstance(payload, dict):
                raise ValueError(f"Invalid prideQC summary object: {summary_path}")
            result = AnalysisResult.from_dict(payload)
            results.append(result)
            refs = {"summary_json": summary_path.relative_to(root).as_posix()}
            prefix = summary_path.name.removesuffix(".summary.json")
            mzqc_path = summary_path.with_name(f"{prefix}.mzQC")
            if mzqc_path.exists():
                refs["mzqc"] = mzqc_path.relative_to(root).as_posix()
            mass_shift_path = summary_path.with_name(f"{prefix}.mass-shifts.tsv")
            if mass_shift_path.exists():
                refs["mass_shifts_tsv"] = mass_shift_path.relative_to(root).as_posix()
            packet_artifacts[result.input_path.name] = refs
        names = [result.input_path.name.casefold() for result in results]
        if len(set(names)) != len(names):
            raise ValueError(
                "Duplicate analyzed input basenames found below results root; "
                "select one accession/run set."
            )
        project_accession = self._cohort_project_accession(
            results, extra_paths=(*summary_paths, sdrf)
        )

        input_validation = self.validator.validate(sdrf)
        document = SDRFDocument.read(sdrf)
        output.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(sdrf, output / "original.sdrf.tsv")

        cohort_results, missing = document.resolve_results(results, aliases)
        if missing:
            raise ValueError(
                "Cohort SDRF refinement requires complete analyzed coverage; missing: "
                + ", ".join(missing)
            )
        semantic_evidence = self._ptm_semantic_evidence()
        synthesis = synthesize_cohort(
            cohort_results,
            semantic_evidence=semantic_evidence,
            project_accession=project_accession,
        )
        write_json(output / "cohort-refinement.json", synthesis.to_dict())
        packet = build_llm_refinement_packet(
            document,
            cohort_results,
            synthesis,
            project_accession=project_accession,
            sdrf_path=sdrf,
            prideqc_version=__version__,
            aliases=aliases,
            source_artifacts=packet_artifacts,
            semantic_evidence=semantic_evidence,
        )
        validate_llm_refinement_packet(packet)
        write_json(output / "llm-refinement-packet.json", packet)
        adjudication_request = build_llm_adjudication_request(packet)
        validate_llm_adjudication_request(adjudication_request)
        write_json(output / "llm-adjudication-request.json", adjudication_request)
        changes = document.annotate(
            cohort_results,
            aliases,
            include_inferred=False,
            overwrite=False,
            include_inferred_fields=COHORT_SDRF_FIELDS,
            overwrite_fields=COHORT_OVERWRITE_FIELDS,
            append_columns=frozenset(
                {
                    "comment[modification parameters]",
                    PUTATIVE_MODIFICATION_COLUMN,
                }
            ),
        )
        document.write(output / "refined.sdrf.tsv")
        self._write_sdrf_change_outputs(changes, output, synthesis.to_dict())
        refined_validation = self.validator.validate(output / "refined.sdrf.tsv")
        validation = {
            "input": input_validation.to_dict(),
            "refined": refined_validation.to_dict(),
        }
        write_json(output / "sdrf-validation.json", validation)
        manifest = {
            "prideqc_version": __version__,
            "mode": "existing-results-sdrf-refinement",
            "results_root": str(root),
            "summary_files": [str(path) for path in summary_paths],
            "summary_count": len(cohort_results),
            "discovered_summary_count": len(summary_paths),
            "project_accession": project_accession,
            "ptm_study_evidence": self.options.ptm_study_evidence,
            "original_sdrf": "original.sdrf.tsv",
            "refined_sdrf": "refined.sdrf.tsv",
            "cohort_refinement": "cohort-refinement.json",
            "llm_refinement_packet": "llm-refinement-packet.json",
            "llm_adjudication_request": "llm-adjudication-request.json",
            "sdrf_changes": "sdrf-changes.tsv",
            "sdrf_refinement_log": "sdrf-refinement.log.txt",
            "experiment_groups": len(synthesis.groups),
            "sdrf_eligible_ptm_families": len(synthesis.ptm_families),
            "ptm_review_families": len(synthesis.ptm_review_families),
            "putative_ptm_families_written": len(synthesis.ptm_review_families),
            "success": refined_validation.valid,
        }
        write_json(output / "manifest.json", manifest)
        return manifest

    def _write_sdrf_change_outputs(
        self,
        changes: list[SDRFChange],
        output: Path,
        cohort: dict[str, Any] | None,
    ) -> None:
        fields = [
            "row",
            "data_file",
            "column",
            "previous",
            "proposed",
            "status",
            "evidence",
            "annotation_field",
            "experiment_group",
            "method",
            "detail",
            "support",
            "total",
        ]
        with atomic_text(output / "sdrf-changes.tsv") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            writer.writerows(asdict(change) for change in changes)

        applied = [
            change
            for change in changes
            if change.status in {"filled", "replaced", "appended"}
        ]
        review = [
            change
            for change in changes
            if change.status in {"conflict", "suggestion", "ambiguous_column", "unmatched"}
        ]
        with atomic_text(output / "sdrf-refinement.log.txt") as handle:
            handle.write("prideQC SDRF refinement log\n")
            handle.write("===========================\n")
            handle.write(f"applied_changes={len(applied)}\n")
            handle.write(f"review_items={len(review)}\n")
            if cohort is not None:
                groups = cohort.get("groups", {})
                ptms = cohort.get("ptm_families", [])
                ptm_review = cohort.get("ptm_review_families", [])
                handle.write(f"experiment_groups={len(groups)}\n")
                handle.write(f"sdrf_eligible_ptm_families={len(ptms)}\n")
                handle.write(f"ptm_review_families={len(ptm_review)}\n")
                for label, summary in groups.items():
                    handle.write(
                        f"group={label} files={summary.get('files', 0)} "
                        f"members={','.join(summary.get('members', []))}\n"
                    )
            handle.write("\n[applied]\n")
            for change in applied:
                handle.write(
                    f"row={change.row} data_file={change.data_file} status={change.status} "
                    f"column={change.column} previous={change.previous!r} "
                    f"proposed={change.proposed!r} field={change.annotation_field} "
                    f"group={change.experiment_group!r} support={change.support}/{change.total} "
                    f"method={change.method}\n"
                )
            handle.write("\n[review_or_not_applied]\n")
            for change in review:
                handle.write(
                    f"row={change.row} data_file={change.data_file} status={change.status} "
                    f"column={change.column} previous={change.previous!r} "
                    f"proposed={change.proposed!r} field={change.annotation_field} "
                    f"group={change.experiment_group!r}\n"
                )
            if cohort is not None:
                handle.write("\n[ptm_review_families]\n")
                for item in cohort.get("ptm_review_families", []):
                    handle.write(
                        f"group={item.get('experiment_group', '')!r} "
                        f"unimod={item.get('unimod_accession', '')} "
                        f"name={item.get('unimod_name', '')!r} "
                        f"mass_da={item.get('median_mass_da')} "
                        f"runs={item.get('family_runs')}/{item.get('group_runs')} "
                        f"raw_prevalence_probability={item.get('raw_prevalence_probability')} "
                        f"semantic_status={item.get('semantic_evidence_status', '')} "
                        f"sdrf_status={item.get('sdrf_status', '')}\n"
                    )

    def _write_tables(self, results: list[AnalysisResult], output: Path) -> None:
        import json

        with atomic_text(output / "metrics.tsv") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["data_file", "metric", "value"])
            for result in results:
                for metric in result.metrics:
                    writer.writerow(
                        [
                            result.input_path.name,
                            metric.key,
                            json.dumps(json_safe(metric.value), allow_nan=False),
                        ]
                    )
        with atomic_text(output / "annotations.tsv") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(
                [
                    "data_file",
                    "field",
                    "value",
                    "evidence",
                    "method",
                    "support",
                    "total",
                    "detail",
                ]
            )
            for result in results:
                for annotation in result.annotations:
                    writer.writerow(
                        [
                            result.input_path.name,
                            annotation.field,
                            json.dumps(json_safe(annotation.value), allow_nan=False),
                            annotation.kind.value,
                            annotation.method,
                            annotation.support,
                            annotation.total,
                            annotation.detail,
                        ]
                    )
        if self.options.estimate_mass_shifts:
            self._write_mass_shift_table(results, output / "mass-shifts.tsv")

    def _write_mass_shift_table(
        self,
        results: list[AnalysisResult],
        path: Path,
    ) -> None:
        """Write one row per recurrent cluster / plausible UniMod candidate."""

        fields = [
            "data_file",
            "cluster_rank",
            "delta_mass_da",
            "cluster_sigma_da",
            "cluster_min_da",
            "cluster_max_da",
            "pair_support",
            "unique_spectrum_support",
            "median_spectral_similarity",
            "classification",
            "confidence",
            "match_tolerance_da",
            "support_fraction_of_accepted_pairs",
            "diagnostic_unimod_candidate_count",
            "diagnostic_suppressed_unimod_candidates",
            "artifact_name",
            "artifact_theoretical_delta_mass_da",
            "artifact_residual_da",
            "candidate_rank",
            "unimod_accession",
            "unimod_name",
            "unimod_theoretical_delta_mass_da",
            "unimod_residual_da",
            "unimod_origins",
            "unimod_term_specificities",
            "unimod_source_classification",
            "unimod_candidate_category",
            "orientation",
        ]
        with atomic_text(path) as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=fields,
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            for result in results:
                annotation = next(
                    (
                        item
                        for item in result.annotations
                        if item.field == "putative_modification_mass_shifts"
                    ),
                    None,
                )
                clusters = annotation.value if annotation and isinstance(annotation.value, list) else []
                for cluster_rank, cluster in enumerate(clusters, start=1):
                    artifact = cluster.get("artifact_candidate") or {}
                    candidates = cluster.get("unimod_candidates") or [None]
                    for candidate_rank, candidate in enumerate(candidates, start=1):
                        candidate = candidate or {}
                        writer.writerow({
                            "data_file": result.input_path.name,
                            "cluster_rank": cluster_rank,
                            "delta_mass_da": cluster.get("delta_mass_da"),
                            "cluster_sigma_da": cluster.get("cluster_sigma_da"),
                            "cluster_min_da": cluster.get("cluster_min_da"),
                            "cluster_max_da": cluster.get("cluster_max_da"),
                            "pair_support": cluster.get("pair_support"),
                            "unique_spectrum_support": cluster.get("unique_spectrum_support"),
                            "median_spectral_similarity": cluster.get(
                                "median_spectral_similarity"
                            ),
                            "classification": cluster.get("classification"),
                            "confidence": cluster.get("confidence"),
                            "match_tolerance_da": cluster.get("match_tolerance_da"),
                            "support_fraction_of_accepted_pairs": cluster.get(
                                "support_fraction_of_accepted_pairs"
                            ),
                            "diagnostic_unimod_candidate_count": cluster.get(
                                "diagnostic_unimod_candidate_count"
                            ),
                            "diagnostic_suppressed_unimod_candidates": cluster.get(
                                "diagnostic_suppressed_unimod_candidates"
                            ),
                            "artifact_name": artifact.get("name"),
                            "artifact_theoretical_delta_mass_da": artifact.get(
                                "theoretical_delta_mass_da"
                            ),
                            "artifact_residual_da": artifact.get("residual_da"),
                            "candidate_rank": candidate_rank if candidate else "",
                            "unimod_accession": candidate.get("unimod_accession"),
                            "unimod_name": candidate.get("name"),
                            "unimod_theoretical_delta_mass_da": candidate.get(
                                "theoretical_delta_mass_da"
                            ),
                            "unimod_residual_da": candidate.get("residual_da"),
                            "unimod_origins": ";".join(candidate.get("origins") or []),
                            "unimod_term_specificities": ";".join(
                                candidate.get("term_specificities") or []
                            ),
                            "unimod_source_classification": candidate.get(
                                "source_classification"
                            ),
                            "unimod_candidate_category": candidate.get(
                                "candidate_category"
                            ),
                            "orientation": cluster.get("orientation"),
                        })
