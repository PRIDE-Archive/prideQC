"""The single analysis pipeline shared by the Python API and command line."""

from __future__ import annotations

import csv
import multiprocessing
import sys
import time
from collections.abc import Iterable
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any

from prideqc import __version__
from prideqc.annotations import DiagnosticIonCollector, TechnicalAnnotator
from prideqc.conversion import ExternalConverter
from prideqc.io import atomic_text, json_safe, write_json
from prideqc.metrics import QCMetricCalculator, RunSummary
from prideqc.models import AnalysisResult, EvidenceCollector, FloatArray, Spectrum, SpectrumReader
from prideqc.mzqc import MzQCWriter
from prideqc.sdrf import SDRFDocument
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

    def __init__(self, reader: SpectrumReader | None = None,
                 calculator: QCMetricCalculator | None = None,
                 annotator: TechnicalAnnotator | None = None, *, estimate_peak_type: bool = False) -> None:
        if reader is None:
            from prideqc.readers import PyOpenMSReader

            reader = PyOpenMSReader(estimate_peak_type=estimate_peak_type)
        self.reader = reader
        self.calculator = calculator or QCMetricCalculator()
        self.annotator = annotator or TechnicalAnnotator()

    def analyze(self, path: Path | str, *, collectors: Iterable[EvidenceCollector] = ()) -> AnalysisResult:
        started = time.perf_counter()
        path = Path(path).resolve()
        collectors = tuple(collectors)
        summary = RunSummary()
        metadata = self.reader.read(path, _AnalysisSink(summary, collectors))
        annotations = self.annotator.annotate(metadata, summary)
        for collector in collectors:
            annotations.extend(collector.annotations())
        return AnalysisResult(path, metadata, self.calculator.calculate(summary), annotations,
                              summary.warnings(), self.reader.engine_version, time.perf_counter() - started)


@dataclass(frozen=True, slots=True)
class WorkflowOptions:
    workers: int = 1
    diagnostics: bool = False
    estimate_peak_type: bool = False
    converter: str | None = None
    converter_executable: str | None = None
    conversion_timeout: float = 3600
    continue_on_error: bool = False
    include_inferred: bool = False
    overwrite_sdrf_values: bool = False
    sdrf_template: str = "ms-proteomics"
    validate_ontology: bool = False
    progress: bool = False

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
        # OpenMS messages.  Flushing makes each completed file visible while
        # the remaining files are still being decoded.
        sys.stderr.write(line + "\n")
        sys.stderr.flush()

    def update(self, outcome: FileOutcome) -> None:
        if not self.enabled:
            return
        self.done += 1
        status = "ok" if outcome.error is None else "failed"
        self._write(self._line(status, outcome.source.name))

    def close(self) -> None:
        return


def _analyze_file(task: tuple[Path, Path, WorkflowOptions]) -> FileOutcome:
    source, output, options = task
    try:
        path = source
        collectors = [DiagnosticIonCollector()] if options.diagnostics else []
        from prideqc.readers import PyOpenMSReader, vendor_format

        reader = PyOpenMSReader(estimate_peak_type=options.estimate_peak_type)
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


class Workflow:
    """Produce per-file QC and batch summaries with explicit failure accounting."""

    def __init__(self, options: WorkflowOptions | None = None, *,
                 validator: SDRFValidator | None = None) -> None:
        self.options = options or WorkflowOptions()
        self.validator = validator or SDRFPipelinesValidator(
            self.options.sdrf_template, ontology=self.options.validate_ontology
        )

    def run_project(
        self, accession: str, output_directory: Path | str, *,
        download_directory: Path | str, filenames: Iterable[str] | None = None,
        sdrf: Path | None = None, aliases: dict[str, str] | None = None,
        download_options: DownloadOptions | None = None,
        repository: PrideRepository | None = None,
        username: str | None = None, password: str | None = None,
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
        if output.exists() and any(output.iterdir()):
            raise FileExistsError(f"Output directory is not empty: {output}")
        document = SDRFDocument.read(sdrf) if sdrf else None
        if filenames is None:
            filenames = [(aliases or {}).get(name, name) for name in document.data_files()] if document else []
        names = selected_files(filenames)
        if self.options.converter is None:
            vendor_names = [
                name for name in names
                if not name.lower().endswith((".mzml", ".mzml.gz"))
            ]
            if vendor_names:
                from prideqc.readers import PyOpenMSReader

                try:
                    reader = PyOpenMSReader()
                    unsupported = [name for name in vendor_names if not reader.supports_direct(Path(name))]
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
            accession, names, download, options=download_options, username=username, password=password,
        )
        manifest = self.run(paths, output, sdrf=sdrf, aliases=aliases)
        manifest["acquisition"] = {
            "accession": accession.strip().upper(),
            "manifest": str(download / "download-manifest.json"),
        }
        write_json(output / "manifest.json", manifest)
        return manifest

    def run(self, paths: Iterable[Path | str], output_directory: Path | str,
            *, sdrf: Path | None = None, aliases: dict[str, str] | None = None) -> dict[str, Any]:
        inputs = [Path(path).resolve(strict=True) for path in paths]
        if not inputs:
            raise ValueError("At least one input file is required.")
        if len({p.name.casefold() for p in inputs}) != len(inputs):
            raise ValueError(
                "Input basenames must be unique (case-insensitive) to prevent output collisions.",
            )
        if len(set(inputs)) != len(inputs):
            raise ValueError("Duplicate input paths are not allowed.")
        document = SDRFDocument.read(sdrf) if sdrf else None
        output = Path(output_directory).resolve()
        if output.exists() and any(output.iterdir()):
            raise FileExistsError(
                f"Output directory is not empty: {output}. Choose a new directory.",
            )
        input_validation = self.validator.validate(sdrf) if sdrf else None
        output.mkdir(parents=True, exist_ok=True)
        if input_validation:
            write_json(output / "sdrf-validation.json", {"input": input_validation.to_dict()})
        tasks = [(source, output, self.options) for source in inputs]
        outcomes: list[FileOutcome] = []
        progress = _Progress(len(inputs), self.options.progress)
        progress.start()
        executor = None
        if self.options.workers > 1:
            executor = ProcessPoolExecutor(max_workers=self.options.workers,
                                           mp_context=multiprocessing.get_context("spawn"))
            # Submit individually so progress reflects completion order rather
            # than waiting behind a slow first file in executor.map().
            stream = (future.result() for future in as_completed(
                [executor.submit(_analyze_file, task) for task in tasks],
            ))
        else:
            stream = map(_analyze_file, tasks)
        try:
            for outcome in stream:
                if outcome.result is not None:
                    try:
                        prefix = output / outcome.source.name
                        # Appending, not replacing, avoids .mzML/.gz stem collisions.
                        MzQCWriter().write(outcome.result, Path(f"{prefix}.mzQC"))
                        write_json(Path(f"{prefix}.summary.json"), outcome.result.to_dict())
                    except Exception as exc:
                        outcome.error = f"OutputError: {exc}"
                        outcome.result = None
                        for suffix in (".mzQC", ".obo", ".summary.json"):
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
                    ),
                )
        finally:
            progress.close()
            if executor:
                executor.shutdown(wait=True, cancel_futures=True)
        # Completion order drives progress, while reports remain deterministic
        # and unprocessed inputs are calculated from the actual completed set.
        positions = {path: index for index, path in enumerate(inputs)}
        outcomes.sort(key=lambda outcome: positions[outcome.source])
        processed = {outcome.source for outcome in outcomes}
        results = [outcome.result for outcome in outcomes if outcome.result is not None]
        manifest: dict[str, Any] = {
            "prideqc_version": __version__, "options": asdict(self.options),
            "files": [{"input": str(o.source), "status": "failed" if o.error else "success",
                       "error": o.error, "mzqc": f"{o.source.name}.mzQC" if o.result else None,
                       "elapsed_seconds": o.result.elapsed_seconds if o.result else None}
                      for o in outcomes],
            "unprocessed": [str(p) for p in inputs if p not in processed],
            "errors": [],
            "sdrf_validation": {"input": input_validation.to_dict()} if input_validation else None,
        }
        try:
            self._write_tables(results, output)
            if document:
                changes = document.annotate(results, aliases, self.options.include_inferred,
                                            self.options.overwrite_sdrf_values)
                document.write(output / "refined.sdrf.tsv")
                with atomic_text(output / "sdrf-changes.tsv") as handle:
                    writer = csv.DictWriter(
                        handle,
                        fieldnames=["row", "data_file", "column", "previous", "proposed", "status", "evidence"],
                        delimiter="\t",
                        lineterminator="\n",
                    )
                    writer.writeheader()
                    writer.writerows(asdict(change) for change in changes)
                try:
                    refined_validation = self.validator.validate(output / "refined.sdrf.tsv")
                except Exception as exc:
                    manifest["sdrf_validation"]["refined_error"] = f"{type(exc).__name__}: {exc}"
                    write_json(output / "sdrf-validation.json", manifest["sdrf_validation"])
                    raise
                manifest["sdrf_validation"]["refined"] = refined_validation.to_dict()
                write_json(output / "sdrf-validation.json", manifest["sdrf_validation"])
                if not refined_validation.valid:
                    manifest["errors"].append(
                        "Refined SDRF has sdrf-pipelines validation findings; see sdrf-validation.json."
                    )
        except Exception as exc:
            manifest["errors"].append(f"{type(exc).__name__}: {exc}")
        manifest["success"] = not (any(o.error for o in outcomes) or manifest["unprocessed"] or manifest["errors"])
        write_json(output / "manifest.json", manifest)
        return manifest

    def _write_tables(self, results: list[AnalysisResult], output: Path) -> None:
        import json

        with atomic_text(output / "metrics.tsv") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["data_file", "metric", "value"])
            for result in results:
                for metric in result.metrics:
                    writer.writerow([result.input_path.name, metric.key,
                                     json.dumps(json_safe(metric.value), allow_nan=False)])
        with atomic_text(output / "annotations.tsv") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(
                ["data_file", "field", "value", "evidence", "method", "support", "total", "detail"],
            )
            for result in results:
                for annotation in result.annotations:
                    writer.writerow([result.input_path.name, annotation.field,
                                     json.dumps(json_safe(annotation.value), allow_nan=False),
                                     annotation.kind.value, annotation.method, annotation.support,
                                     annotation.total, annotation.detail])
