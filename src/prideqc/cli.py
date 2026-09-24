"""Command line entry point with lazy imports for analysis and remote acquisition."""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path
from typing import Any

from prideqc import __version__


def _validation_arguments(command: argparse.ArgumentParser) -> None:
    command.add_argument("--sdrf-template", default="ms-proteomics", help="sdrf-pipelines template")
    command.add_argument(
        "--validate-ontology", action="store_true",
        help="Also validate ontology terms (requires uv sync --extra ontology)",
    )


def _download_arguments(command: argparse.ArgumentParser) -> None:
    command.add_argument(
        "--file",
        action="append",
        dest="remote_files",
        help="Exact PRIDE filename; repeat for a subset",
    )
    command.add_argument("--protocol", choices=["ftp", "aspera", "s3", "globus"], default="ftp")
    command.add_argument(
        "--no-checksum-check",
        action="store_true",
        help="Disable pridepy checksum checking",
    )
    command.add_argument("--aspera-maximum-bandwidth", default="100M")
    command.add_argument("--username", help="Username for private PRIDE access")
    command.add_argument(
        "--password-env",
        default="PRIDE_PASSWORD",
        help="Environment variable holding the password",
    )


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(
        prog="prideqc", description="Per-file mass spectrometry QC and technical SDRF evidence",
    )
    result.add_argument("--version", action="version", version=f"prideQC {__version__}")
    commands = result.add_subparsers(dest="command", required=True)
    analyze = commands.add_parser(
        "analyze",
        help="Generate mzQC and annotation reports from local or PRIDE files",
    )
    analyze.add_argument("files", nargs="*", type=Path)
    analyze.add_argument(
        "-o",
        "--output-dir",
        required=True,
        type=Path,
        help="Output directory (must be empty unless --overwrite is given)",
    )
    analyze.add_argument("--sdrf", type=Path)
    analyze.add_argument("--accession", help="Download selected PRIDE files before analysis")
    analyze.add_argument(
        "--download-dir",
        type=Path,
        help="New or empty directory separate from QC output",
    )
    _download_arguments(analyze)
    _validation_arguments(analyze)
    analyze.add_argument(
        "--file-map",
        type=Path,
        help="JSON object: exact SDRF filename to analyzed filename",
    )
    analyze.add_argument(
        "--workers",
        type=int,
        default=1,
        help="Independent file processes; default 1",
    )
    analyze.add_argument(
        "--diagnostics",
        action="store_true",
        help="Screen centroid MS2/MS3 for candidate reporter/oxonium signatures",
    )
    analyze.add_argument(
        "--estimate-mass-error",
        action="store_true",
        help="Estimate repeat-spectrum precursor/fragment mass-error precision",
    )
    analyze.add_argument(
        "--estimate-mass-shifts",
        action="store_true",
        help=(
            "Screen centroid MS2 for recurrent related-spectrum neutral mass shifts and "
            "annotate mass-compatible OpenMS/UniMod candidates"
        ),
    )
    analyze.add_argument(
        "--estimate-peak-type",
        action="store_true",
        help="Also estimate centroid/profile type with OpenMS",
    )
    analyze.add_argument(
        "--include-inferred",
        action="store_true",
        help="Allow acquisition-mode heuristic suggestions to fill SDRF cells",
    )
    analyze.add_argument(
        "--overwrite-sdrf-values",
        action="store_true",
        help="Replace conflicting SDRF values and record each change",
    )
    analyze.add_argument(
        "--refine-sdrf-qc",
        action="store_true",
        help=(
            "Write confidence-gated cohort tolerance/PTM evidence to standard SDRF columns; "
            "preserve original modifications and record every proposal/change"
        ),
    )
    analyze.add_argument(
        "--ptm-study-evidence",
        type=Path,
        help=(
            "Optional independent PTM semantic-evidence TSV; exact UniMod support is required "
            "before recurrent mass-shift evidence is written as an SDRF modification"
        ),
    )
    analyze.add_argument(
        "--project-accession",
        help="Explicit ProteomeXchange accession for local cohort refinement provenance",
    )
    analyze.add_argument(
        "--continue-on-error",
        action="store_true",
        help="Process remaining local inputs after analysis failure; exit code still nonzero",
    )
    analyze.add_argument(
        "--overwrite",
        action="store_true",
        help="Clear an existing output directory before writing results",
    )
    analyze.add_argument(
        "--no-progress",
        dest="progress",
        action="store_false",
        help="Disable the live per-file progress display",
    )
    analyze.add_argument(
        "--converter",
        choices=["thermorawfileparser", "msconvert"],
        help="Optional installed vendor converter",
    )
    analyze.add_argument(
        "--converter-executable",
        help="Executable name or path (no shell command)",
    )
    analyze.add_argument("--conversion-timeout", type=float, default=3600)
    validate = commands.add_parser(
        "check-sdrf", aliases=["validate-sdrf"], help="Validate SDRF with sdrf-pipelines templates",
    )
    validate.add_argument("file", type=Path)
    validate.add_argument(
        "--json",
        type=Path,
        dest="report",
        help="Write validation results as JSON",
    )
    _validation_arguments(validate)
    refine = commands.add_parser(
        "refine-sdrf-qc",
        help="Refine one SDRF from previously generated prideQC per-file summaries",
    )
    refine.add_argument(
        "--results-root",
        required=True,
        type=Path,
        help="Directory containing one accession/run set of *.summary.json artifacts",
    )
    refine.add_argument("--sdrf", required=True, type=Path)
    refine.add_argument(
        "-o",
        "--output-dir",
        required=True,
        type=Path,
        help="New or empty cohort-refinement output directory",
    )
    refine.add_argument(
        "--file-map",
        type=Path,
        help="Optional JSON object: exact SDRF filename to analyzed filename",
    )
    refine.add_argument(
        "--ptm-study-evidence",
        type=Path,
        help=(
            "Optional independent PTM semantic-evidence TSV using pxd_accession, "
            "unimod_accession and evidence_status columns"
        ),
    )
    refine.add_argument(
        "--project-accession",
        help="Explicit ProteomeXchange accession when older summaries do not contain it",
    )
    refine.add_argument(
        "--overwrite",
        action="store_true",
        help="Clear an existing output directory before writing refinement artifacts",
    )
    _validation_arguments(refine)

    llm = commands.add_parser(
        "llm",
        help="Manage and run the optional local LLM used for SDRF adjudication",
    )
    llm_commands = llm.add_subparsers(dest="llm_command", required=True)
    llm_setup = llm_commands.add_parser(
        "setup",
        help="Download and verify the pinned llama.cpp runtime and default GGUF model",
    )
    llm_setup.add_argument("--cache-dir", type=Path)
    llm_setup.add_argument(
        "--force",
        action="store_true",
        help="Redownload and replace the managed runtime and model",
    )
    llm_status = llm_commands.add_parser(
        "status",
        help="Show whether the managed local runtime and model are installed",
    )
    llm_status.add_argument("--cache-dir", type=Path)
    llm_adjudicate = llm_commands.add_parser(
        "adjudicate",
        help="Adjudicate a prideQC request with the managed local llama.cpp model",
    )
    llm_adjudicate.add_argument("--request", required=True, type=Path)
    llm_adjudicate.add_argument(
        "-o",
        "--output",
        type=Path,
        help="Validated decision JSON (default: beside request as llm-refinement-decisions.json)",
    )
    llm_adjudicate.add_argument(
        "--audit-output",
        type=Path,
        help="Raw model/audit JSON (default: beside request as llm-model-run.json)",
    )
    llm_adjudicate.add_argument("--cache-dir", type=Path)
    llm_adjudicate.add_argument("--server-path", type=Path)
    llm_adjudicate.add_argument("--model-path", type=Path)
    llm_adjudicate.add_argument("--startup-timeout", type=float, default=180.0)
    llm_adjudicate.add_argument("--request-timeout", type=float, default=900.0)
    llm_adjudicate.add_argument("--context-size", type=int, default=8192)

    fetch = commands.add_parser("fetch", help="Download an explicit PRIDE selection using pridepy")
    fetch.add_argument("accession")
    fetch.add_argument(
        "-o",
        "--output-dir",
        required=True,
        type=Path,
        help="New or empty download directory",
    )
    fetch.add_argument(
        "--sdrf",
        type=Path,
        help="Select filenames from SDRF when --file is omitted",
    )
    _download_arguments(fetch)
    return result


def _download_settings(arguments: argparse.Namespace) -> dict[str, Any]:
    from prideqc.pride import DownloadOptions

    password = os.environ.get(arguments.password_env) if arguments.username else None
    if arguments.username and not password:
        raise ValueError(f"Set the password environment variable {arguments.password_env!r}.")
    return {
        "options": DownloadOptions(
            protocol=arguments.protocol, checksum_check=not arguments.no_checksum_check,
            aspera_maximum_bandwidth=arguments.aspera_maximum_bandwidth,
        ),
        "username": arguments.username, "password": password,
    }


def _analyze(arguments: argparse.Namespace) -> int:
    # Set before importing NumPy/OpenMS; explicit user settings take precedence.
    for variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ.setdefault(variable, "1")
    from prideqc.pipeline import Workflow, WorkflowOptions

    aliases = None
    if arguments.file_map:
        aliases = json.loads(arguments.file_map.read_text(encoding="utf-8"))
        if not isinstance(aliases, dict) or not all(
            isinstance(key, str) and isinstance(value, str) for key, value in aliases.items()
        ):
            raise ValueError("--file-map must contain a JSON object of string-to-string mappings.")
        if not arguments.sdrf:
            raise ValueError("--file-map requires --sdrf.")
    if arguments.converter_executable and not arguments.converter:
        raise ValueError("--converter-executable requires --converter.")
    workflow = Workflow(WorkflowOptions(
        workers=arguments.workers, diagnostics=arguments.diagnostics,
        estimate_mass_error=arguments.estimate_mass_error,
        estimate_mass_shifts=arguments.estimate_mass_shifts,
        estimate_peak_type=arguments.estimate_peak_type, converter=arguments.converter,
        converter_executable=arguments.converter_executable,
        conversion_timeout=arguments.conversion_timeout, continue_on_error=arguments.continue_on_error,
        include_inferred=arguments.include_inferred, overwrite_sdrf_values=arguments.overwrite_sdrf_values,
        refine_sdrf_qc=arguments.refine_sdrf_qc,
        ptm_study_evidence=(
            str(arguments.ptm_study_evidence) if arguments.ptm_study_evidence else None
        ),
        project_accession=arguments.project_accession or arguments.accession,
        sdrf_template=arguments.sdrf_template, validate_ontology=arguments.validate_ontology,
        progress=arguments.progress, overwrite=arguments.overwrite,
    ))
    if arguments.accession:
        if arguments.files or not arguments.download_dir:
            raise ValueError(
                "--accession requires --download-dir and cannot be combined with local files.",
            )
        settings = _download_settings(arguments)
        manifest = workflow.run_project(
            arguments.accession, arguments.output_dir, download_directory=arguments.download_dir,
            filenames=arguments.remote_files, sdrf=arguments.sdrf, aliases=aliases,
            download_options=settings["options"], username=settings["username"], password=settings["password"],
        )
    else:
        if not arguments.files:
            raise ValueError("Provide local files or --accession with --file/--sdrf.")
        if (arguments.download_dir or arguments.remote_files or arguments.username
                or arguments.no_checksum_check or arguments.protocol != "ftp"
                or arguments.aspera_maximum_bandwidth != "100M" or arguments.password_env != "PRIDE_PASSWORD"):
            raise ValueError("PRIDE download options require --accession.")
        manifest = workflow.run(
            arguments.files,
            arguments.output_dir,
            sdrf=arguments.sdrf,
            aliases=aliases,
        )
    succeeded = sum(item["status"] == "success" for item in manifest["files"])
    total = len(manifest["files"]) + len(manifest.get("unprocessed", []))
    print(f"Analyzed {succeeded}/{total} files. Results: {arguments.output_dir}")
    for item in manifest["files"]:
        if item["error"]:
            print(f"{item['input']}: {item['error']}", file=sys.stderr)
    for error in manifest["errors"]:
        print(error, file=sys.stderr)
    return 0 if manifest["success"] else 1


def _refine_sdrf_qc(arguments: argparse.Namespace) -> int:
    from prideqc.pipeline import Workflow, WorkflowOptions

    aliases = None
    if arguments.file_map:
        aliases = json.loads(arguments.file_map.read_text(encoding="utf-8"))
        if not isinstance(aliases, dict) or not all(
            isinstance(key, str) and isinstance(value, str) for key, value in aliases.items()
        ):
            raise ValueError("--file-map must contain a JSON object of string-to-string mappings.")
    workflow = Workflow(
        WorkflowOptions(
            refine_sdrf_qc=True,
            ptm_study_evidence=(
                str(arguments.ptm_study_evidence) if arguments.ptm_study_evidence else None
            ),
            project_accession=arguments.project_accession,
            overwrite=arguments.overwrite,
            sdrf_template=arguments.sdrf_template,
            validate_ontology=arguments.validate_ontology,
        )
    )
    manifest = workflow.refine_existing_sdrf(
        arguments.results_root,
        arguments.output_dir,
        sdrf=arguments.sdrf,
        aliases=aliases,
    )
    print(
        f"Refined {manifest['summary_count']} files into "
        f"{manifest['experiment_groups']} experiment group(s). "
        f"Results: {arguments.output_dir}"
    )
    return 0 if manifest["success"] else 1


def _llm(arguments: argparse.Namespace) -> int:
    if arguments.llm_command == "setup":
        from prideqc.local_llm import setup_local_llm

        status = setup_local_llm(
            arguments.cache_dir, force=arguments.force, show_progress=True
        )
        print(json.dumps(status, indent=2, sort_keys=True))
        return 0
    if arguments.llm_command == "status":
        from prideqc.local_llm import local_llm_status

        status = local_llm_status(arguments.cache_dir)
        print(json.dumps(status, indent=2, sort_keys=True))
        return 0 if status["runtime"]["ready"] and status["model"]["ready"] else 1
    if arguments.llm_command == "adjudicate":
        from prideqc.io import write_json
        from prideqc.refinement_adjudication import validate_llm_adjudication_request
        from prideqc.refinement_model import LocalLlamaCppAdapter

        request = json.loads(arguments.request.read_text(encoding="utf-8"))
        if not isinstance(request, dict):
            raise ValueError("--request must contain one JSON object")
        validate_llm_adjudication_request(request)
        if arguments.context_size < 2048:
            raise ValueError("--context-size must be at least 2048")
        if arguments.startup_timeout <= 0 or arguments.request_timeout <= 0:
            raise ValueError("LLM timeouts must be positive")
        output = arguments.output or arguments.request.with_name("llm-refinement-decisions.json")
        audit_output = arguments.audit_output or arguments.request.with_name("llm-model-run.json")
        def report_progress(index: int, total: int, decision_id: str) -> None:
            print(f"Adjudicating {index}/{total}: {decision_id}")

        adapter = LocalLlamaCppAdapter(
            cache_dir=arguments.cache_dir,
            server_path=arguments.server_path,
            model_path=arguments.model_path,
            startup_timeout=arguments.startup_timeout,
            request_timeout=arguments.request_timeout,
            context_size=arguments.context_size,
            progress=report_progress,
        )
        result = adapter.adjudicate(request)
        write_json(output, result.decisions)
        write_json(audit_output, result.audit)
        print(f"Validated LLM decisions: {output}")
        print(f"LLM audit: {audit_output}")
        return 0
    raise ValueError(f"Unknown llm command: {arguments.llm_command}")


def main(argv: list[str] | None = None) -> int:
    arguments = parser().parse_args(argv)
    try:
        if arguments.command in {"check-sdrf", "validate-sdrf"}:
            from prideqc.io import write_json
            from prideqc.validation import SDRFPipelinesValidator

            report = SDRFPipelinesValidator(
                arguments.sdrf_template, ontology=arguments.validate_ontology,
            ).validate(arguments.file)
            if arguments.report:
                write_json(arguments.report, report.to_dict())
            if report.valid:
                print(f"SDRF valid for {report.template} (ontology validation: {report.ontology}).")
            for issue in report.issues:
                print(issue, file=sys.stderr)
            return 0 if report.valid else 1
        if arguments.command == "refine-sdrf-qc":
            return _refine_sdrf_qc(arguments)
        if arguments.command == "llm":
            return _llm(arguments)
        if arguments.command == "fetch":
            from prideqc.pride import PrideRepository
            from prideqc.sdrf import SDRFDocument

            names = arguments.remote_files
            if names is None and arguments.sdrf:
                names = SDRFDocument.read(arguments.sdrf).data_files()
            paths = PrideRepository().download_files(
                arguments.accession, names or [], arguments.output_dir, **_download_settings(arguments),
            )
            print(f"Downloaded {len(paths)} files to {arguments.output_dir}.")
            return 0
        return _analyze(arguments)
    except (OSError, ValueError, RuntimeError) as exc:
        print(f"prideqc: {exc}", file=sys.stderr)
        return 2
