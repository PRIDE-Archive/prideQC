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
        estimate_peak_type=arguments.estimate_peak_type, converter=arguments.converter,
        converter_executable=arguments.converter_executable,
        conversion_timeout=arguments.conversion_timeout, continue_on_error=arguments.continue_on_error,
        include_inferred=arguments.include_inferred, overwrite_sdrf_values=arguments.overwrite_sdrf_values,
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
