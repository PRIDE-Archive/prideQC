#!/usr/bin/env python3
"""Run semantic-first v4 separately for the original prideQC GT001-GT008 cases.

Each benchmark case is scoped to one exact data file. This is intentional:
PXD017618 mixes top-down and bottom-up runs, while PXD003772 mixes label-free
and TMT runs. Accession-level aggregation would erase the original benchmark
design.

The evaluator creates symlink-only case roots pointing to existing frozen
*.mass-shifts.tsv files; it never edits or copies the science outputs.
"""

from __future__ import annotations

import argparse
import csv
import json
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path


def read_manifest(path: Path) -> dict[str, list[dict[str, str]]]:
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "gt_id", "pxd_accession", "raw_file", "sdrf_data_file",
            "unimod_accession", "evidence_status", "evidence_source", "evidence_note",
        }
        missing = required - set(reader.fieldnames or ())
        if missing:
            raise SystemExit(
                "manifest missing required columns: " + ", ".join(sorted(missing))
            )
        for row in reader:
            gt = (row.get("gt_id") or "").strip()
            if not gt:
                raise SystemExit("manifest row missing gt_id")
            grouped[gt].append(row)
    return dict(sorted(grouped.items()))


def table_data_file(path: Path) -> str:
    try:
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                value = (row.get("data_file") or "").strip()
                if value:
                    return Path(value).name
    except (OSError, UnicodeError, csv.Error):
        return ""
    return ""


def index_mass_shift_tables(result_root: Path) -> list[tuple[Path, set[str]]]:
    indexed: list[tuple[Path, set[str]]] = []
    for path in sorted(result_root.glob("**/*.mass-shifts.tsv")):
        names = {path.name.removesuffix(".mass-shifts.tsv")}
        observed = table_data_file(path)
        if observed:
            names.add(observed)
        indexed.append((path.resolve(), names))
    return indexed


def resolve_case_table(
    indexed: list[tuple[Path, set[str]]],
    *,
    accession: str,
    raw_file: str,
    sdrf_data_file: str,
) -> Path:
    wanted = {Path(raw_file).name, Path(sdrf_data_file).name}
    matches: list[Path] = []
    for path, names in indexed:
        if accession.upper() not in {part.upper() for part in path.parts}:
            continue
        if wanted & names:
            matches.append(path)
    unique = sorted(set(matches))
    if len(unique) != 1:
        detail = "\n".join(str(path) for path in unique[:20])
        raise SystemExit(
            f"{accession} {raw_file}: expected exactly one mass-shift TSV, "
            f"found {len(unique)}"
            + (f"\n{detail}" if detail else "")
        )
    return unique[0]


def write_case_evidence(path: Path, rows: list[dict[str, str]]) -> None:
    fields = [
        "pxd_accession",
        "unimod_accession",
        "evidence_status",
        "evidence_source",
        "evidence_note",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def read_review(path: Path, gt_id: str, raw_file: str) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    for row in rows:
        row["gt_id"] = gt_id
        row["benchmark_data_file"] = raw_file
    return rows


def write_summary(path: Path, rows: list[dict[str, str]]) -> None:
    preferred = [
        "gt_id",
        "benchmark_data_file",
        "pxd_accession",
        "unimod_accession",
        "unimod_name",
        "semantic_evidence_status",
        "raw_confirmation_status",
        "raw_family_mass_da",
        "raw_mass_residual_da",
        "raw_family_runs",
        "successful_runs",
        "raw_family_prevalence",
        "raw_prevalence_probability",
        "mass_identity_ambiguous",
        "mass_compatible_alternative_accessions",
        "sdrf_status",
    ]
    extras = sorted({key for row in rows for key in row} - set(preferred))
    fields = preferred + extras
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--result-root", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument(
        "--v4-script",
        type=Path,
        default=Path(__file__).with_name("calibrate_mass_shift_probability_v4.py"),
    )
    parser.add_argument("--catalog-source", choices=("openms", "observed"), default="openms")
    parser.add_argument("--family-tolerance-da", type=float, default=0.02)
    parser.add_argument("--raw-match-window-da", type=float, default=0.02)
    parser.add_argument("--minimum-raw-prevalence", type=float, default=0.10)
    parser.add_argument(
        "--raw-prevalence-probability-threshold", type=float, default=0.95
    )
    parser.add_argument(
        "--minimum-family-runs",
        type=int,
        default=1,
        help=(
            "GT001-GT008 are one-file benchmark cases, so the evaluation default "
            "is 1. Production v4 remains 2."
        ),
    )
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    result_root = args.result_root.resolve()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    groups = read_manifest(args.manifest.resolve())
    indexed = index_mass_shift_tables(result_root)
    if not indexed:
        raise SystemExit(f"no *.mass-shifts.tsv files found below {result_root}")

    resolved: list[dict[str, str]] = []
    combined: list[dict[str, str]] = []
    for gt_id, rows in groups.items():
        accession = rows[0]["pxd_accession"].strip().upper()
        raw_file = rows[0]["raw_file"].strip()
        sdrf_data_file = rows[0]["sdrf_data_file"].strip()
        if any(row["pxd_accession"].strip().upper() != accession for row in rows):
            raise SystemExit(f"{gt_id}: mixed accessions in manifest")
        source = resolve_case_table(
            indexed,
            accession=accession,
            raw_file=raw_file,
            sdrf_data_file=sdrf_data_file,
        )
        resolved.append(
            {
                "gt_id": gt_id,
                "pxd_accession": accession,
                "raw_file": raw_file,
                "source_mass_shift_tsv": str(source),
            }
        )
        print(f"{gt_id}\t{accession}\t{raw_file}\t{source}")

        if args.dry_run:
            continue

        case_root = output / "cases" / gt_id
        input_dir = case_root / "input" / accession
        v4_out = case_root / "v4"
        if case_root.exists():
            shutil.rmtree(case_root)
        input_dir.mkdir(parents=True)
        v4_out.mkdir(parents=True)

        link = input_dir / source.name
        link.symlink_to(source)

        evidence = case_root / "study-evidence.tsv"
        write_case_evidence(evidence, rows)

        command = [
            sys.executable,
            str(args.v4_script.resolve()),
            "--result-root", str(case_root / "input"),
            "--output-dir", str(v4_out),
            "--study-evidence-tsv", str(evidence),
            "--catalog-source", args.catalog_source,
            "--family-tolerance-da", str(args.family_tolerance_da),
            "--raw-match-window-da", str(args.raw_match_window_da),
            "--minimum-raw-prevalence", str(args.minimum_raw_prevalence),
            "--raw-prevalence-probability-threshold",
            str(args.raw_prevalence_probability_threshold),
            "--minimum-family-runs", str(args.minimum_family_runs),
        ]
        completed = subprocess.run(
            command,
            check=False,
            text=True,
            capture_output=True,
        )
        (case_root / "v4.stdout.txt").write_text(completed.stdout, encoding="utf-8")
        (case_root / "v4.stderr.txt").write_text(completed.stderr, encoding="utf-8")
        if completed.returncode != 0:
            raise SystemExit(
                f"{gt_id}: v4 failed with rc={completed.returncode}\n"
                f"{completed.stderr}"
            )
        review = v4_out / "mass-shift-probability-v4.sdrf-review.tsv"
        combined.extend(read_review(review, gt_id, raw_file))

    resolved_path = output / "original-gt8-resolved-inputs.tsv"
    with resolved_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["gt_id","pxd_accession","raw_file","source_mass_shift_tsv"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(resolved)

    if args.dry_run:
        print(f"resolved_inputs={resolved_path}")
        print(f"cases={len(resolved)}")
        return

    summary = output / "original-gt8-v4-summary.tsv"
    write_summary(summary, combined)
    status_counts: dict[str, int] = defaultdict(int)
    for row in combined:
        status_counts[row.get("sdrf_status", "")] += 1
    metadata = {
        "evaluation": "original-prideqc-GT001-GT008-case-scoped-v4",
        "result_root": str(result_root),
        "manifest": str(args.manifest.resolve()),
        "case_count": len(groups),
        "semantic_candidate_count": len(combined),
        "minimum_family_runs": args.minimum_family_runs,
        "note": (
            "One-file benchmark diagnostic. minimum_family_runs=1 is evaluation-only; "
            "production v4 remains at 2."
        ),
        "sdrf_status_counts": dict(sorted(status_counts.items())),
    }
    metadata_path = output / "original-gt8-v4-metadata.json"
    metadata_path.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    print(f"summary={summary}")
    print(f"metadata={metadata_path}")
    print(f"cases={len(groups)}")
    print(f"semantic_candidates={len(combined)}")
    for status, count in sorted(status_counts.items()):
        print(f"sdrf_status_{status.replace('-', '_')}={count}")


if __name__ == "__main__":
    main()
