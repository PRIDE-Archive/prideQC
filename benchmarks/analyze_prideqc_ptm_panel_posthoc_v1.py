#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import statistics
from collections import Counter, defaultdict
from pathlib import Path

TARGETS = {
    "PXD000138": [
        ("phosphorylation", 79.966331, 96),
    ],
    "PXD041271": [
        ("phosphorylation", 79.966331, 67),
    ],
    "PXD027496": [
        ("acetylation", 42.010565, 11),
    ],
    "PXD012814": [
        ("succinylation", 100.016044, 30),
    ],
    "PXD022005": [
        ("succinylation", 100.016044, 41),
    ],
    "PXD014574": [
        ("mono-methylation", 14.015650, 28),
        ("di-methylation", 28.031300, 28),
        ("tri-methylation", 42.046950, 28),
    ],
    "PXD022367": [
        ("ubiquitin-diglycine-remnant", 114.042927, 169),
    ],
}


def as_float(value: str | None) -> float | None:
    if value in (None, ""):
        return None
    try:
        return float(value)
    except ValueError:
        return None


def as_int(value: str | None) -> int | None:
    if value in (None, ""):
        return None
    try:
        return int(float(value))
    except ValueError:
        return None


def first_nonempty(row: dict[str, str], *names: str) -> str:
    for name in names:
        value = row.get(name, "")
        if value not in ("", None):
            return str(value)
    return ""


def median_or_blank(values):
    values = [v for v in values if v is not None]
    return statistics.median(values) if values else ""


def range_or_blank(values):
    values = [v for v in values if v is not None]
    if not values:
        return "", ""
    return min(values), max(values)


def cluster_from_group(rank: int, rows: list[dict[str, str]]) -> dict:
    row = rows[0]

    candidate_accessions = sorted(
        {
            first_nonempty(
                item,
                "unimod_accession",
                "candidate_accession",
                "unimod_candidate_accession",
            )
            for item in rows
            if first_nonempty(
                item,
                "unimod_accession",
                "candidate_accession",
                "unimod_candidate_accession",
            )
        }
    )
    candidate_names = sorted(
        {
            first_nonempty(
                item,
                "unimod_name",
                "candidate_name",
                "unimod_candidate_name",
            )
            for item in rows
            if first_nonempty(
                item,
                "unimod_name",
                "candidate_name",
                "unimod_candidate_name",
            )
        }
    )

    return {
        "rank": rank,
        "delta_mass_da": as_float(row.get("delta_mass_da")),
        "match_tolerance_da": as_float(row.get("match_tolerance_da")),
        "pair_support": as_int(row.get("pair_support")),
        "unique_spectrum_support": as_int(row.get("unique_spectrum_support")),
        "support_fraction": as_float(
            first_nonempty(
                row,
                "support_fraction_of_accepted_pairs",
                "support_fraction",
            )
        ),
        "median_spectral_similarity": as_float(
            row.get("median_spectral_similarity")
        ),
        "classification": row.get("classification", ""),
        "confidence": row.get("confidence", ""),
        "candidate_count": as_int(
            first_nonempty(
                row,
                "diagnostic_unimod_candidate_count",
                "unimod_candidate_count",
            )
        ),
        "candidate_accessions": candidate_accessions,
        "candidate_names": candidate_names,
    }


def read_clusters(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        groups: dict[int, list[dict[str, str]]] = defaultdict(list)

        for row in reader:
            rank = as_int(row.get("cluster_rank"))
            if rank is None:
                continue
            groups[rank].append(row)

    clusters = [
        cluster_from_group(rank, rows)
        for rank, rows in sorted(groups.items())
    ]
    return [
        cluster
        for cluster in clusters
        if cluster["delta_mass_da"] is not None
    ]


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Post-hoc PTM recovery audit for the frozen prideQC PTM panel."
    )
    parser.add_argument(
        "--result-root",
        required=True,
        type=Path,
        help="Frozen run root, e.g. .../results/prideqc-v23_1-ptm-panel-full-v1",
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        type=Path,
        help="Output directory outside the frozen result tree",
    )
    args = parser.parse_args()

    result_root = args.result_root.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    per_run_rows = []

    for accession, targets in TARGETS.items():
        paths = sorted(
            (result_root / accession).glob(
                "repository-metadata/**/*.raw.mass-shifts.tsv"
            )
        )

        for path in paths:
            clusters = read_clusters(path)
            raw_name = path.name.removesuffix(".mass-shifts.tsv")

            for target_name, target_mass, expected_manifest_runs in targets:
                nearest = None
                compatible = []

                for cluster in clusters:
                    residual = abs(cluster["delta_mass_da"] - target_mass)
                    item = (residual, cluster)

                    if nearest is None or residual < nearest[0]:
                        nearest = item

                    tolerance = cluster["match_tolerance_da"]
                    if tolerance is not None and residual <= tolerance:
                        compatible.append(item)

                compatible.sort(key=lambda item: item[0])

                if compatible:
                    residual, chosen = compatible[0]
                    is_compatible = True
                elif nearest is not None:
                    residual, chosen = nearest
                    is_compatible = False
                else:
                    residual = None
                    chosen = {
                        "rank": None,
                        "delta_mass_da": None,
                        "match_tolerance_da": None,
                        "pair_support": None,
                        "unique_spectrum_support": None,
                        "support_fraction": None,
                        "median_spectral_similarity": None,
                        "classification": "",
                        "confidence": "",
                        "candidate_count": None,
                        "candidate_accessions": [],
                        "candidate_names": [],
                    }
                    is_compatible = False

                per_run_rows.append(
                    {
                        "pxd_accession": accession,
                        "raw_file": raw_name,
                        "target_name": target_name,
                        "target_mass_da": target_mass,
                        "expected_manifest_runs": expected_manifest_runs,
                        "compatible": int(is_compatible),
                        "compatible_cluster_count": len(compatible),
                        "cluster_rank": chosen["rank"],
                        "observed_delta_mass_da": chosen["delta_mass_da"],
                        "residual_da": residual,
                        "match_tolerance_da": chosen["match_tolerance_da"],
                        "pair_support": chosen["pair_support"],
                        "unique_spectrum_support": chosen[
                            "unique_spectrum_support"
                        ],
                        "support_fraction": chosen["support_fraction"],
                        "median_spectral_similarity": chosen[
                            "median_spectral_similarity"
                        ],
                        "classification": chosen["classification"],
                        "confidence": chosen["confidence"],
                        "candidate_count": chosen["candidate_count"],
                        "candidate_accessions": " | ".join(
                            chosen["candidate_accessions"]
                        ),
                        "candidate_names": " | ".join(
                            chosen["candidate_names"]
                        ),
                        "source_tsv": str(path),
                    }
                )

    per_run_path = output_dir / "ptm-posthoc-v1.per-run.tsv"
    per_run_fields = [
        "pxd_accession",
        "raw_file",
        "target_name",
        "target_mass_da",
        "expected_manifest_runs",
        "compatible",
        "compatible_cluster_count",
        "cluster_rank",
        "observed_delta_mass_da",
        "residual_da",
        "match_tolerance_da",
        "pair_support",
        "unique_spectrum_support",
        "support_fraction",
        "median_spectral_similarity",
        "classification",
        "confidence",
        "candidate_count",
        "candidate_accessions",
        "candidate_names",
        "source_tsv",
    ]

    with per_run_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=per_run_fields,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(per_run_rows)

    grouped = defaultdict(list)
    for row in per_run_rows:
        grouped[
            (
                row["pxd_accession"],
                row["target_name"],
                row["target_mass_da"],
                row["expected_manifest_runs"],
            )
        ].append(row)

    summary_rows = []

    for (
        accession,
        target_name,
        target_mass,
        expected_manifest_runs,
    ), rows in sorted(grouped.items()):
        successful_runs = len(rows)
        hits = [row for row in rows if row["compatible"] == 1]
        hit_count = len(hits)

        ranks = [as_int(str(row["cluster_rank"])) for row in hits]
        residuals = [as_float(str(row["residual_da"])) for row in hits]
        pairs = [as_int(str(row["pair_support"])) for row in hits]
        spectra = [
            as_int(str(row["unique_spectrum_support"]))
            for row in hits
        ]
        support_fractions = [
            as_float(str(row["support_fraction"]))
            for row in hits
        ]
        similarities = [
            as_float(str(row["median_spectral_similarity"]))
            for row in hits
        ]

        rank_min, rank_max = range_or_blank(ranks)
        classification_counts = Counter(
            row["classification"] or "(blank)" for row in hits
        )
        confidence_counts = Counter(
            row["confidence"] or "(blank)" for row in hits
        )

        summary_rows.append(
            {
                "pxd_accession": accession,
                "target_name": target_name,
                "target_mass_da": target_mass,
                "expected_manifest_runs": expected_manifest_runs,
                "successful_science_runs": successful_runs,
                "compatible_runs": hit_count,
                "coverage_successful_fraction": (
                    hit_count / successful_runs if successful_runs else 0.0
                ),
                "coverage_manifest_fraction": (
                    hit_count / expected_manifest_runs
                    if expected_manifest_runs
                    else 0.0
                ),
                "median_rank": median_or_blank(ranks),
                "rank_min": rank_min,
                "rank_max": rank_max,
                "median_residual_da": median_or_blank(residuals),
                "median_pair_support": median_or_blank(pairs),
                "median_unique_spectrum_support": median_or_blank(spectra),
                "median_support_fraction": median_or_blank(
                    support_fractions
                ),
                "median_spectral_similarity": median_or_blank(similarities),
                "classification_counts": "; ".join(
                    f"{key}={value}"
                    for key, value in sorted(classification_counts.items())
                ),
                "confidence_counts": "; ".join(
                    f"{key}={value}"
                    for key, value in sorted(confidence_counts.items())
                ),
            }
        )

    summary_path = output_dir / "ptm-posthoc-v1.summary.tsv"
    summary_fields = [
        "pxd_accession",
        "target_name",
        "target_mass_da",
        "expected_manifest_runs",
        "successful_science_runs",
        "compatible_runs",
        "coverage_successful_fraction",
        "coverage_manifest_fraction",
        "median_rank",
        "rank_min",
        "rank_max",
        "median_residual_da",
        "median_pair_support",
        "median_unique_spectrum_support",
        "median_support_fraction",
        "median_spectral_similarity",
        "classification_counts",
        "confidence_counts",
    ]

    with summary_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=summary_fields,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(summary_rows)

    print(f"result_root={result_root}")
    print(f"per_run={per_run_path}")
    print(f"summary={summary_path}")
    print(f"per_run_rows={len(per_run_rows)}")
    print(f"summary_rows={len(summary_rows)}")


if __name__ == "__main__":
    main()
