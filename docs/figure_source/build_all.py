"""Regenerate the six public prideQC figures from synthetic data and source code.

Normal Read the Docs builds use the committed SVG files; this script runs only when
maintainers deliberately regenerate them. It imports prideQC implementation helpers
and OpenMS/UniMod catalog data; review output before committing changed SVG files.

Examples (from the repository root):

    python docs/figure_source/build_all.py
    python docs/figure_source/build_all.py --preview-dir docs/_build/figure-previews
    python docs/figure_source/build_all.py --figures 1 2 3
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import fig1_overview
import fig2_qc_metrics
import fig3_mass_error
import fig4_mass_shift
import fig5_cohort_sdrf
import fig6_reporting

FIGURES = {
    1: ("figure1_overview", fig1_overview),
    2: ("figure2_qc_metrics", fig2_qc_metrics),
    3: ("figure3_mass_error", fig3_mass_error),
    4: ("figure4_mass_shift", fig4_mass_shift),
    5: ("figure5_cohort_sdrf", fig5_cohort_sdrf),
    6: ("figure6_reporting", fig6_reporting),
}

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_SVG_DIR = ROOT / "docs" / "_static" / "figures"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_SVG_DIR)
    parser.add_argument("--preview-dir", type=Path, default=None)
    parser.add_argument(
        "--figures", type=int, nargs="+", choices=FIGURES.keys(), default=list(FIGURES)
    )
    args = parser.parse_args(argv)

    preview = None
    if args.preview_dir is not None:
        from render import render  # CairoSVG only needed for preview generation

        preview = render
        args.preview_dir.mkdir(parents=True, exist_ok=True)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    bad = 0
    for n in args.figures:
        name, module = FIGURES[n]
        try:
            fig = module.build()
        except Exception as exc:
            raise RuntimeError(
                f"Unable to build {name}. The figure uses selected prideQC helpers "
                "and figures 4–5 require the installed OpenMS/UniMod catalogue "
                "(pyopenms). Ensure the figure requirements are installed and "
                "the local prideQC checkout matches the scientific semantics."
            ) from exc
        problems = fig.lint(name)
        for problem in problems:
            print(problem, file=sys.stderr)
        bad += len(problems)
        svg_path = args.output_dir / f"{name}.svg"
        fig.save(svg_path)
        if preview is not None:
            preview(svg_path, args.preview_dir / f"{name}.png", scale=1.5)
        print(
            f"{name}: {fig.w:g}×{fig.h:g}, {len(fig.texts)} text items, {len(problems)} layout warnings"
        )
    if bad:
        print(f"Layout lint: {bad} issues; review the SVGs before committing.", file=sys.stderr)
        return 1
    print("All selected figure layouts passed lint.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
