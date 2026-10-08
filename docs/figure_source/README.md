# prideQC scientific figure generators (maintainer-only)

The six editable SVGs in `docs/_static/figures/` are **checked-in publication assets**.
Read the Docs uses these SVGs directly; it never executes this folder, imports
pyOpenMS or generates scientific data during a documentation build.

The sources in this directory were developed for the prideQC documentation using
synthetic input and selected **implementation/private helper** functions from the
prideQC checkout (`metrics`, `mass_error`, `mass_shift` and `cohort`). Figure 2 runs
`RunSummary` / `QCMetricCalculator` on a synthetic run. Figures 4–5 use the
installed OpenMS/UniMod catalogue, so outputs can change if the catalogue changes.
All examples are synthetic; no results represent a real PRIDE dataset.

## Regenerate intentionally

From the root of the prideQC checkout, using `uv`:

```bash
# Separate from the Read the Docs Sphinx build environment.
uv venv .venv-figures --python 3.12
uv pip install --python .venv-figures/bin/python -r docs/figure_source/requirements.txt

# For a full regeneration; pyOpenMS/OpenMS is needed by Figures 4 and 5.
.venv-figures/bin/python docs/figure_source/build_all.py \
  --preview-dir docs/_build/figure-previews

# Regenerate only the fully synthetic overview/QC/mass-error/reporting figures:
.venv-figures/bin/python docs/figure_source/build_all.py --figures 1 2 3 6
```

Default SVG output is `docs/_static/figures/` regardless of current working
directory. Optional PNG previews go under the explicitly specified
`--preview-dir`; previews are **not** committed.

The font-width linter needs locally installed *Liberation Sans* and *Liberation
Mono* metrics, normally supplied by your OS. Set
`PRIDEQC_FIGURE_FONT_DIR` if those fonts live in another directory. No font
files are included here or redistributed with prideQC.

If you want to run against a different local source tree, set
`PRIDEQC_SRC=/path/to/prideQC/src`. Otherwise the generator resolves the
checkout containing this script. This is for maintainers; it is not a public
prideQC API.

## Review before commit

1. Inspect SVG and PNG previews visually at the final docs width.
2. Check generator layout-lint output; it exits nonzero when overlap issues
   are detected.
3. Confirm numeric labels/chemistry against the pinned source and catalogue.
4. Commit only intentional changed SVG files and source changes.
5. Run the strict Sphinx build after integration.

Figure captions and intended placements are recorded in `captions.md`.
