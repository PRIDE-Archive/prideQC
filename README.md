# prideQC

<!-- prideQC project badges -->
[![Python CI](https://github.com/PRIDE-Archive/prideQC/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/PRIDE-Archive/prideQC/actions/workflows/ci.yml)
[![Docker](https://github.com/PRIDE-Archive/prideQC/actions/workflows/container-docker.yml/badge.svg?branch=main)](https://github.com/PRIDE-Archive/prideQC/actions/workflows/container-docker.yml)
[![Apptainer SIF](https://github.com/PRIDE-Archive/prideQC/actions/workflows/container-sif.yml/badge.svg?branch=main)](https://github.com/PRIDE-Archive/prideQC/actions/workflows/container-sif.yml)
[![Documentation](https://img.shields.io/badge/docs-Read%20the%20Docs-8CA1AF?logo=readthedocs&logoColor=white)](https://prideqc.readthedocs.io/en/latest/)
[![Python 3.11–3.13](https://img.shields.io/badge/python-3.11%20%7C%203.12%20%7C%203.13-3776AB?logo=python&logoColor=white)](pyproject.toml)
[![License: Apache-2.0 AND MIT](https://img.shields.io/badge/license-Apache--2.0%20AND%20MIT-2C7A7B)](LICENSE)

<!-- After the first PyPI publication, optionally add a dynamic release badge:
[![PyPI version](https://img.shields.io/pypi/v/prideQC?logo=pypi&logoColor=white)](https://pypi.org/project/prideQC/)
-->

**prideQC** computes quality-control metrics from raw mass-spectrometry data and uses that measured evidence to refine technical SDRF metadata.

It is designed around a simple separation of responsibilities:

1. **Measure QC evidence from the data.**
2. **Summarize and report that evidence.**
3. **Refine metadata only when the study-level evidence supports it.**

## Documentation

- Documentation: <https://prideqc.readthedocs.io/>
- Source: <https://github.com/PRIDE-Archive/prideQC>
- PRIDE Archive: <https://www.ebi.ac.uk/pride/>

Until the hosted documentation is enabled, the same material is available under `docs/` in this repository.

## Highlights

prideQC can:

- analyze local mzML and supported raw/vendor inputs;
- compute per-file mass-spectrometry QC metrics;
- estimate precursor and fragment measurement precision;
- scout recurrent neutral mass shifts and mass-compatible modification candidates;
- write structured JSON, TSV, and mzQC-compatible outputs;
- aggregate results for interactive pmultiqc / MultiQC reporting;
- synthesize QC evidence across a study/cohort;
- propose conservative, auditable SDRF refinements;
- validate refined SDRFs with `sdrf-pipelines`.

## Quick start

```bash
prideqc analyze sample1.mzML sample2.mzML \
  --output-dir results \
  --estimate-mass-error \
  --estimate-mass-shifts \
  --diagnostics
```

Refine an SDRF from existing QC results:

```bash
prideqc refine-sdrf-qc \
  --results-root results \
  --sdrf study.sdrf.tsv \
  --output-dir refinement
```

Validate an SDRF:

```bash
prideqc check-sdrf study.sdrf.tsv
```

## Reporting

prideQC's mzQC-compatible outputs can be aggregated into interactive HTML reports with pmultiqc / MultiQC:

```bash
multiqc --module mzqc \
  --outdir report \
  results/
```

See the documentation for the complete QC, reporting, and SDRF-refinement workflow.

## Local documentation build

```bash
uv venv .venv-docs --python 3.12
uv pip install --python .venv-docs/bin/python -r docs/requirements.txt
.venv-docs/bin/sphinx-build -W --keep-going -b html docs docs/_build/html
```

The six scientific SVG figures are checked into `docs/_static/figures/`; a
normal Sphinx / Read the Docs build does **not** need pyOpenMS, matplotlib, or
to regenerate the figures. For intentional figure regeneration, see the
maintainer-only `docs/figure_source/README.md`.

## License

See `LICENSE`, `NOTICE`, and `licenses/`.
