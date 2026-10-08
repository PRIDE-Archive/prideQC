# CLI reference

## `prideqc analyze`

Compute QC metrics and evidence from local or PRIDE-hosted mass-spectrometry files.

### Core arguments

- `files...` — local input files.
- `--output-dir PATH` — required output directory.
- `--sdrf PATH` — optional SDRF used for mapping/refinement.
- `--accession PXD...` — download selected PRIDE files before analysis.
- `--download-dir PATH` — separate download location.
- `--workers N` — independent file workers.
- `--overwrite` — clear an existing output directory.

### Evidence switches

- `--diagnostics`
- `--estimate-mass-error`
- `--estimate-mass-shifts`
- `--estimate-peak-type`

### SDRF refinement switches

- `--include-inferred`
- `--overwrite-sdrf-values`
- `--refine-sdrf-qc`
- `--ptm-study-evidence PATH`
- `--project-accession PXD...`

## `prideqc refine-sdrf-qc`

Build accession/study-level refinement outputs from previously generated `*.summary.json` artifacts.

### Core arguments

- `--results-root PATH`
- `--sdrf PATH`
- `--output-dir PATH`
- `--file-map PATH`
- `--ptm-study-evidence PATH`
- `--project-accession PXD...`

## `prideqc check-sdrf`

Validate an SDRF via `sdrf-pipelines`.

### Core arguments

- `file`
- `--json PATH`
- `--sdrf-template TEMPLATE`
- `--validate-ontology`

## `prideqc llm`

Manage the optional local adjudication runtime.

### Subcommands

- `prideqc llm setup`
- `prideqc llm status`
- `prideqc llm adjudicate`

## Version and help

```bash
prideqc --version
prideqc --help
prideqc analyze --help
prideqc refine-sdrf-qc --help
```
