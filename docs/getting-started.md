# Getting started

## Installation

### From PyPI

```bash
python -m pip install prideQC
```

or with `uv`:

```bash
uv tool install prideQC
```

### Development install

```bash
git clone https://github.com/PRIDE-Archive/prideQC.git
cd prideQC
uv sync --active --group dev
```

## Analyze raw data

A typical analysis enables the main evidence collectors:

```bash
prideqc analyze sample1.mzML sample2.mzML \
  --output-dir results \
  --estimate-mass-error \
  --estimate-mass-shifts \
  --diagnostics
```

Useful options include:

- `--estimate-mass-error` — estimate precursor and fragment measurement precision.
- `--estimate-mass-shifts` — detect recurrent related-spectrum neutral mass shifts.
- `--diagnostics` — screen MS2/MS3 spectra for diagnostic-ion evidence.
- `--estimate-peak-type` — estimate centroid/profile representation where needed.
- `--workers N` — process independent files in separate worker processes.

## Example output tree

For two input mzML files with mass-error and mass-shift estimation enabled, an analysis directory can look like:

```text
results/
├── manifest.json
├── metrics.tsv
├── annotations.tsv
├── mass-shifts.tsv
├── sample1.mzML.mzQC
├── sample1.mzML.summary.json
├── sample1.mzML.mass-shifts.tsv
├── sample2.mzML.mzQC
├── sample2.mzML.summary.json
└── sample2.mzML.mass-shifts.tsv
```

If an SDRF is supplied and cohort refinement is enabled, additional artifacts can include:

```text
results/
├── original.sdrf.tsv
├── refined.sdrf.tsv
├── sdrf-validation.json
├── sdrf-changes.tsv
├── cohort-refinement.json
├── llm-refinement-packet.json
├── llm-adjudication-request.json
└── manifest.json
```

The exact set depends on the enabled analysis/refinement options.

## Refine an SDRF from existing QC results

If per-file summaries already exist, refinement can be run separately:

```bash
prideqc refine-sdrf-qc \
  --results-root results \
  --sdrf study.sdrf.tsv \
  --output-dir refinement
```

Optional inputs include:

- `--file-map` when SDRF filenames differ from analyzed filenames;
- `--ptm-study-evidence` when independent PTM semantic evidence should gate modification write-back;
- `--project-accession` when the accession should be supplied explicitly for provenance.

## Validate an SDRF

```bash
prideqc check-sdrf refinement/refined.sdrf.tsv
```

With ontology validation enabled:

```bash
uv sync --active --extra ontology
prideqc check-sdrf refinement/refined.sdrf.tsv --validate-ontology
```

## Next steps

- Read [QC workflow](qc-workflow.md) for how the raw-data evidence is computed.
- Read [Measurement precision & mass shifts](measurement-evidence.md) for the two main spectral-evidence algorithms.
- Read [Reporting](reporting.md) for pmultiqc / MultiQC integration.
- Read [SDRF refinement](sdrf-refinement.md) for the study-level write-back model.
