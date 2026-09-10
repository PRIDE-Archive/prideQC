# prideQC

A unified Python workflow for **per-file mass spectrometry quality metrics,
technical SDRF evidence, and mzQC output**, integrating the useful analysis
contracts of [rawQC](https://github.com/bigbio/rawQC) and
[TechSDRF](https://github.com/bigbio/techsdrf).

Version 0.2 uses the `prideQC` distribution/repository name, the `prideqc` Python
package and the `prideqc` command. This integrated implementation has its own
API and does not claim complete behavioral compatibility with both upstream CLIs. See [migration](docs/migration.md)
for every feature retained, changed, or deferred, and
[validation](docs/validation.md) for what has actually been executed.

## Quick start with uv

Python 3.11–3.13; Python 3.12 is the development default.

```bash
cd prideQC
uv sync
uv run prideqc analyze data/*.mzML -o results/run-001
```

Refine an existing SDRF and process independent files concurrently:

```bash
uv run prideqc analyze data/*.mzML \
  --sdrf study.sdrf.tsv \
  --workers 2 \
  --continue-on-error \
  -o results/run-002
```

The destination must be new or empty. Each successful input receives its own
`.mzQC`, `.obo` vocabulary and `.summary.json`, plus shared `metrics.tsv`,
`annotations.tsv`, and `manifest.json`. Supplying an SDRF also produces
`refined.sdrf.tsv`, `sdrf-changes.tsv` and `sdrf-validation.json`. Keep each mzQC with its companion OBO.
Output filenames retain the entire input basename, e.g. `sample.mzML.mzQC`.
Use `--overwrite` when intentionally rerunning into an existing results
directory; it clears that directory's contents before analysis. Without this
flag, prideqc refuses to touch a non-empty destination.

The implementation was developed in an environment that blocks PyPI, so it
does **not** include a fabricated or stale `uv.lock`. On the first networked
checkout, `uv sync` resolves and creates it. Review and commit that lock, then
use `uv sync --locked`. The upstream lock is incompatible with its own source
metadata (it locks pyOpenMS 3.4 while the uploaded project requires >=3.5), so
it has deliberately not been copied.

## Dependencies and processing

The four direct runtime dependencies are **pyOpenMS, NumPy, sdrf-pipelines and
pridepy**. The latter two are maintained by the collaborating teams and own
SDRF validation and PRIDE acquisition. Imports are lazy: local QC without an
SDRF imports neither library; CLI help imports none of the scientific stack.

| Dependency | Responsibility |
| --- | --- |
| `pyopenms>=3.5,<4` | Native spectrum decoding and optional peak-type estimation |
| `numpy>=1.26,<3` | Compiled numerical reductions and exact quantiles |
| `sdrf-pipelines>=0.1.6,<0.2` | SDRF parsing for validation, templates and optional ontology checks |
| `pridepy>=0.0.16,<0.1` | PRIDE API access and file transfers through its current `Client` API |

To try a current OpenMS development build, configure uv/pip to use the
OpenMS package index (`https://pypi.openms.de/simple/pyopenms/`).

The base install uses sdrf-pipelines' structural/template validation. Add
`uv sync --extra ontology` for its ontology dependencies and explicitly request
`--validate-ontology` to use them. RunAssessor and Param-Medic remain excluded.
There is no custom PRIDE HTTP client or copied SDRF rule/ontology engine.
The small TSV editor preserves repeated columns, original values and evidence
conflicts; it does not decide SDRF conformance.

The CLI, JSON output, process orchestration and converter invocation use the
standard library. These four packages bring transitive dependencies, including
pandas and transport libraries; lazy imports reduce startup and worker overhead,
not installation size. Inspect `uv tree --no-dev` after resolving the environment.
No deployment-size or real-file speedup claim has been measured.

One pyOpenMS `MzMLFile.transform` pass decodes spectra and distributes them to
QC and optional evidence collectors. A bounded header read preserves exact
instrument CV terms. Peak arrays are discarded after each spectrum. Compact
scalar arrays retain exact quantiles: memory scales with **scan count plus the
largest spectrum**, not total peak count. This is not constant-memory QC.
Already ordered RT arrays avoid sorting; other arrays are sorted stably.

Defaults use one process per file sequentially. `--workers N` uses separate
processes, with small final summaries returned to the parent. CLI defaults
OpenMP/BLAS threads to one unless already configured. Worker count multiplies
memory and can saturate storage; select it using representative benchmarks.
For independent local files, `--workers 3` (or the number of files, bounded by
available CPU/RAM and storage bandwidth) analyzes them concurrently. The
dependency-free progress indicator prints an initial line and one flushed line
per completed/failed file; fail-fast runs also report files that were not
processed. Use `--no-progress` in non-interactive logs.

## Input formats and vendor RAW

Native analysis supports indexed or non-indexed `.mzML` and `.mzML.gz`.
Gzip input is expanded to a temporary disk file and cleaned up after analysis.

Thermo `.raw` and Bruker `.d` inputs are read directly when the installed
pyOpenMS build exposes its native `ThermoRawFile` or `BrukerTimsFile` reader.
Capability detection, rather than only
the version string, also works with nightly/development builds. Bruker
`.d.zip` files are extracted to a temporary directory when the native reader
accepts only a directory; the extraction is removed after analysis. Older
builds fall back to an explicitly selected external converter:

```bash
uv run prideqc analyze run.raw \
  --converter thermorawfileparser -o results/thermo

uv run prideqc analyze run.d \
  --converter msconvert -o results/bruker
```

If no native reader is available and no converter is selected, prideqc reports
the installed pyOpenMS version and suggests upgrading to a newer/current
development build or installing ThermoRawFileParser/msconvert. Development
builds are published at `https://pypi.openms.de/simple/pyopenms/`. The
dependency constraint remains `pyopenms>=3.5,<4` so environments can choose a
stable or development build; no converter is installed as a Python dependency.

`--converter-executable /path/to/program` selects a specific executable.
Use a wrapper executable for a converter that requires Mono or a launcher;
this argument is not a shell command. Conversion has a timeout, logs its output,
and must produce exactly one nonempty mzML. Converted data are retained under
`converted/` so input provenance remains resolvable. External conversion still
depends on the converter and platform. No centroiding filter is added by
prideqc; converter defaults still apply. Directory input and sidecar files
must be provided as the converter requires.

## SDRF annotations

Instrument identity comes from the mzML CV, without model substring guessing.
All instrument configurations, analyzers, ionization terms, serial numbers,
MS-level counts, scan windows, isolation widths, precursor charges, polarity,
FAIMS settings and MS2/MS3 fragmentation evidence are retained in the reports.
One unambiguous instrument and fully observed MS2 dissociation/energy values
can fill the corresponding SDRF columns.

Existing nonmissing values are preserved. Disagreements appear in
`sdrf-changes.tsv`; `--overwrite-sdrf-values` explicitly enables replacement.
Duplicate modification columns, repeated sample rows, unmatched files and
leading metadata comments are preserved. `not applicable` is treated as an
existing assertion, not a missing value.

Matching uses exact basenames, original source names embedded in mzML, or an
explicit map. No implicit stem matching merges unrelated files. For conversions
without source metadata:

```json
{"run.raw": "run.mzML"}
```

Save that as `file-map.json` and add `--file-map file-map.json`. Ambiguous input
basenames/aliases are rejected. Multiple instrument values remain in the report
without being collapsed into an invented single instrument.

DDA/DIA isolation-width heuristics are labeled **inferred**, require substantial
coverage, and abstain for mixed/intermediate/sparse data. They cannot reliably
separate DDA from PRM or narrow-window DIA. Only `--include-inferred` allows these
suggestions into the SDRF. The PRIDE accessions follow the current
[SDRF acquisition guidance](https://github.com/bigbio/proteomics-sample-metadata/blob/master/sdrf-proteomics/README.adoc).

Search tolerances are **unavailable** unless a future calibrated estimator is
provided; an isolation width or instrument category is not a search tolerance.
No enzyme, sample label, fixed modification or enrichment is invented.

`--diagnostics` adds a centroid MS2/MS3 reporter/oxonium screen for TMT-family,
iTRAQ-family and glycan signatures. This returns counts and thresholds,
not confirmed PTMs or plex/channel assignments, and never auto-fills SDRF.
`--estimate-peak-type` additionally runs OpenMS PeakTypeEstimator for scans with
more than ten peaks, on the same data pass.

Supplying `--sdrf` validates the input and refined output through
`sdrf_pipelines.sdrf.sdrf.read_sdrf(...).validate_sdrf(...)`. The default template
is `ms-proteomics`; use `--sdrf-template dia-acquisition` when appropriate.
Input findings are recorded and refinement can repair them. Any finding in the
refined output makes the workflow unsuccessful while retaining its per-file QC
and change reports. Validator/import failures are errors, never silent success.
Reports record the upstream version, template, ontology setting and findings.
The validity flag conservatively requires no upstream findings.

```bash
uv run prideqc check-sdrf study.sdrf.tsv --json validation.json
uv sync --extra ontology
uv run prideqc validate-sdrf study.sdrf.tsv --validate-ontology
```

`check-sdrf` and `validate-sdrf` are aliases. They now perform upstream template
validation, so an incomplete TSV that passed the old structural check may fail.
Ontology checks can need external services or upstream caches. Structural-only
validation is the default. Existing workflow conversions are available directly
through sdrf-pipelines, installed into the same uv environment:

```bash
uv run parse_sdrf convert-openms -s results/run-002/refined.sdrf.tsv
uv run parse_sdrf convert-msstats -s results/run-002/refined.sdrf.tsv -o msstats.csv
```

## PRIDE acquisition

Download exactly the filenames listed in an SDRF and run the same analysis path:

```bash
uv run prideqc analyze --accession PXD008644 \
  --sdrf study.sdrf.tsv --download-dir downloads/PXD008644 \
  --converter thermorawfileparser -o results/PXD008644
```

Omit the converter for an mzML-only selection. Repeated SDRF sample rows download
each file once. Add repeated `--file` arguments to restrict the selection;
`--file-map` can map SDRF cells to selected repository filenames. Download and
QC directories must be new/empty, separate and non-nested. Selections must be
exact basenames; folder-based vendor data and sidecars need separate acquisition.

Download without analysis, or select an explicit subset:

```bash
uv run prideqc fetch PXD008644 --sdrf study.sdrf.tsv -o downloads/fetch-001
uv run prideqc fetch PXD008644 --file run1.raw --file run2.raw -o downloads/fetch-002
```

The adapter uses `pridepy.download.client.Client.download_file_by_name`. It
passes protocol selection and checksum checking to pridepy, which owns its
transport/retry behavior. Checksums are requested by default;
`--no-checksum-check` disables that request. `--protocol` selects `ftp`, `aspera`,
`s3` or `globus`; available transports depend on the upstream installation.
A nonempty staged file is moved into place after the client returns. A download
failure stops acquisition before analysis, retains completed files, removes the
failed temporary file and writes `download-manifest.json`. This wrapper requires
a fresh destination; use pridepy's own CLI for advanced resumable acquisition.

For private projects, set `PRIDE_PASSWORD` in the environment and pass
`--username`. Use `--password-env NAME` to choose a different environment
variable. Credentials are not included in prideQC manifests. Broader metadata
search and non-PRIDE repository discovery remain available via the installed
`uv run pridepy --help`; prideQC's integrated selection accepts PRIDE `PXD` IDs.

## Python API

```python
from pathlib import Path

from prideqc.annotations import DiagnosticIonCollector
from prideqc.mzqc import MzQCWriter
from prideqc.pipeline import Analyzer, Workflow, WorkflowOptions

result = Analyzer().analyze(
    Path("run.mzML"), collectors=[DiagnosticIonCollector()]
)
MzQCWriter().write(result, Path("run.mzQC"))

manifest = Workflow(WorkflowOptions(workers=2)).run(
    ["a.mzML", "b.mzML"], "results/api", sdrf=Path("study.sdrf.tsv")
)
assert manifest["success"]

# Or acquire an explicit PRIDE subset before QC:
manifest = Workflow().run_project(
    "PXD008644", "results/project", download_directory="downloads/project",
    filenames=["run.mzML"], sdrf=Path("study.sdrf.tsv"),
)
```

`Analyzer` accepts a `SpectrumReader` implementation and an alternative
`QCMetricCalculator` or `TechnicalAnnotator`. Optional algorithms implement
`EvidenceCollector`; create a fresh collector for each run. Retain scalar
evidence or bounded samples, not incoming peak arrays. There are no entry-point
plugin frameworks. `SDRFValidator` and the injectable pridepy client allow
validation/acquisition to be tested independently of scan processing.

## Metrics, semantics and provenance

The metric catalog is in [docs/metrics.tsv](docs/metrics.tsv). Quantile vectors
and tables are serialized as complete values instead of repeated scalar CV
terms. Summary JSON/TSV explicitly uses `null` for undefined values; unavailable
top-level metrics are omitted from mzQC. Empty counts remain zero.

Some rawQC numerical definitions intentionally change: exact TIC boundary
integration, invalid peak exclusion, and removal of MS2-TIC substitution for
missing precursor intensity. Negative/unknown charges are excluded from mean
and range but remain in the fraction denominator. ID-free proxies do not
borrow ID-based accessions. [Migration notes](docs/migration.md) describe this.

mzQC contains the analyzed input URI, prideqc/OpenMS versions, CV declarations
and instrument terms. Custom terms have deterministic `QCPRIDE` accessions and
a per-file OBO. These local URIs identify the original output location; when
moving or publishing results, update the local CV/input URIs appropriately.
Schema validation and full ontology semantic validation are distinct gates.
This initial implementation does not claim complete semantic certification.

## Development

```bash
uv sync --group dev
uv run python -m unittest discover -s tests -v
uv run ruff check .
uv run ruff format --check .
uv run mypy
uv build
uv tree --no-dev
```

The CI matrix checks Python 3.11–3.13 with minimum/current versions of
pyOpenMS, sdrf-pipelines and pridepy. `PRIDE_QC_REQUIRE_INTEGRATION=1` prevents missing integration
dependencies from silently turning a CI run green. The test suite uses the
standard library's unittest framework. See [benchmarks](benchmarks/README.md)
for repeatable measurements; no real-file speedup has yet been measured.

Exit codes: **0** complete workflow or valid SDRF; **1** file/output failure,
incomplete batch or SDRF findings; **2** invalid arguments/setup or acquisition failure. `--continue-on-error` still returns 1 on
partial failure. Successful per-file outputs remain available after failures.
Already running parallel workers may finish when fail-fast stops accepting
further results; the manifest identifies all inputs without reported results.

## Attribution

Derived contracts and selected algorithms: rawQC (MIT, BigBio 2025), TechSDRF
(Apache-2.0, PRIDE Team). Changes and retained notices are in [NOTICE](NOTICE)
and [licenses](licenses/). The mzQC schema fixture is from HUPO-PSI, CC BY 4.0.
