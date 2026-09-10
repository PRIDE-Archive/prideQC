# Validation record

Run date: 2026-09-09. Environment: Linux x86-64, CPython 3.12.14. The source
archives are recorded with SHA-256 in `source-manifest.json`.

## Executed successfully

- **88 tests discovered: 81 passed, 7 skipped.** The executed tests cover exact
  TIC integration, irregular/duplicate/invalid retention times, level-specific
  scan rates, fixed-window frequency, sample CV, empty/invalid peaks, zero
  intensity, charge denominators, MS3, FAIMS, chromatograms, ID-free proxy
  formulas, releasing peak arrays, single-pass collector dispatch, evidence
  thresholds, SDRF duplication/matching/conflicts, atomic output, provenance,
  enum conversion, fail-fast/partial failures, and CLI behavior. The 0.2 update
  adds delegation tests for SDRF templates/ontology flags, input and refined
  validation reports, failed validator behavior, PRIDE selection/deduplication,
  the current Client import path, download failure cleanup, protocol/checksum
  forwarding, credential handling, acquisition-to-QC orchestration, and
  feature-detected native Thermo/Bruker reader paths including `.d.zip`
  extraction cleanup, explicit output-directory overwrite, and flushed progress
  reporting. These local adapter tests use doubles; they do not claim
  live service validation.
- Python compilation for source and tests.
- Source-distribution and wheel build using uv with the already installed
  setuptools build backend, in offline mode without build isolation.
- Wheel installed into a clean project venv with dependencies deliberately
  excluded; the installed `prideqc --version` returns `prideQC 0.2.0`, and
  `import prideqc` resolves version 0.2.0. Wheel metadata includes all four direct
  dependencies and the optional ontology extra.
  This verifies packaging/entry point, **not native processing**.
- A catalog covers 81 default MS1/MS2 summary keys; further levels expand the
  same definitions. Missing measurements are explicit. Every one of the 90
  registered rawQC source metric keys has a documented migration target/status.
- Three independent synthetic runs measured the real reduction/statistics code
  without file decoding. Results are retained as JSON, with Python/NumPy/platform
  metadata. These historical 0.1 measurements predate the maintained-library
  integrations and remain illustrative, not a real-file performance claim.

| Synthetic input | Wall time | Peak RSS |
| --- | ---: | ---: |
| 10,000 scans × 1,000 peaks | 0.171 s | 27.6 MiB |
| 10,000 scans × 10,000 peaks | 0.437 s | 27.5 MiB |
| 100,000 scans × 1,000 peaks | 1.728 s | 35.9 MiB |

RSS is the process lifetime high-water mark, not an allocation delta. Small
differences across fresh processes are noise. The synthetic stream reuses peak
buffers and specifically checks that prideqc does not retain total peak data.
It does not represent mzML decoding, raw conversion, or storage throughput.

## Not executed here

PyPI returned HTTP 403 under this environment's network restrictions. No
pyOpenMS, sdrf-pipelines, pridepy, jsonschema, Ruff or mypy distribution was
present locally. An installation attempt for the restored dependencies was
blocked by the same registry restriction. Offline uv
resolution also failed because required distributions were not cached.

Consequently, the following are **pending**, not verified successes:

- Four native integration tests: indexed/non-indexed mzML, gzip input,
  serial/parallel equivalence and corrupt-file handling.
- One official JSON Schema test with URI/date-time format validation.
- Two installed-library checks: actual sdrf-pipelines parser/validator findings
  and an autospecced pridepy Client signature check without network access.
- Live PRIDE transfers, private-project access, protocol fallback and checksums.
- Ontology-enabled SDRF validation with the upstream optional dependencies.
- Ruff lint/format and mypy execution.
- A fresh resolved and committed uv.lock.
- Full ontology-semantic mzQC conformance and cross-tool numeric parity on
  real data.
- Actual Thermo/Bruker vendor-file decoding, native-reader behavior across
  pyOpenMS builds, and ThermoRawFileParser/msconvert invocation remain
  platform-dependent and require representative vendor fixtures.
- Real-file performance comparison against rawQC/TechSDRF and process scaling.

The CI matrix installs required dependencies and sets
`PRIDE_QC_REQUIRE_INTEGRATION=1`, which makes absent pyOpenMS, jsonschema,
sdrf-pipelines or pridepy a test failure rather than silently accepting skips. It targets Python 3.11–3.13 and
minimum/current dependency combinations: pyOpenMS 3.5.0, sdrf-pipelines 0.1.6
and pridepy 0.0.16 versus currently resolvable versions within declared bounds. CI has been
provided but has not been run on a hosted service in this task.

No upstream test suite is claimed to pass unchanged: the public APIs differ,
and selected numerical definitions intentionally change. New regression tests
are derived from the audited upstream edge cases and independent analytical
oracles. The complete upstream metric mapping supports further parity testing.

## Schema provenance

`tests/data/mzqc_schema.json` is the unmodified mzQC 1.0.0 JSON schema present in
the supplied rawQC archive. Its `$id` points to the HUPO-PSI source. It is used
as an offline fixture; default QC downloads no schema or ontology. Optional
SDRF ontology validation uses sdrf-pipelines and may access its external services.

- [HUPO-PSI mzQC specification and schema](https://github.com/HUPO-PSI/mzQC)
- [PSI-MS controlled vocabulary](https://github.com/HUPO-PSI/psi-ms-CV)
- [pyOpenMS MzMLFile streaming interface](https://pyopenms.readthedocs.io/en/release-3.5.0/apidocs/_autosummary/pyopenms/pyopenms.MzMLFile.html)
- [OpenMS pyOpenMS development package index](https://pypi.openms.de/simple/pyopenms/)
- [SDRF-Proteomics technical annotation guidance](https://github.com/bigbio/proteomics-sample-metadata/blob/master/sdrf-proteomics/README.adoc)

Schema validity checks structure and JSON value types. It does not establish
that an accession, canonical label, unit, category and mathematical definition
agree. The latter remains a separate release gate, explicitly documented in
`migration.md`.

## Integration API sources

The supplied TechSDRF reader/refiner uses
`sdrf_pipelines.sdrf.sdrf.read_sdrf(...).validate_sdrf(...)`; its PRIDE downloader
uses the former `pridepy.files.files.Files`. Current upstream documentation
records the replacement of `Files` with `pridepy.download.client.Client` in
0.0.16, with the public download methods retained. We target that released API,
not the unreleased 0.0.17 version visible in the upstream master pyproject.
These documents were checked on 2026-09-09:

- [pridepy Python API and download options](https://github.com/PRIDE-Archive/pridepy/blob/master/docs/usage.md)
- [Published pridepy distribution](https://pypi.org/project/pridepy/)
- [sdrf-pipelines validation and conversion commands](https://github.com/bigbio/sdrf-pipelines/blob/main/COMMANDS.md)
- [sdrf-pipelines base/ontology installation and published release](https://pypi.org/project/sdrf-pipelines/)

The blocked package install prevents execution against these releases here.
The installed-library CI gates are required before releasing 0.2.

## Mypy and NumPy stub compatibility

The package supports Python 3.11 through 3.13 at runtime. The quality job runs
mypy under Python 3.12 because current NumPy 2.3 type stubs use Python 3.12's
`type` statement. `tool.mypy.python_version = "3.12"` selects the syntax target
used to parse those stubs; it does not raise the package's runtime floor or
permit Python 3.12-only syntax in the source. Ruff continues to target `py311`,
and the source/test grammar checks remain Python 3.11-compatible.

## FAIMS fixture correction

The user's pyOpenMS 3.5.0 diagnostic showed no `MS:1001581` CV terms in the
original writer-generated test file. Reloaded spectra had the missing-value
sentinel (`-1.0`) and unspecified drift-time unit (`0`). Therefore the failing
`[-45.0]` assertion was being applied to a fixture that did not contain the
expected voltage. This observation does not establish a reader defect or which
later OpenMS releases change the writer's behavior.

The integration tests now copy explicit synthetic input fixtures from
`tests/data/faims-nonindexed.mzML` and `tests/data/faims-indexed.mzML`, with
`MS:1001581 = -45 V` on the two MS1 scans. No OpenMS writer generates these
inputs at test time. The four scans retain the original RTs (0, 10, 20, 30 s),
two peaks per scan (m/z 100 and 200; intensities 10 and 20), polarity and MS2
precursor/fragmentation evidence. These files are artificial test data, not an
experimental acquisition.

The original FAIMS metric assertion remains. An additional native-reader
assertion checks values and FAIMS units before the prideQC adapter is called;
the gzip test also checks the negative voltage. Three dependency-free checks
verify the explicit voltage/unit and peak bytes, index offsets/checksum, and
equivalent contents in both representations. All three passed locally. The
native integration rerun still requires the user's environment or CI.

The indexed fixture follows the byte-offset and checksum definition in the
[HUPO-PSI indexed mzML schema](https://github.com/HUPO-PSI/mzML/blob/81e0145f65dd7abf56a9bea51bfdf66ed7767905/schema/schema_1.1/mzML1.1.1_idx.xsd).
The two data files are included in the source distribution via `MANIFEST.in`.
No production reader, metric definition, dependency bound or lockfile changes
are needed for this correction.
