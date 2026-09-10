# Integration decisions

The supplied snapshots are rawQC 0.1.1 and TechSDRF 1.0.0. The source archives
are the implementation baseline; no upstream commit identifier was supplied.
Archive SHA-256 values are recorded in `source-manifest.json`. The original
archives remain separate; this repository is a new integrated implementation.

## Naming and dependency update (0.2)

The repository/distribution is `prideQC`, the import is `prideqc`, and the command
is `prideqc`. The initial artifact used `pride-qc` / `pride_qc`; update imports,
scripts and downstream `pride_qc_version` manifest consumers to `prideqc_version`.
No legacy import alias is installed. Metric keys and QCPRIDE accessions remain
stable. Historical benchmark JSON retains the original software labels.

This update reverses the initial decision to defer SDRF validation and PRIDE
acquisition: both libraries are maintained by the collaborating teams. It keeps
the evidence-aware editor because duplicate columns, cell preservation and
annotation conflict handling are specific to this workflow. It does not copy
SDRF validation rules or implement another PRIDE HTTP/FTP transport.

## Functional scope

| Upstream feature | prideQC implementation | Status |
| --- | --- | --- |
| Local mzML QC | Shared pyOpenMS streaming pass, exact scalar summaries | Implemented |
| Scan/peak counts, RT/rates, TIC, density, charges, precursor distributions | Typed calculator; expanded MS3+ coverage; vector/table output | Implemented, selected definitions corrected |
| Polarity, FAIMS, chromatogram counts/time ranges | Shared scan/chromatogram collectors | Implemented |
| Annotated and estimated profile/centroid type | Annotated counts always; estimator opt-in in same pass | Implemented |
| Instrument identity, analyzer and ionization | Preserve instrumentConfiguration CV terms, all configurations | Implemented; no fragile name lookup |
| Analyzer resolution, scan-filter-specific detector type | Not reliably preserved by this header adapter | Deferred; no fabricated value |
| MS2/MS3 dissociation and collision energy | Per-level evidence, units, missing-observation coverage | Implemented |
| Scan windows, isolation width, charge range | Metrics and technical summaries | Implemented; width labeled Th/m/z |
| DDA/DIA classification | Width heuristic with minimum coverage, explicit inference, abstention | Implemented conservatively |
| SDRF read/write/comparison/report | Standard-library TSV, exact file mapping, duplicate columns, conflict report | Implemented |
| SDRF template/ontology validation | Delegate input/output checks to sdrf-pipelines; templates configurable, ontology opt-in | Restored in 0.2 |
| Per-file mzQC | Dedicated JSON and OBO per input, provenance, no pymzqc runtime | Implemented; schema gate supplied, full semantic audit pending |
| TMT/iTRAQ reporter and glycan diagnostic evidence | Optional centroid-only screen with counts/thresholds | Implemented as candidate evidence |
| Specific label plex/channel assignment | No assignment from ambiguous diagnostic evidence | Deferred |
| Phospho neutral loss, mass-shift PTM pairing and null-model significance | Not carried into the lightweight core | Deferred pending validated native algorithms |
| Param-Medic/RunAssessor/Crux empirical tolerances | No stale wrapper dependencies; absent estimate is explicit | Deferred pending calibrated estimator |
| Instrument-default search tolerance guesses | Do not equate nominal class/resolution with measured tolerance | Removed |
| Thermo RAW and Bruker `.d`/`.d.zip` input | Feature-detected native pyOpenMS readers; temporary extraction for `.d.zip`; explicit converter fallback | Implemented; real vendor files still require platform fixtures |
| PRIDE downloads | Exact SDRF/explicit subset through pridepy Client, then shared QC workflow | Restored in 0.2 |
| Broader repository discovery, search and workflow conversions | Use the installed pridepy and parse_sdrf CLIs | Reuse upstream functionality directly |
| Heatmaps, notebooks and online demo downloads | Consume metrics.tsv with separate visualization tooling | Outside core |
| Accepted identification IDs/ID-based metrics | Current mission uses raw, ID-free data | Not exposed in this API |

`rawqc-metric-migration.tsv` maps all 90 registered upstream metric keys to new
keys/tables or an explicit removal. This mapping concerns data products rather
than compatibility of the old Python function names or command-line options.

## Dependency decisions

| Dependency | Decision and reason |
| --- | --- |
| pyopenms >=3.5,<4 | Kept for mzML/native decoding and optional peak-type estimation; feature-detected Thermo/Bruker readers are used when exposed |
| numpy >=1.26,<3 | Kept for compiled reductions/quantiles; declared directly because imported directly |
| RunAssessor, Param-Medic | Removed; no vendored copy or subprocess invocation in core |
| SciPy, pyteomics | No core functionality needs an additional numerical/parser stack |
| pymzqc, pronto/ontology graph stack | Replace runtime serialization with stdlib JSON, local definitions and separate conformance tests |
| sdrf-pipelines >=0.1.6,<0.2 | Required; owns SDRF templates and validation, ontology extra optional |
| pridepy >=0.0.16,<0.1 | Required; current Client API owns PRIDE downloads; legacy Files API removed upstream |
| pandas, click, requests, transport and plotting libraries | No extra direct imports; transitive dependencies of the maintained packages remain |
| External vendor converters | Explicitly selected executable, not a Python package dependency |
| jsonschema, ruff, mypy | Development group only |
| pytest | Standard-library unittest is sufficient for this repository |

Do not conceal transitive dependencies by omitting imports from pyproject.toml.
`uv tree --no-dev` is the source of truth after dependency resolution.

## Scientific changes requiring attention

- **Precursor intensity:** only recorded positive first-precursor intensity is
  used. The OpenMS zero/unset default is missing. MS2 fragment TIC is a different
  measurement and is not substituted. Compare only positive recorded cases when
  assessing parity with rawQC's precursor statistics.
- **TIC quartile integration:** split each piecewise-linear trapezoid at the RT
  quartile boundary and integrate it exactly. For a ramp from 0 to 20 over 10
  seconds, the four areas are 6.25, 18.75, 31.25 and 43.75. Interpolating
  cumulative area linearly gives 25 each and has the wrong boundary behavior.
- **Missing and invalid data:** finite paired RT/TIC values define temporal
  calculations. All scans still count. Invalid/negative peak values are excluded
  with counts and warnings. Undefined values remain null, not zero.
- **Charge:** the first precursor supplies per-scan distributions. Nonpositive
  charges are unknown; they remain in fraction denominators but not known-charge
  means or extrema. Multiple-precursor scans are counted and reported.
- **Thresholds:** signal jumps/falls retain rawQC's inclusive >=10 / <=0.1
  convention, explicitly defined with local accessions. Ratios from zero are
  undefined; a fall to zero from a positive value counts.
- **ID-free proxies:** IQR rates, middle/shortest-half TIC and precursor dynamic
  range do not claim ID-based PSI-MS terms.
- **CV mappings:** quantile tuples have one metric value. Conservative standard
  mappings coexist with a defined local CV. The TechSDRF source had contradictory
  instrument tables, so exact input CV terms replace substring lookup. DDA uses
  PRIDE:0000627 and DIA PRIDE:0000450 from current SDRF guidance.
- **Acquisition inference:** widths alone cannot prove acquisition type. A
  minimum of 100 valid MS2 windows, >=90% coverage and >=90% consensus is needed
  even to emit a candidate. Intermediate/mixed/sparse runs abstain.
- **No dataset-wide majority override:** each SDRF row receives evidence from
  its own matched file. Heterogeneous files retain distinct technical details.

## Extending the package

`SpectrumReader` and `SpectrumSink` decouple file decoding. `Analyzer` owns the
per-file lifecycle and shares each event with `RunSummary` and optional
`EvidenceCollector` instances. `QCMetricCalculator` performs pure calculations;
`TechnicalAnnotator` labels evidence. `MzQCWriter` and `SDRFDocument` own output
formats. `SDRFPipelinesValidator` delegates conformance to upstream rules through
the `SDRFValidator` protocol. `PrideRepository` delegates transfers to pridepy.
`Workflow` owns acquisition, files, conversion, processes and failure manifests.

A replacement mass-error estimator should supply an `EvidenceCollector` with a
documented calibration domain, bounded retained state, error units, uncertainty,
support counts and an abstention rule. Validate it against independently known
mass errors and benchmark it before allowing its suggestions into SDRF. The
same interface can support identification-aware modules later without adding
them to the default dependency graph.

## Release work remaining

1. Resolve with network access, inspect transitive dependencies, commit uv.lock.
2. Run the supplied minimum/current pyOpenMS CI matrix, schema validation,
   formatter and type checker; resolve any binding/platform differences. Run the
   sdrf-pipelines validation and installed pridepy Client signature tests as well.
3. Audit canonical CV labels/value shapes against an explicitly pinned ontology
   snapshot and perform full mzQC semantic validation in addition to JSON Schema.
4. Run representative public and collaborator raw/mzML files, including vendor
   conversion, mixed configurations, DIA and MS3. Review output with SDRF owners.
5. Measure real-file throughput/RSS against both upstream workflows, then choose
   worker defaults for deployment. The included synthetic data are not a speedup claim.

These are release-validation steps, not silently completed claims.
