# SDRF refinement

The second major role of prideQC is to use the QC evidence computed from raw data to improve **technical SDRF metadata**.

```{figure} _static/figures/figure5_cohort_sdrf.svg
:alt: Synthetic run grouping, supported-tolerance maximum-versus-mean comparison, PTM evidence gate matrix, and validated SDRF write-back outcomes.
:width: 100%
:class: doc-figure

**Synthetic cohort-refinement illustration.** Compatible runs provide shared tolerance proposals only after coverage and unit checks; using the supported maximum covers the noisiest run. Mass-shift proposals additionally pass identity, category and independent semantic-evidence gates. Ineligible cases remain unchanged or undergo review.
```

## Why refinement is a separate stage

A single run can be noisy or atypical. A technical SDRF value often describes a larger experimental group or an entire study. prideQC therefore separates:

- **run-level measurement**, where QC evidence is computed;
- **cohort-level synthesis**, where compatible evidence is combined;
- **metadata write-back**, where policy and validation determine whether a proposal is safe.

## Typical refinement targets

Current refinement logic can support technical fields such as:

- precursor mass tolerance;
- fragment mass tolerance;
- modification parameters;
- annotation provenance.

## Cohort synthesis

For each relevant set of runs, prideQC asks whether the evidence is sufficiently complete and compatible to support a shared recommendation.

For tolerance recommendations, the study-level value is intentionally conservative: the refinement layer chooses a common setting that still covers the supported runs rather than silently making every run use a different tolerance.

## Modification proposals

Recurrent mass-shift evidence can contribute to a modification proposal, but mass compatibility alone is insufficient. prideQC can require independent PTM semantic/study evidence before a candidate is eligible for standard SDRF modification write-back.

## Preserve, propose, abstain

A useful refinement system must be able to say **no**.

prideQC therefore distinguishes between:

- evidence strong enough to support a proposal;
- existing metadata that should be preserved;
- ambiguous or incomplete cases that should remain unchanged or go to review.

## Measured run-level QC proposals (standards-pending)

The three measurements below are computed directly from each raw data file. They
are **observed QC measurements**, not experimental search parameters and not
cohort-level tolerance recommendations.

| prideQC metric | PSI-MS term | Unit | Proposed, **not yet approved**, SDRF header |
|---|---|---|---|
| `NumberOfSpectra_MS1` | `MS:4000059` | count (`UO:0000189`) | `comment[number of MS1 spectra]` |
| `NumberOfSpectra_MS2` | `MS:4000060` | count (`UO:0000189`) | `comment[number of MS2 spectra]` |
| `ChromatographyDuration` | `MS:4000067` | seconds (`UO:0000010`) | `comment[MS run duration]` |

The counts include empty spectra; zero is a valid count. Duration is the
maximum minus minimum **finite observed spectrum retention time**, over all
measured MS levels. A one-RT scan can have zero duration, while a file without
any finite spectrum RT has **no duration measurement**. Duration is not
necessarily the programmed LC gradient time. Do not substitute the SDRF
method-setting fields `comment[ms min rt]` or `comment[ms max rt]`.

The proposed headers are absent from the current SDRF-Proteomics `TERMS.tsv`
registry and the `ms-proteomics` v1.1.0 template. Accordingly, prideQC does
**not** write them into the refined SDRF. During the existing `--refine-sdrf-qc`
workflow, prideQC instead produces `sdrf-qc-proposals.tsv` and
`sdrf-qc-proposals.json` as **read-only, per-row proposals**.

For example, a matching `sample.raw` with 0 MS1 scans, 48 MS2 scans, and
120 seconds of observed RT evidence yields proposal records such as:

```text
metric_key               proposed_column                    value   unit
NumberOfSpectra_MS1      comment[number of MS1 spectra]      0       count
NumberOfSpectra_MS2      comment[number of MS2 spectra]      48      count
ChromatographyDuration   comment[chromatography duration]    120 s   second
```

These are illustrative values, not an SDRF output header/value grammar approved
by the standard. The actual TSV contains additional linkage, ontology, status,
existing-value, and method/provenance columns. Each matching SDRF row gets its
own proposal even when multiple rows represent the same data file. Counts and
RT spans are never averaged across files. Unmatched files and missing/invalid
measurements remain explicitly unresolved rather than being replaced with zero.

Existing values are never overwritten by the proposal exporter. Proposed
headers already present in the input are inspected for matches, conflicts,
missing cells, or duplicate columns, but the original is unchanged. The
standard refined SDRF continues through the existing `sdrf-pipelines`
validation step; neither ontology nor schema checks are relaxed. Once the
SDRF term registry and value grammar are accepted upstream, a separate
reviewed change can enable guarded opt-in write-back via the existing
`SDRFDocument.annotate()` mechanism.

Putative PTM candidates are a separate scientific-evidence lane: an
unregistered miscellaneous-parameter representation needs upstream review,
but it must not be silently reinterpreted as a confirmed
`comment[modification parameters]` assertion.

## Audit and validation artifacts

An SDRF refinement directory can contain:

```text
refinement/
├── original.sdrf.tsv
├── refined.sdrf.tsv
├── sdrf-changes.tsv
├── sdrf-validation.json
├── cohort-refinement.json
├── sdrf-qc-proposals.tsv
├── sdrf-qc-proposals.json
├── llm-refinement-packet.json
├── llm-adjudication-request.json
└── manifest.json
```

The original SDRF is retained, proposed changes are recorded, and the refined output is validated with `sdrf-pipelines` before it is treated as a successful refinement.

## Optional adjudication

When deterministic evidence is informative but does not justify an automatic decision, prideQC can create a bounded evidence packet for human or local-LLM adjudication. The adjudication layer can accept, reject, or abstain on the candidate values already prepared by prideQC; it does not invent new measurements.
