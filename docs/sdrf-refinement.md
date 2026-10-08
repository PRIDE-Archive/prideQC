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

## Audit and validation artifacts

An SDRF refinement directory can contain:

```text
refinement/
├── original.sdrf.tsv
├── refined.sdrf.tsv
├── sdrf-changes.tsv
├── sdrf-validation.json
├── cohort-refinement.json
├── llm-refinement-packet.json
├── llm-adjudication-request.json
└── manifest.json
```

The original SDRF is retained, proposed changes are recorded, and the refined output is validated with `sdrf-pipelines` before it is treated as a successful refinement.

## Optional adjudication

When deterministic evidence is informative but does not justify an automatic decision, prideQC can create a bounded evidence packet for human or local-LLM adjudication. The adjudication layer can accept, reject, or abstain on the candidate values already prepared by prideQC; it does not invent new measurements.
