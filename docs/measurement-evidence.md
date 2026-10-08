# Measurement precision & mass-shift evidence

Two prideQC analyses are especially useful for turning spectral patterns into technical QC evidence: **mass-error estimation** and **recurrent neutral-mass-shift scouting**.

Neither analysis starts from peptide identifications. Both use repeated patterns in the measured spectra themselves.

## Mass-error estimation

```{figure} _static/figures/figure3_mass_error.svg
:alt: Synthetic repeated-precursor and fragment mass errors shown as pairwise distributions, robust precision estimation, and gated six-sigma tolerance suggestions.
:width: 100%
:class: doc-figure

**Synthetic identification-free precision example.** Repeated observations yield pairwise errors; robust median/MAD estimates are converted to single-measurement precision by dividing by √2. A supported six-sigma suggestion depends on evidence sufficiency and resolution regime. These are not recovered historical search tolerances.
```

### Concept

If the same underlying ion is measured repeatedly, independent measurement noise causes the observed masses to move slightly around the underlying value. Differences between repeated observations therefore provide a way to estimate **measurement precision**.

prideQC uses repeat-observation evidence to estimate robust precursor and fragment error distributions. The key ideas are:

1. identify sufficiently comparable repeated observations;
2. calculate pairwise mass differences;
3. estimate a robust center and spread while limiting the influence of outliers/background matches;
4. convert pairwise spread to an approximate single-measurement precision;
5. emit a tolerance suggestion only when support and confidence requirements are met.

The resulting fields can include:

- estimated precursor mass error in ppm / Da;
- estimated fragment mass error in ppm / Da;
- suggested precursor search tolerance;
- suggested fragment search tolerance;
- support counts, robust inlier diagnostics, and resolution-regime information.

### OpenMS implementation

prideQC uses OpenMS and pyOpenMS for mass-spectrometry data handling and several QC algorithms. When available, mass-error estimation is performed by the native OpenMS `IDFreeMassErrorEstimator` and adapted into prideQC's structured evidence model.

This implementation detail is useful for reproducibility, but the user-facing interpretation remains the same: the output describes **measurement precision evidence**, not a reconstruction of the original historical search parameters.

- OpenMS source: <https://github.com/OpenMS/OpenMS>
- pyOpenMS documentation: <https://pyopenms.readthedocs.io/>

## Recurrent neutral-mass-shift scouting

```{figure} _static/figures/figure4_mass_shift.svg
:alt: Related synthetic MS2 spectra with shared-bin candidate retrieval, unchanged and shifted fragment matches, recurrent neutral-mass clusters, and tight modification-candidate mass windows.
:width: 100%
:class: doc-figure

**Synthetic mass-shift scouting example.** Related MS2 pairs support recurrent neutral-mass differences. Clustering tolerance and chemical mass-annotation tolerance are distinct; OpenMS/UniMod compatibility neither identifies a modification nor localizes it to a residue.
```

### Concept

Two MS2 spectra can have a similar fragment pattern while their precursors differ by a reproducible neutral mass. Repeated occurrences of the same shift create a **mass-shift family**.

prideQC:

1. searches for related spectra with sufficient fragment-pattern overlap;
2. calculates their neutral-mass differences;
3. clusters recurrent shifts;
4. records support, cluster width, spectral similarity, and other diagnostics;
5. compares the observed cluster mass with compatible OpenMS / UniMod entries.

### Mass compatibility is not identity

A mass match does **not** prove a PTM identity. Different modifications can have similar or identical masses, and some recurrent shifts may be technical artifacts rather than biological modifications.

For that reason, prideQC keeps candidate ambiguity visible and requires additional study/semantic evidence before a mass-shift candidate can become a standard SDRF modification proposal.

## Why both analyses matter for metadata refinement

Mass-error estimation can support cohort-wide precursor/fragment tolerance recommendations.

Mass-shift scouting can surface modification-related evidence that may be missing from the deposited SDRF.

In both cases, the QC result remains evidence first; [SDRF refinement](sdrf-refinement.md) applies a separate cohort-level policy before changing metadata.
