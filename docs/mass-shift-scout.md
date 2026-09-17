# Identification-free recurrent mass-shift scout (v22.1)

The mass-shift scout is an opt-in QC lane for finding recurring precursor
neutral-mass differences between **fragment-related MS2 spectra**. It is intended
to surface modification-like chemistry quickly without making peptide, protein,
or site-localization claims.

## Scientific boundary

A reported mass shift is a hypothesis. A UniMod match means only that the
observed recurrent delta is mass-compatible with one or more modification
records exposed by the pinned OpenMS runtime. The scout does not identify a
sequence, localize a residue, or prove that a PTM is present.

Exact isotope-selection and common positive-mode adduct substitutions are
classified before UniMod annotation. v22.1 additionally recognizes broader
isotope-adjacent satellite families so they do not dominate the QC-facing PTM
list. Those suppressed structures remain counted in diagnostics.

The default UniMod display is deliberately narrower than the diagnostic mass
match. Decoy records and amino-acid substitutions are hidden from the normal QC
candidate list, while biological PTM-like, sample-preparation/chemical, and
other mass-compatible chemistry can remain visible. The broad bounded candidate
set remains embedded in the structured diagnostic evidence for each reported
cluster.

Because evidence comes from pairs of related spectra with different precursor
masses, the scout is most sensitive to variable or co-occurring modified and
reference-like forms. A fixed modification carried by every corresponding
peptide can be invisible when no related reference form exists in the run.

## Spectrum normalization

v22.1 preserves native centroid MS2 exactly as before and adds profile support:

```text
MS2 spectrum
  -> native centroid -------------------------------> fingerprint
  -> native profile -> OpenMS PeakPickerHiRes -----> fingerprint
  -> unknown + OpenMS type estimate
       -> estimated centroid -----------------------> fingerprint
       -> estimated profile -> PeakPickerHiRes ----> fingerprint
       -> unresolved -------------------------------> abstain
```

Profile centroiding is ephemeral: no input data are rewritten and only the
bounded top-peak fingerprint is retained. The reader already enables OpenMS
`PeakTypeEstimator` when `--estimate-mass-shifts` is active, so unknown spectrum
type can be resolved when the OpenMS runtime has enough evidence.

Diagnostics distinguish native centroid, explicit profile, estimated centroid,
estimated profile, peak-pick failures, and unresolved spectrum type.

## Runtime architecture

The implementation shares the existing single spectrum pass:

```text
normalized MS2 peaks
  -> strongest peaks
  -> bounded coarse-fragment inverted index
  -> small related-spectrum candidate set
  -> unchanged + delta/half-delta fragment matching
  -> charge-aware neutral precursor delta
  -> recurrent one-dimensional delta clustering
  -> isotope/adduct + isotope-adjacent classification
  -> tight OpenMS ModificationsDB / UniMod annotation
  -> bounded QC-facing report + diagnostic accounting
```

The candidate index is bounded in both retained spectra and posting-list size,
avoiding all-vs-all spectral comparisons. Peak arrays are not persisted beyond
the bounded fingerprint representation.

## v19 precision and annotation tolerance

If `--estimate-mass-error` is enabled at the same time, the frozen v19 precursor
precision estimate is consumed **read-only** to widen mass-shift clustering when
supported. v22.1 intentionally decouples this from UniMod annotation: the
UniMod lookup window remains at the configured tight fixed tolerance and does
not inherit a wider run-specific clustering radius.

This prevents noisy runs from turning a broad clustering tolerance into a broad
chemical-database lookup. v19 itself is unchanged.

## QC-facing reporting

The scout retains raw recurrent-cluster accounting but bounds the normal report:

- exact isotope/adduct artifacts are reported as QC artifacts;
- isotope-adjacent satellite clusters are suppressed from the normal PTM list
  and counted diagnostically;
- UniMod decoys and amino-acid substitutions are suppressed from the default
  candidate display;
- candidate-bearing clusters are support-ranked and bounded;
- unmatched/mass-compatible-other clusters remain visible in a bounded
  high-support unknown list;
- diagnostic counters expose raw vs reported cluster counts and suppression
  reasons.

This is reporting hardening, not peptide identification. A cluster can still
have several legitimate mass-compatible candidates.

## Running

```bash
prideqc analyze \
  --estimate-mass-error \
  --estimate-mass-shifts \
  -o results \
  input.raw
```

Using `--estimate-mass-error` is recommended but not required. Without enough
v19 precursor precision support, the scout uses the fixed fallback cluster
window.

Outputs include:

- `<input>.mass-shifts.tsv` for a per-file workflow;
- `mass-shifts.tsv` for aggregate workflow output;
- `putative_modification_mass_shifts` in `annotations.tsv` / summary JSON;
- `mass_shift_scout_diagnostics` with normalization, raw/reportable-cluster and
  suppression accounting;
- the structured candidate table as a local-vocabulary mzQC quality metric.

`mass-shifts.tsv` contains one row per reported recurrent cluster and visible
UniMod candidate. It includes diagnostic candidate counts and the QC candidate
category. A reported cluster with no visible database candidate receives one row
with blank UniMod columns.

## pyOpenMS / OpenMS use

The scout uses OpenMS through pyOpenMS for:

- proton mass and C13-C12 mass-difference constants;
- `PeakTypeEstimator` in the streaming reader for unknown spectrum types;
- `PeakPickerHiRes` for ephemeral profile-MS2 centroiding;
- `ModificationsDB` / UniMod accessions, monoisotopic delta masses,
  specificity/origin metadata, and source classification.

The bounded candidate index and shortlist scorer continue to operate directly on
NumPy arrays. `SpectrumAlignment` remains a benchmark candidate only if profiling
shows shortlist scoring to be a material runtime bottleneck.
