# Identification-free recurrent mass-shift scout (v22)

The v22 scout is an opt-in QC lane for finding recurring precursor neutral-mass
differences between **fragment-related centroid MS2 spectra**. It is intended to
surface modification-like chemistry quickly without making peptide, protein, or
site-localization claims.

## Scientific boundary

A reported mass shift is a hypothesis. A UniMod match means only that the
observed recurrent delta is mass-compatible with one or more modification
records exposed by the pinned OpenMS runtime. The scout does not identify a
sequence, localize a residue, or prove that a PTM is present.

Isotope-selection and common positive-mode adduct substitutions are classified
before UniMod annotation. Unmatched recurrent shifts remain visible rather than
being forced to a database entry. When multiple UniMod records fit the same
mass window, all compatible candidates are preserved up to a bounded cap.

Because evidence comes from pairs of related spectra with different precursor
masses, the initial scout is most sensitive to variable or co-occurring modified
and reference-like forms. A fixed modification carried by every corresponding
peptide can be invisible when no related reference form exists in the run. That
limitation is reported explicitly rather than being filled in from protocol or
search-engine expectations.

## Runtime architecture

The implementation shares the existing single spectrum pass:

```text
MS2 spectrum
  -> strongest centroid peaks
  -> bounded coarse-fragment inverted index
  -> small related-spectrum candidate set
  -> unchanged + delta/half-delta fragment matching
  -> charge-aware neutral precursor delta
  -> recurrent one-dimensional delta clustering
  -> isotope/adduct classification
  -> OpenMS ModificationsDB / UniMod annotation
```

The candidate index is bounded in both the number of retained spectra and the
posting-list size. It avoids all-vs-all spectral comparisons. Peak arrays are
not persisted beyond the bounded top-peak representation.

If `--estimate-mass-error` is enabled at the same time, the frozen v19 precursor
precision estimate is consumed **read-only** to widen the cluster/UniMod mass
window when supported. The v19 estimator, thresholds, and multiplier are not
retuned or modified by v22.

## Running

```bash
prideqc analyze \
  --estimate-mass-error \
  --estimate-mass-shifts \
  -o results \
  input.raw
```

Using `--estimate-mass-error` is recommended but not required. Without enough
v19 precursor precision support, the scout uses conservative fixed fallback
mass windows.

The new outputs are:

- `<input>.mass-shifts.tsv` for a per-file workflow;
- `mass-shifts.tsv` for the aggregate workflow output;
- `putative_modification_mass_shifts` in `annotations.tsv` / summary JSON;
- the same structured candidate table as a local-vocabulary mzQC quality metric.

`mass-shifts.tsv` contains one row per recurrent cluster and mass-compatible
UniMod candidate. A cluster with no database candidate still receives one row
with blank UniMod columns.

## pyOpenMS / OpenMS use

The scout uses OpenMS through pyOpenMS for:

- the proton mass and C13-C12 mass-difference constants;
- the modification catalogue in `ModificationsDB`;
- UniMod accessions, monoisotopic delta masses, specificity/origin metadata, and
  source classification.

The bounded candidate index and shortlist scoring operate directly on the NumPy
arrays already produced by the streaming reader, avoiding repeated construction
of `MSSpectrum` objects at the Python/C++ boundary. `SpectrumAlignment` remains
a benchmark candidate for shortlisted pairs, not a requirement of the initial
v22 path.

Probe the exact SIF before the first cluster experiment:

```bash
singularity exec "$SIF" \
  python /prideqc/scripts/benchmarks/probe_mass_shift_pyopenms.py
```

## Initial scope

The initial lane accepts native centroid MS2 and unknown scans for which the
existing OpenMS peak-type estimator reports centroid. Native profile spectra are
explicitly excluded. A later iteration may benchmark `PeakPickerHiRes` for
profile MS2, but profile support must not change the frozen v19 behavior.
