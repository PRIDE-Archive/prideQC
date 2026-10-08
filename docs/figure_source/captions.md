# prideQC figures: captions and placement

All figures are editable SVG (real `<text>`, no embedded raster), 1400 px wide viewBox, white background,
one palette and one Arial/Helvetica font stack. All data shown are **synthetic**; where a number is printed it was
produced by prideQC's own code (`metrics.py`, `mass_error.py`, `mass_shift.py`, `cohort.py`) and the OpenMS/UniMod
catalogue on the synthetic input. Regenerate intentionally with `python docs/figure_source/build_all.py` from the prideQC checkout; the source imports the local prideQC implementation. Read the Docs uses committed SVGs and never runs the generators.

---

## Figure 1: prideQC overview
**File:** `docs/_static/figures/figure1_overview.svg` | **Place:** `docs/index.md` (top of the landing page, after the one-paragraph summary)

**Caption:** prideQC reads each mzML or supported vendor RAW file once and dispatches every spectrum to evidence collectors
(scan and chromatography, acquisition, mass-error precision, diagnostic ions, recurrent mass shifts). Results are written
as mzQC, TSV tables and summary JSON. mzQC files can be aggregated by pmultiqc / MultiQC into an interactive cross-run
report; per-run evidence can be synthesised across a cohort and, where it passes confidence and scope gates, used to
refine and validate SDRF metadata. Uncertain evidence abstains.

**Alt text:** Left-to-right workflow from raw files through a streaming QC engine to structured outputs, which feed a reporting branch and a gated SDRF refinement branch.

---

## Figure 2: How QC metrics are computed
**File:** `docs/_static/figures/figure2_qc_metrics.svg` | **Place:** `docs/qc-workflow.md`, section on QC metrics

**Caption:** A synthetic 90-minute DDA run is reduced to compact metrics. (A) TIC trace and a 20 s zoom of MS1 and MS2 scans,
MS2 coloured by precursor charge. (B) TIC area by trapezoidal integration and TIC coefficient of variation. (C) Base-peak
intensity per MS1 scan. (D) Scan rate, fastest 60 s window and MS1 cycle time. (E) Precursor charge distribution.
(F) Isolation-window widths against the 15 Th boundary. (G) Repeated-target behaviour, contrasted with a fixed target
list. Values are computed by prideQC's metric code on the synthetic scans.

**Alt text:** Raw scan stream with six small quantitative panels and formula callouts.

---

## Figure 3: Identification-free mass-error precision
**File:** `docs/_static/figures/figure3_mass_error.svg` | **Place:** `docs/measurement-evidence.md`, section on mass-error precision

**Caption:** Repeated observations of the same precursor (same charge, within 120 s, within 20 ppm of the cluster mean) and
likely repeated MS2 spectra give pairwise mass errors (A, B). A robust median/MAD scale, converted to a Gaussian-equivalent
sigma (1.4826 x MAD) and divided by the square root of two, gives single-measurement precision (C). Six times that precision
is the supported tolerance suggestion when support gates pass (D). For fragments, the resolution regime decides the unit
(ppm or Da); intermediate or censored cases abstain (E). High-resolution fragment precision in the OpenMS estimator comes
from a Gaussian plus uniform-background fit. The result measures precision, not the historical search parameters.

**Alt text:** Four connected quantitative panels (repeat observations, error histogram, MAD cumulative curve, tolerance region) and a regime decision strip.

---

## Figure 4: Recurrent mass-shift scouting
**File:** `docs/_static/figures/figure4_mass_shift.svg` | **Place:** `docs/measurement-evidence.md`, section on mass-shift evidence (or a dedicated page if it is split out)

**Caption:** Two related synthetic MS2 spectra are reduced to their strongest peaks (1), indexed by shared 1 Da bins (2) and
compared for unchanged fragment matches and matches shifted by the full precursor difference or by half of it (3). The
accepted pair yields a neutral mass difference (4). Many accepted pairs form recurrent families (5); the recurrence window
is calibrated from mass-error precision, whereas the chemical annotation window is a tighter, fixed +/-0.020 Da. Families are
checked for isotope and adduct explanations and looked up in the OpenMS/UniMod catalogue (6). Mass compatibility does not
establish modification identity or site localization.

**Alt text:** Mirrored spectra, shared-bin barcode, three match tests, accepted-pair gates, a delta-mass histogram with two windows, and a catalogue number line.

---

## Figure 5: Cohort synthesis and SDRF refinement
**File:** `docs/_static/figures/figure5_cohort_sdrf.svg` | **Place:** `docs/sdrf-refinement.md`, top of the page

**Caption:** Runs are first partitioned into inferred experiment groups; gates apply within a group (A). The tolerance
branch requires a supported estimate in at least 80% of runs, one compatible unit, takes the maximum supported per-run value
(a shared tolerance must cover the noisiest run) and rounds upward (B). The modification branch passes a mass-shift family
through prevalence, unique-identity, category and independent-evidence gates; families that fail stop, and families that
only lack independent evidence are held as putative entries (C). Accepted proposals are applied to the SDRF, validated and
logged (D).

**Alt text:** Run grouping scatter, a five-step tolerance rail with a maximum-versus-mean chart, a gate matrix for five example mass-shift families, and an original-versus-refined SDRF table with audit files.

---

## Figure 6: Cross-run reporting with pmultiqc / MultiQC
**File:** `docs/_static/figures/figure6_reporting.svg` | **Place:** `docs/reporting.md`, top of the page

**Caption:** One mzQC file per input run is aggregated by pmultiqc / MultiQC into a single interactive HTML report with
run-to-run distributions, outlier runs, acquisition characteristics and comparison views. The reporting layer displays values
that prideQC has already computed; it does not repeat the scientific calculations, and SDRF refinement reads each run's
summary JSON rather than this report.

**Alt text:** Seven mzQC files feeding a collector, then a browser-style report with four example plots.
