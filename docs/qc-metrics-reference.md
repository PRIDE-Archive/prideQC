# QC metrics and experimental annotations

This page is generated from the tracked `docs/metrics.tsv` measurement
registry and `docs/metric-presentation.tsv` display policy. The latter
selects which existing mzQC metrics are prominently displayed; it does
**not** change acquisition, calculation, CV identity, serialization or merging.
Experimental method and label annotations are recorded separately.

**Registered QC metric definitions: 30.**
Presentation tiers: **10 core**, **10 diagnostic**,
and **10 extended**. A registration is not a promise
that a value exists for every run. Metrics absent or undefined in the source
mzQC must remain absent/undefined; do not fill missing values with zero.

## Recommended PRIDE website overview

These 10 **core display candidates** cover run structure and acquisition.
Show a value or visualization only when the run has the corresponding
stored mzQC metric. For SRM/MRM or other MS1-free acquisitions, MS1 and
MS1-delimited cycle entries may be inapplicable; this is not QC failure.

The overview uses readable names; the full registry keys, CV accessions,
value shapes and definitions remain available in the sections below.

| Core metric | Applies to | Interpretation / limitation |
| --- | --- | --- |
| Run / MS levels | all | Distinct MS levels recorded; not a quality score. |
| Run time (s) | all | RT span from spectra, not LC separation chemistry. |
| MS1 spectra | ms1 | Includes scans with no detected peaks. |
| MS1 peaks per scan (quartiles) | ms1 | Valid stored peaks per scan, not identified molecular features. |
| MS2 spectra | ms2 | Includes scans with no detected peaks. |
| MS2 peaks per scan (quartiles) | ms2 | Valid stored peaks per scan, not identified molecular features. |
| DIA isolation windows | ms2 | Distinct window geometries; not rounded precursor-target count. DIA evidence and complete metadata required. |
| DIA window widths (min/max) | ms2 | Minimum and maximum distinct window widths in m/z; not individual scalar extrema. |
| MS1/MS2 scan ratio | ms1_ms2 | Undefined when there are no MS2 scans; not a performance score. |
| MS2 scans per acquisition cycle | cycles | MS1-delimited cycles only; not applicable to MS1-free targeted runs. |

Counts include empty scans where stated by the registry. Charge
statistics and DIA windows are reported only when supported. Quantile
and tuple metrics remain multivalued; the overview must not silently collapse
them into an undocumented scalar or recompute them from mzQC.

Acquisition method (DDA/DIA/PRM/SRM), labelling candidates, ionization,
instrument and LC annotations belong in the experimental-context section,
not in a numeric QC score. No universal composite QC score is defined.

## Experimental and technical annotation fields

| Field | Information | Evidence policy |
| --- | --- | --- |
| `acquisition_method` | DDA/DIA/PRM/SRM run classification | Explicit header CV or DDA/DIA inference; targeted scans abstain |
| `acquisition_method_evidence` | PRIDE CV term, support and classifier | Structured mzQC; observed or inferred |
| `acquisition_method_conflict` | Contradictory RAW-header and spectral evidence | No automatic classification/writeback |
| `labeling_strategy` | Verified sample labelling strategy | Unavailable from RAW alone; no SDRF writeback |
| `labeling_strategy_candidates` | TMT/iTRAQ reporter evidence | Optional, exploratory, non-exclusive; structured mzQC |
| `ionization` | Ion source PSI-MS CV terms | Observed when present in instrument configuration |
| `chromatography_type` | LC separation type | Unavailable unless external method metadata is supplied |
| `instrument` | Instrument PSI-MS CV | Observed where supplied |
| `mass_analyzers` | Analyzer PSI-MS CV | Observed where supplied |
| `dissociation_ms2` | MS2 dissociation methods | Observed when scan metadata are complete |
| `collision_energy_ms2` | Recorded collision energies | Observed where present |
| `tmt_reporter_evidence` | TMT centroid MS2/MS3 channel evidence | Optional experimental collector |
| `tmt_profile_reporter_evidence` | TMT profile MS2/MS3 evidence | Optional experimental collector |
| `signature_TMT_family` | Diagnostic TMT reporter-ion screen | Optional, candidate-only |
| `signature_iTRAQ_family` | Diagnostic iTRAQ reporter-ion screen | Optional, candidate-only |
| `estimated_precursor_mass_error_ppm` | Precursor mass precision | Optional estimator |
| `estimated_fragment_mass_error_ppm` | Fragment mass precision | Optional estimator |
| `putative_modification_mass_shifts` | Neutral mass-shift hypotheses | Optional; not peptide IDs |

The PRIDE ontology groups acquisition types under `PRIDE:0000659`
and sample labels under `PRIDE:0000514`.

Acquisition classifications use `PRIDE:0000450` (DIA) and
`PRIDE:0000627` (DDA) only with sufficient spectral evidence.
PRM (`PRIDE:0000629`) and SRM/MRM (`PRIDE:0000630`) are not
automatically called from narrow-isolation scans. Exact PRIDE acquisition
terms recorded in mzML fileContent may establish PRM or SRM/MRM,
unless they conflict with spectral evidence. TMT (`PRIDE:0000314`)
and iTRAQ (`PRIDE:0000313`) may be **candidates**, never sample-label
assignments. SILAC (`PRIDE:0000315`) cannot be inferred from
ordinary precursor charge/isotope information alone. No RAW-only
absence of signals is interpreted as label-free.

## Registered quantitative QC metrics

Tier and applicability are editorial display metadata, not additional
scientific measurements or controlled vocabulary terms.
Browse the compact three-column indexes by category, then expand
a category's definitions to see full CV accessions, value shapes
and the original metric descriptions.

### Acquisition cycles (6)

| Registry metric | Tier | Applies to |
| --- | --- | --- |
| `AcquisitionCycle_Count` | diagnostic | cycles |
| `AcquisitionCycle_MS2Count_Quantiles` | core | cycles |
| `AcquisitionCycle_MS2Count_Mode` | diagnostic | cycles |
| `AcquisitionCycle_TargetCoverageFraction` | diagnostic | cycles |
| `AcquisitionCycle_TargetEligibleFraction` | extended | cycles |
| `AcquisitionCycle_UniqueTargetCount_Quantiles` | extended | cycles |

:::{dropdown} Full definitions and CV accessions — Acquisition cycles

- **`AcquisitionCycle_Count`** — diagnostic; scalar; applies to cycles
  - CV: `QCPRIDE:E5D44EF02DD0C112FCCD`
  - Definition: Count of MS1-delimited cycles containing at least one MS2 scan. MS2 scans before the first MS1 are ignored; a trailing cycle is included when it contains MS2 scans.

- **`AcquisitionCycle_MS2Count_Quantiles`** — core; tuple; applies to cycles
  - CV: `QCPRIDE:CB0E144250E9279F1725`
  - Definition: 25th, 50th, 75th percentiles of MS2 scans per observed MS1-delimited cycle.

- **`AcquisitionCycle_MS2Count_Mode`** — diagnostic; scalar; applies to cycles
  - CV: `QCPRIDE:38064FE0BC760341F5CE`
  - Definition: Most frequent MS2 scan count per observed MS1-delimited cycle; lower count wins ties.

- **`AcquisitionCycle_TargetCoverageFraction`** — diagnostic; scalar; applies to cycles
  - CV: `QCPRIDE:1E7E0C9468F1B503778C`
  - Definition: Fraction of MS2 scans inside observed MS1-delimited cycles with a finite positive first-precursor m/z paired with a positive isolation width.

- **`AcquisitionCycle_TargetEligibleFraction`** — extended; scalar; applies to cycles
  - CV: `QCPRIDE:0EA9DB6C9708D5F09586`
  - Definition: Fraction of observed MS1-delimited cycles eligible for target-set/order comparison.

- **`AcquisitionCycle_UniqueTargetCount_Quantiles`** — extended; tuple; applies to cycles
  - CV: `QCPRIDE:EAE7AB20278D155810C1`
  - Definition: 25th, 50th, 75th percentiles of unique first-precursor isolation targets per eligible cycle after 0.1 Th quantization.

:::

### Chromatography (1)

| Registry metric | Tier | Applies to |
| --- | --- | --- |
| `ChromatographyDuration` | core | all |

:::{dropdown} Full definitions and CV accessions — Chromatography

- **`ChromatographyDuration`** — core; scalar; applies to all
  - CV: `MS:4000053`
  - Definition: Last minus first finite spectrum RT, over all MS levels.

:::

### Ion intensity and TIC (1)

| Registry metric | Tier | Applies to |
| --- | --- | --- |
| `RT_TIC_Quantiles` | extended | ms1 |

:::{dropdown} Full definitions and CV accessions — Ion intensity and TIC

- **`RT_TIC_Quantiles`** — extended; tuple; applies to ms1
  - CV: `MS:4000183`
  - Definition: Four intervals at first scan-wise cumulative MS1 TIC quartile crossings, divided by the MS1 time span; this uses accumulation, not integration.

:::

### Ion mobility / FAIMS (3)

| Registry metric | Tier | Applies to |
| --- | --- | --- |
| `FAIMS_CV_Values` | diagnostic | faims |
| `FAIMS_CV_Count` | diagnostic | faims |
| `FAIMS_CV_Range` | diagnostic | faims |

:::{dropdown} Full definitions and CV accessions — Ion mobility / FAIMS

- **`FAIMS_CV_Values`** — diagnostic; tuple; applies to faims
  - CV: `QCPRIDE:E69917C56F3DCCD6D0EE`
  - Definition: Distinct finite FAIMS compensation voltages, sorted, in volts; not drift times.

- **`FAIMS_CV_Count`** — diagnostic; scalar; applies to faims
  - CV: `QCPRIDE:0FF9AAA5ADA010E8801C`
  - Definition: Number of distinct finite FAIMS compensation voltages.

- **`FAIMS_CV_Range`** — diagnostic; tuple; applies to faims
  - CV: `QCPRIDE:3C721D04BCC72F40C1A0`
  - Definition: Minimum and maximum finite FAIMS compensation voltages, in volts.

:::

### Other run-level QC (4)

| Registry metric | Tier | Applies to |
| --- | --- | --- |
| `ObservedMzRange_MS1` | diagnostic | ms1 |
| `AvgCycleTime_MS1` | diagnostic | ms1 |
| `ObservedMzRange_MS2` | diagnostic | ms2 |
| `MS1_to_MS2_Ratio` | core | ms1_ms2 |

:::{dropdown} Full definitions and CV accessions — Other run-level QC

- **`ObservedMzRange_MS1`** — diagnostic; tuple; applies to ms1
  - CV: `QCPRIDE:2A7C646C00FBCB34F4F8`
  - Definition: Min/max valid peak m/z; this is an observed range, not the instrument scan window, MS1.

- **`AvgCycleTime_MS1`** — diagnostic; scalar; applies to ms1
  - CV: `QCPRIDE:2C6256D0D438E08FF913`
  - Definition: Mean strictly positive adjacent sorted MS1 RT differences, in seconds.

- **`ObservedMzRange_MS2`** — diagnostic; tuple; applies to ms2
  - CV: `QCPRIDE:47E4D2D8B1F1F4D61C19`
  - Definition: Min/max valid peak m/z; this is an observed range, not the instrument scan window, MS2.

- **`MS1_to_MS2_Ratio`** — core; scalar; applies to ms1_ms2
  - CV: `QCPRIDE:3E6B36C0810EEDA24D7A`
  - Definition: MS1 scan count divided by MS2 scan count; undefined when denominator is zero.

:::

### Precursors and isolation (6)

| Registry metric | Tier | Applies to |
| --- | --- | --- |
| `DIAIsolationWindowCount` | core | ms2 |
| `DIAIsolationWindowWidths` | core | ms2 |
| `ChargeMean` | extended | ms2 |
| `ChargeMedian` | diagnostic | ms2 |
| `ChargeRatio_3over2` | extended | ms2 |
| `ChargeRatio_4over2` | extended | ms2 |

:::{dropdown} Full definitions and CV accessions — Precursors and isolation

- **`DIAIsolationWindowCount`** — core; scalar; applies to ms2
  - CV: `MS:4000194`
  - Definition: Number of distinct DIA MS2 isolation windows by low/high m/z boundaries and available FAIMS and ion mobility window settings; emitted only for complete, DIA-supported observations.

- **`DIAIsolationWindowWidths`** — core; tuple; applies to ms2
  - CV: `MS:4000195`
  - Definition: Minimum and maximum width in m/z across distinct DIA MS2 isolation windows; fixed-width DIA repeats the same value twice.

- **`ChargeMean`** — extended; scalar; applies to ms2
  - CV: `MS:4000173`
  - Definition: Mean positive first-precursor MS2 charge; unknown charges excluded.

- **`ChargeMedian`** — diagnostic; scalar; applies to ms2
  - CV: `MS:4000175`
  - Definition: Median positive first-precursor MS2 charge; unknown charges excluded.

- **`ChargeRatio_3over2`** — extended; scalar; applies to ms2
  - CV: `MS:4000169`
  - Definition: Count of MS2 charge 3 divided by charge 2; unavailable unless both are present (rawQC convention).

- **`ChargeRatio_4over2`** — extended; scalar; applies to ms2
  - CV: `MS:4000171`
  - Definition: Count of MS2 charge 4 divided by charge 2; unavailable unless both are present (rawQC convention).

:::

### Scan timing (4)

| Registry metric | Tier | Applies to |
| --- | --- | --- |
| `FastestFrequency_MS1` | extended | ms1 |
| `RT_MS1_Quantiles` | extended | ms1 |
| `FastestFrequency_MS2` | extended | ms2 |
| `RT_MS2_Quantiles` | extended | ms2 |

:::{dropdown} Full definitions and CV accessions — Scan timing

- **`FastestFrequency_MS1`** — extended; scalar; applies to ms1
  - CV: `MS:4000065`
  - Definition: Maximum scan count in a closed 60-second window divided by 60, including short runs.

- **`RT_MS1_Quantiles`** — extended; tuple; applies to ms1
  - CV: `MS:4000184`
  - Definition: Four consecutive scan-RT quartile widths divided by the MS1 time span.

- **`FastestFrequency_MS2`** — extended; scalar; applies to ms2
  - CV: `MS:4000066`
  - Definition: Maximum scan count in a closed 60-second window divided by 60, including short runs.

- **`RT_MS2_Quantiles`** — extended; tuple; applies to ms2
  - CV: `MS:4000185`
  - Definition: Four consecutive scan-RT quartile widths divided by the MS2 time span.

:::

### Spectra and peak summaries (5)

| Registry metric | Tier | Applies to |
| --- | --- | --- |
| `NumberOfMSLevels` | core | all |
| `NumberOfSpectra_MS1` | core | ms1 |
| `PeakDensity_MS1_Quantiles` | core | ms1 |
| `NumberOfSpectra_MS2` | core | ms2 |
| `PeakDensity_MS2_Quantiles` | core | ms2 |

:::{dropdown} Full definitions and CV accessions — Spectra and peak summaries

- **`NumberOfMSLevels`** — core; scalar; applies to all
  - CV: `QCPRIDE:C782CFB84B18715F8D67`
  - Definition: Number of distinct observed MS levels.

- **`NumberOfSpectra_MS1`** — core; scalar; applies to ms1
  - CV: `MS:4000059`
  - Definition: Count of all MS1 spectra, including empty scans.

- **`PeakDensity_MS1_Quantiles`** — core; tuple; applies to ms1
  - CV: `MS:4000061`
  - Definition: 25th, 50th, 75th percentiles of valid peaks per MS1 spectrum; empty scans included.

- **`NumberOfSpectra_MS2`** — core; scalar; applies to ms2
  - CV: `MS:4000060`
  - Definition: Count of all MS2 spectra, including empty scans.

- **`PeakDensity_MS2_Quantiles`** — core; tuple; applies to ms2
  - CV: `MS:4000062`
  - Definition: 25th, 50th, 75th percentiles of valid peaks per MS2 spectrum; empty scans included.

:::

## Provenance and caveats

The selected set of calculated QC metrics remains in mzQC per run
and in losslessly consolidated accession mzQC. A website or
MongoDB projection may select **stored values only**, never generate
report-only scientific measurements. Grouped spectrum counts, MS2/MS3
representation, exploratory reporter-ion evidence and denominators
must retain their original definitions. A missing measurement is
not a zero value, low score, or proof that an experimental method
was absent.

`NumberOfSpectra_MS1`, `NumberOfSpectra_MS2`,
`DIAIsolationWindowWidths`, `ChromatographyDuration`,
and the supported DIA window statistics remain in the scientific registry.
Ion-source CV metadata are observed only if the reader supplies it.
RT/duration do not reveal LC separation chemistry. Label/channel
annotations require verified submitter/SDRF evidence for writeback.

See also [QC workflow](qc-workflow.md) and [measurement evidence](measurement-evidence.md).
