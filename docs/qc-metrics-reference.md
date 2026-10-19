# QC metrics and experimental annotations

This page is generated from the tracked `docs/metrics.tsv` measurement
registry and `docs/metric-presentation.tsv` display policy. The latter
selects which existing mzQC metrics are prominently displayed; it does
**not** change acquisition, calculation, CV identity, serialization or merging.
Experimental method and label annotations are recorded separately.

**Registered QC metric definitions: 104.**
Presentation tiers: **18 core**, **39 diagnostic**,
and **47 extended**. A registration is not a promise
that a value exists for every run. Metrics absent or undefined in the source
mzQC must remain absent/undefined; do not fill missing values with zero.

## Recommended PRIDE website overview

These 18 **core display candidates** cover run structure, signal and acquisition.
Show a value or visualization only when the run has the corresponding
stored mzQC metric. For SRM/MRM or other MS1-free acquisitions, MS1 and
MS1-delimited cycle entries may be inapplicable; this is not QC failure.

The overview uses readable names; the full registry keys, CV accessions,
value shapes and definitions remain available in the sections below.

| Core metric | Applies to | Interpretation / limitation |
| --- | --- | --- |
| Run / MS levels | all | Distinct MS levels recorded; not a quality score. |
| Total spectra | all | Counts every MS level including empty scans. |
| Run time (s) | all | RT span from spectra, not LC separation chemistry. |
| MS1 spectra | ms1 | Includes scans with no detected peaks. |
| Empty MS1 scans | ms1 | Zero-TIC scans; interpret alongside the MS1 count. |
| MS1 peaks per scan (quartiles) | ms1 | Valid stored peaks per scan, not identified molecular features. |
| Median MS1 TIC | ms1 | Arbitrary instrument units; do not compare without calibration. |
| MS1 TIC variability (CV) | ms1 | Undefined if the mean TIC is zero or too few scans are available. |
| MS1 scan rate | ms1 | Scans per minute over the MS1 RT span. |
| MS2 spectra | ms2 | Includes scans with no detected peaks. |
| Empty MS2 scans | ms2 | Zero-TIC scans; interpret alongside the MS2 count. |
| MS2 peaks per scan (quartiles) | ms2 | Valid stored peaks per scan, not identified molecular features. |
| Median MS2 TIC | ms2 | Arbitrary instrument units; do not compare without calibration. |
| MS2 scan rate | ms2 | Scans per minute over the MS2 RT span. |
| Median MS2 isolation width | ms2 | Requires positive recorded isolation offsets; absent is not zero. |
| MS1/MS2 scan ratio | ms1_ms2 | Undefined when there are no MS2 scans; not a performance score. |
| MS2 scans per acquisition cycle | cycles | MS1-delimited cycles only; not applicable to MS1-free targeted runs. |
| MS2 precursor charge fractions | ms2 | Denominator is all MS2 scans; unknown charge is retained. |

The counts include empty scans where stated by the registry. MS1/MS2
TIC values have arbitrary instrument-dependent intensity units and are
**not comparable across uncalibrated runs or instrument types**. Charge
fractions retain unknowns in the all-MS2 denominator. Quantile and table
metrics remain multivalued; the overview must not silently collapse
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

### Acquisition cycles (12)

| Registry metric | Tier | Applies to |
| --- | --- | --- |
| `AcquisitionCycle_Count` | diagnostic | cycles |
| `AcquisitionCycle_MS2Count_Quantiles` | core | cycles |
| `AcquisitionCycle_MS2Count_Mode` | diagnostic | cycles |
| `AcquisitionCycle_MS2Count_ModalFraction` | extended | cycles |
| `AcquisitionCycle_TargetCoverageFraction` | diagnostic | cycles |
| `AcquisitionCycle_TargetEligibleCount` | extended | cycles |
| `AcquisitionCycle_TargetEligibleFraction` | extended | cycles |
| `AcquisitionCycle_UniqueTargetCount_Quantiles` | extended | cycles |
| `AcquisitionCycle_TargetSetDistinctCount` | extended | cycles |
| `AcquisitionCycle_TargetSetModalFraction` | diagnostic | cycles |
| `AcquisitionCycle_TargetOrderDistinctCount` | extended | cycles |
| `AcquisitionCycle_TargetOrderModalFraction` | extended | cycles |

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

- **`AcquisitionCycle_MS2Count_ModalFraction`** — extended; scalar; applies to cycles
  - CV: `QCPRIDE:F6EA9ED8D480A916D62F`
  - Definition: Fraction of observed MS1-delimited cycles having the modal MS2 scan count.

- **`AcquisitionCycle_TargetCoverageFraction`** — diagnostic; scalar; applies to cycles
  - CV: `QCPRIDE:1E7E0C9468F1B503778C`
  - Definition: Fraction of MS2 scans inside observed MS1-delimited cycles with a finite positive first-precursor m/z paired with a positive isolation width.

- **`AcquisitionCycle_TargetEligibleCount`** — extended; scalar; applies to cycles
  - CV: `QCPRIDE:172D73B0C593B7E3F376`
  - Definition: Count of observed MS1-delimited cycles where every MS2 scan has target metadata eligible for target-set/order comparison.

- **`AcquisitionCycle_TargetEligibleFraction`** — extended; scalar; applies to cycles
  - CV: `QCPRIDE:0EA9DB6C9708D5F09586`
  - Definition: Fraction of observed MS1-delimited cycles eligible for target-set/order comparison.

- **`AcquisitionCycle_UniqueTargetCount_Quantiles`** — extended; tuple; applies to cycles
  - CV: `QCPRIDE:EAE7AB20278D155810C1`
  - Definition: 25th, 50th, 75th percentiles of unique first-precursor isolation targets per eligible cycle after 0.1 Th quantization.

- **`AcquisitionCycle_TargetSetDistinctCount`** — extended; scalar; applies to cycles
  - CV: `QCPRIDE:3DAD05F0A9B6418D738F`
  - Definition: Number of distinct per-cycle isolation-target sets among eligible cycles after 0.1 Th quantization.

- **`AcquisitionCycle_TargetSetModalFraction`** — diagnostic; scalar; applies to cycles
  - CV: `QCPRIDE:C1EE0988D7C3648932F1`
  - Definition: Fraction of eligible cycles whose quantized isolation-target set equals the modal set.

- **`AcquisitionCycle_TargetOrderDistinctCount`** — extended; scalar; applies to cycles
  - CV: `QCPRIDE:971C73C2545A365A4D7E`
  - Definition: Number of distinct ordered per-cycle isolation-target sequences among eligible cycles after 0.1 Th quantization.

- **`AcquisitionCycle_TargetOrderModalFraction`** — extended; scalar; applies to cycles
  - CV: `QCPRIDE:4A5FB670C1ACFB1FBDB2`
  - Definition: Fraction of eligible cycles whose quantized ordered isolation-target sequence equals the modal sequence.

:::

### Chromatography (5)

| Registry metric | Tier | Applies to |
| --- | --- | --- |
| `ChromatographyDuration` | core | all |
| `NumberOfChromatograms` | diagnostic | chromatograms |
| `NumberOfChromatogramDataPoints` | diagnostic | chromatograms |
| `ChromatogramTypes` | diagnostic | chromatograms |
| `ChromatogramRTRange` | extended | chromatograms |

:::{dropdown} Full definitions and CV accessions — Chromatography

- **`ChromatographyDuration`** — core; scalar; applies to all
  - CV: `MS:4000053`
  - Definition: Last minus first finite spectrum RT, over all MS levels.

- **`NumberOfChromatograms`** — diagnostic; scalar; applies to chromatograms
  - CV: `QCPRIDE:3A03550E7F73CC8BC2E7`
  - Definition: Number of stored chromatograms.

- **`NumberOfChromatogramDataPoints`** — diagnostic; scalar; applies to chromatograms
  - CV: `QCPRIDE:166BE0366D9F74CDD72A`
  - Definition: Total stored chromatogram array lengths, not chromatographic peaks.

- **`ChromatogramTypes`** — diagnostic; table; applies to chromatograms
  - CV: `QCPRIDE:F6028B02FCB099FDD16A`
  - Definition: Counts grouped by OpenMS ChromatogramType enum value; engine version is recorded.

- **`ChromatogramRTRange`** — extended; tuple; applies to chromatograms
  - CV: `QCPRIDE:DA479CE356AAFB988099`
  - Definition: Minimum and maximum finite chromatogram time, in seconds.

:::

### Ion intensity and TIC (22)

| Registry metric | Tier | Applies to |
| --- | --- | --- |
| `BasePeak_All_Max` | extended | all |
| `BasePeak_MS1_Mean` | extended | ms1 |
| `TIC_MS1_Median` | core | ms1 |
| `TIC_MS1_CV` | core | ms1 |
| `TIC_MS1_Area` | diagnostic | ms1 |
| `TIC_MS1_Area_RTQuantiles` | extended | ms1 |
| `TIC_MS1_QuartileLogRatios` | extended | ms1 |
| `TIC_MS1_ChangeQuartileLogRatios` | extended | ms1 |
| `TIC_MS1_SignalJump10x_Count` | diagnostic | ms1 |
| `TIC_MS1_SignalFall10x_Count` | diagnostic | ms1 |
| `MedianTIC_in_RT_MS1_IQR` | extended | ms1 |
| `TIC_MS1_MedianInHalfRange` | diagnostic | ms1 |
| `RT_TIC_Quantiles` | extended | ms1 |
| `BasePeak_MS2_Mean` | extended | ms2 |
| `TIC_MS2_Median` | core | ms2 |
| `TIC_MS2_CV` | diagnostic | ms2 |
| `TIC_MS2_Area` | diagnostic | ms2 |
| `TIC_MS2_Area_RTQuantiles` | extended | ms2 |
| `TIC_MS2_QuartileLogRatios` | extended | ms2 |
| `TIC_MS2_ChangeQuartileLogRatios` | extended | ms2 |
| `TIC_MS2_SignalJump10x_Count` | diagnostic | ms2 |
| `TIC_MS2_SignalFall10x_Count` | diagnostic | ms2 |

:::{dropdown} Full definitions and CV accessions — Ion intensity and TIC

- **`BasePeak_All_Max`** — extended; scalar; applies to all
  - CV: `QCPRIDE:0DFA7AA5174921CAC620`
  - Definition: Largest valid intensity across all spectra, including MS3 and higher; arbitrary intensity units.

- **`BasePeak_MS1_Mean`** — extended; scalar; applies to ms1
  - CV: `QCPRIDE:D452DD102903397E2819`
  - Definition: Mean of per-spectrum maximum valid intensity, excluding scans with no valid peaks, MS1.

- **`TIC_MS1_Median`** — core; scalar; applies to ms1
  - CV: `QCPRIDE:C9C317E0A64A03006F7D`
  - Definition: Median sum of valid peak intensities per scan, MS1.

- **`TIC_MS1_CV`** — core; scalar; applies to ms1
  - CV: `QCPRIDE:00C4C4D556C5B163B0FE`
  - Definition: Sample standard deviation (ddof=1) divided by mean TIC, MS1.

- **`TIC_MS1_Area`** — diagnostic; scalar; applies to ms1
  - CV: `QCPRIDE:01859CE82042A7057DDF`
  - Definition: Trapezoidal TIC integral against RT, stable time order, intensity times seconds, MS1.

- **`TIC_MS1_Area_RTQuantiles`** — extended; tuple; applies to ms1
  - CV: `QCPRIDE:59BA94235C553AF4FDC6`
  - Definition: Four exact integrals of piecewise-linear TIC split at scan-RT quartiles, intensity times seconds, MS1.

- **`TIC_MS1_QuartileLogRatios`** — extended; tuple; applies to ms1
  - CV: `QCPRIDE:604EE74E7854759E05FC`
  - Definition: Natural logs of TIC Q2/Q1, Q3/Q2, Q4/Q3 (Q4=max), MS1.

- **`TIC_MS1_ChangeQuartileLogRatios`** — extended; tuple; applies to ms1
  - CV: `QCPRIDE:B4CB91D742E1BE7A0B5E`
  - Definition: Natural logs of absolute adjacent TIC change Q2/Q1, Q3/Q2, Q4/Q3, MS1.

- **`TIC_MS1_SignalJump10x_Count`** — diagnostic; scalar; applies to ms1
  - CV: `QCPRIDE:CC18B4EE12E3DE4D10CC`
  - Definition: Adjacent sorted-TIC ratios >=10, with strictly positive previous TIC, MS1.

- **`TIC_MS1_SignalFall10x_Count`** — diagnostic; scalar; applies to ms1
  - CV: `QCPRIDE:968FA11A94C7458E4F1C`
  - Definition: Adjacent sorted-TIC ratios <=0.1, with strictly positive previous TIC; zero next TIC counts, MS1.

- **`MedianTIC_in_RT_MS1_IQR`** — extended; scalar; applies to ms1
  - CV: `QCPRIDE:384CE930503A31FD6181`
  - Definition: Median TIC in middle two scan-index quartile groups (recycle-and-sort assignment), MS1. ID-free proxy; not an identified-peptide metric.

- **`TIC_MS1_MedianInHalfRange`** — diagnostic; tuple; applies to ms1
  - CV: `QCPRIDE:866968003763CDDE9CF5`
  - Definition: Median TIC in the shortest window containing ceil(n/2) finite-RT MS1 scans; first window wins ties. ID-free proxy.

- **`RT_TIC_Quantiles`** — extended; tuple; applies to ms1
  - CV: `MS:4000183`
  - Definition: Four intervals at first scan-wise cumulative MS1 TIC quartile crossings, divided by the MS1 time span; this uses accumulation, not integration.

- **`BasePeak_MS2_Mean`** — extended; scalar; applies to ms2
  - CV: `QCPRIDE:25C2A0E5AF8883957F5E`
  - Definition: Mean of per-spectrum maximum valid intensity, excluding scans with no valid peaks, MS2.

- **`TIC_MS2_Median`** — core; scalar; applies to ms2
  - CV: `QCPRIDE:D249949F6A3F6A3EBF3E`
  - Definition: Median sum of valid peak intensities per scan, MS2.

- **`TIC_MS2_CV`** — diagnostic; scalar; applies to ms2
  - CV: `QCPRIDE:8AF7F9131561A60EDB62`
  - Definition: Sample standard deviation (ddof=1) divided by mean TIC, MS2.

- **`TIC_MS2_Area`** — diagnostic; scalar; applies to ms2
  - CV: `QCPRIDE:9FD64B0DA7A59E8514E2`
  - Definition: Trapezoidal TIC integral against RT, stable time order, intensity times seconds, MS2.

- **`TIC_MS2_Area_RTQuantiles`** — extended; tuple; applies to ms2
  - CV: `QCPRIDE:E191EB22593B10C83302`
  - Definition: Four exact integrals of piecewise-linear TIC split at scan-RT quartiles, intensity times seconds, MS2.

- **`TIC_MS2_QuartileLogRatios`** — extended; tuple; applies to ms2
  - CV: `QCPRIDE:80B75AE92157AF8E1090`
  - Definition: Natural logs of TIC Q2/Q1, Q3/Q2, Q4/Q3 (Q4=max), MS2.

- **`TIC_MS2_ChangeQuartileLogRatios`** — extended; tuple; applies to ms2
  - CV: `QCPRIDE:760CEDD06B184EE701FF`
  - Definition: Natural logs of absolute adjacent TIC change Q2/Q1, Q3/Q2, Q4/Q3, MS2.

- **`TIC_MS2_SignalJump10x_Count`** — diagnostic; scalar; applies to ms2
  - CV: `QCPRIDE:4BE7DF39FCE21F9E3583`
  - Definition: Adjacent sorted-TIC ratios >=10, with strictly positive previous TIC, MS2.

- **`TIC_MS2_SignalFall10x_Count`** — diagnostic; scalar; applies to ms2
  - CV: `QCPRIDE:DCF0DBC8D6F83087ED6E`
  - Definition: Adjacent sorted-TIC ratios <=0.1, with strictly positive previous TIC; zero next TIC counts, MS2.

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

### Other run-level QC (9)

| Registry metric | Tier | Applies to |
| --- | --- | --- |
| `EmptyScans_MS1` | core | ms1 |
| `ObservedMzRange_MS1` | diagnostic | ms1 |
| `ScanWindow_MS1` | diagnostic | ms1 |
| `AvgCycleTime_MS1` | diagnostic | ms1 |
| `EmptyScans_MS2` | core | ms2 |
| `ObservedMzRange_MS2` | diagnostic | ms2 |
| `ScanWindow_MS2` | diagnostic | ms2 |
| `MzRange_MS2` | extended | ms2 |
| `MS1_to_MS2_Ratio` | core | ms1_ms2 |

:::{dropdown} Full definitions and CV accessions — Other run-level QC

- **`EmptyScans_MS1`** — core; scalar; applies to ms1
  - CV: `QCPRIDE:6F59F45176612D8C6992`
  - Definition: MS1 scans with no valid positive intensity (zero TIC).

- **`ObservedMzRange_MS1`** — diagnostic; tuple; applies to ms1
  - CV: `QCPRIDE:2A7C646C00FBCB34F4F8`
  - Definition: Min/max valid peak m/z; this is an observed range, not the instrument scan window, MS1.

- **`ScanWindow_MS1`** — diagnostic; tuple; applies to ms1
  - CV: `QCPRIDE:05EB1DDF8CCEE51BE48C`
  - Definition: Envelope of recorded scan-window bounds in Th (m/z), MS1.

- **`AvgCycleTime_MS1`** — diagnostic; scalar; applies to ms1
  - CV: `QCPRIDE:2C6256D0D438E08FF913`
  - Definition: Mean strictly positive adjacent sorted MS1 RT differences, in seconds.

- **`EmptyScans_MS2`** — core; scalar; applies to ms2
  - CV: `QCPRIDE:3BE75D0FF4AD59C42832`
  - Definition: MS2 scans with no valid positive intensity (zero TIC).

- **`ObservedMzRange_MS2`** — diagnostic; tuple; applies to ms2
  - CV: `QCPRIDE:47E4D2D8B1F1F4D61C19`
  - Definition: Min/max valid peak m/z; this is an observed range, not the instrument scan window, MS2.

- **`ScanWindow_MS2`** — diagnostic; tuple; applies to ms2
  - CV: `QCPRIDE:BBD580CB1663521C4976`
  - Definition: Envelope of recorded scan-window bounds in Th (m/z), MS2.

- **`MzRange_MS2`** — extended; tuple; applies to ms2
  - CV: `QCPRIDE:3432120DEC9FAE1739E7`
  - Definition: Range of positive finite first-precursor m/z, not fragment peak m/z, MS2.

- **`MS1_to_MS2_Ratio`** — core; scalar; applies to ms1_ms2
  - CV: `QCPRIDE:3E6B36C0810EEDA24D7A`
  - Definition: MS1 scan count divided by MS2 scan count; undefined when denominator is zero.

:::

### Precursors and isolation (26)

| Registry metric | Tier | Applies to |
| --- | --- | --- |
| `PrecursorMz_MS2_Median` | diagnostic | ms2 |
| `PrecursorIntensity_MS2_Quantiles` | extended | ms2 |
| `PrecursorIntensity_MS2_Mean` | extended | ms2 |
| `PrecursorIntensity_MS2_Sd` | extended | ms2 |
| `PrecursorIntensity_MS2_MissingCount` | diagnostic | ms2 |
| `IsolationWidth_MS2_Median` | core | ms2 |
| `IsolationWidth_MS2_Quantiles` | diagnostic | ms2 |
| `IsolationWidth_MS2_Min` | extended | ms2 |
| `IsolationWidth_MS2_Max` | extended | ms2 |
| `IsolationWidth_MS2_Count` | diagnostic | ms2 |
| `IsolationWidth_MS2_FractionLe15` | extended | ms2 |
| `IsolationWidth_MS2_FractionGe15` | extended | ms2 |
| `IsolationPrecursorMz_MS2_Count` | extended | ms2 |
| `IsolationPrecursorMz_MS2_UniqueCount` | extended | ms2 |
| `IsolationPrecursorMz_MS2_UniqueFraction` | diagnostic | ms2 |
| `IsolationPrecursorMz_MS2_MaxRepeatFraction` | extended | ms2 |
| `IsolationPrecursorMz_MS2_RepeatedFraction` | extended | ms2 |
| `MultiplePrecursors_MS2_Count` | diagnostic | ms2 |
| `ExtentPrecursorIntensity_95over5_MS2` | extended | ms2 |
| `MS2_PrecursorCharge_Fractions` | core | ms2 |
| `ChargeMean` | extended | ms2 |
| `ChargeMedian` | diagnostic | ms2 |
| `ChargeMin` | extended | ms2 |
| `ChargeMax` | extended | ms2 |
| `ChargeRatio_3over2` | extended | ms2 |
| `ChargeRatio_4over2` | extended | ms2 |

:::{dropdown} Full definitions and CV accessions — Precursors and isolation

- **`PrecursorMz_MS2_Median`** — diagnostic; scalar; applies to ms2
  - CV: `QCPRIDE:0DB185A646DB806298FC`
  - Definition: Median positive finite first-precursor m/z, MS2.

- **`PrecursorIntensity_MS2_Quantiles`** — extended; tuple; applies to ms2
  - CV: `QCPRIDE:8FB2EFF7D173D0D4C3A2`
  - Definition: Quartiles of recorded positive first-precursor intensity, MS2; no fragment-TIC substitution.

- **`PrecursorIntensity_MS2_Mean`** — extended; scalar; applies to ms2
  - CV: `QCPRIDE:36D2085F9869CE667625`
  - Definition: Mean recorded positive first-precursor intensity, MS2.

- **`PrecursorIntensity_MS2_Sd`** — extended; scalar; applies to ms2
  - CV: `QCPRIDE:AA281BF4225DD99B5942`
  - Definition: Sample SD of recorded positive first-precursor intensity (ddof=1), MS2.

- **`PrecursorIntensity_MS2_MissingCount`** — diagnostic; scalar; applies to ms2
  - CV: `QCPRIDE:37ABD88F90B17EAA841D`
  - Definition: Scans lacking a positive finite first-precursor intensity, including the OpenMS unset default zero, MS2.

- **`IsolationWidth_MS2_Median`** — core; scalar; applies to ms2
  - CV: `QCPRIDE:652589C44AB88BC96F40`
  - Definition: Median positive first-precursor lower+upper isolation offsets in Th (m/z), MS2.

- **`IsolationWidth_MS2_Quantiles`** — diagnostic; tuple; applies to ms2
  - CV: `QCPRIDE:482962258BA982FED6BB`
  - Definition: 25th, 50th, 75th percentiles of positive first-precursor isolation widths in Th (m/z), MS2.

- **`IsolationWidth_MS2_Min`** — extended; scalar; applies to ms2
  - CV: `QCPRIDE:228CD009EF66CB6342FD`
  - Definition: Minimum positive first-precursor isolation width in Th (m/z), MS2.

- **`IsolationWidth_MS2_Max`** — extended; scalar; applies to ms2
  - CV: `QCPRIDE:E0ADF25D530854F3543D`
  - Definition: Maximum positive first-precursor isolation width in Th (m/z), MS2.

- **`IsolationWidth_MS2_Count`** — diagnostic; scalar; applies to ms2
  - CV: `QCPRIDE:0F163ECB26874F0C65C7`
  - Definition: Count of first-precursor isolation widths with a finite positive value, MS2.

- **`IsolationWidth_MS2_FractionLe15`** — extended; scalar; applies to ms2
  - CV: `QCPRIDE:C093A5D8209394150C2B`
  - Definition: Fraction of recorded positive first-precursor isolation widths <=15 Th (m/z), MS2.

- **`IsolationWidth_MS2_FractionGe15`** — extended; scalar; applies to ms2
  - CV: `QCPRIDE:DE0B3E119A726F58C6C1`
  - Definition: Fraction of recorded positive first-precursor isolation widths >=15 Th (m/z), MS2.

- **`IsolationPrecursorMz_MS2_Count`** — extended; scalar; applies to ms2
  - CV: `QCPRIDE:5C5E15B79F1D0F5E5C1F`
  - Definition: Count of finite positive first-precursor m/z values paired with positive isolation widths, MS2.

- **`IsolationPrecursorMz_MS2_UniqueCount`** — extended; scalar; applies to ms2
  - CV: `QCPRIDE:C98D3B2F6D4A5E4D3E1B`
  - Definition: Count of unique first-precursor m/z values after rounding to 0.1 Th, MS2.

- **`IsolationPrecursorMz_MS2_UniqueFraction`** — diagnostic; scalar; applies to ms2
  - CV: `QCPRIDE:B7B2C4A7F8F6D1A4C6E9`
  - Definition: Unique first-precursor m/z target count divided by paired isolation-width observations, after 0.1 Th rounding, MS2.

- **`IsolationPrecursorMz_MS2_MaxRepeatFraction`** — extended; scalar; applies to ms2
  - CV: `QCPRIDE:AAE8B64F7E3D2B6C4D9A`
  - Definition: Largest fraction of paired isolation-width observations sharing one rounded (0.1 Th) first-precursor m/z target, MS2.

- **`IsolationPrecursorMz_MS2_RepeatedFraction`** — extended; scalar; applies to ms2
  - CV: `QCPRIDE:8D2E4B9C7A1F6E5B3C0D`
  - Definition: Fraction of paired isolation-width observations whose rounded (0.1 Th) first-precursor m/z target occurs at least twice, MS2.

- **`MultiplePrecursors_MS2_Count`** — diagnostic; scalar; applies to ms2
  - CV: `QCPRIDE:82C9F41E70F7A9F93FCC`
  - Definition: Scans with multiple precursors; first precursor supplies precursor-distribution metrics, MS2.

- **`ExtentPrecursorIntensity_95over5_MS2`** — extended; scalar; applies to ms2
  - CV: `QCPRIDE:BD6FAF3B9FAECC3FF579`
  - Definition: 95th/5th percentile ratio of positive recorded first-precursor MS2 intensities; ID-free proxy, no fragment-TIC fallback.

- **`MS2_PrecursorCharge_Fractions`** — core; table; applies to ms2
  - CV: `QCPRIDE:EF3BCB4B9E83C77AAC5B`
  - Definition: Counts and fractions of first-precursor charges, denominator all MS2 scans; missing or nonpositive charge is unknown.

- **`ChargeMean`** — extended; scalar; applies to ms2
  - CV: `QCPRIDE:E303A6C9571AC1DD9AD6`
  - Definition: Mean positive first-precursor MS2 charge; unknown charges excluded.

- **`ChargeMedian`** — diagnostic; scalar; applies to ms2
  - CV: `QCPRIDE:3E06B5777AC79EC2AAC6`
  - Definition: Median positive first-precursor MS2 charge; unknown charges excluded.

- **`ChargeMin`** — extended; scalar; applies to ms2
  - CV: `QCPRIDE:87C7AAEB9CFD9A6651E1`
  - Definition: Minimum positive first-precursor MS2 charge.

- **`ChargeMax`** — extended; scalar; applies to ms2
  - CV: `QCPRIDE:7BEF7C9E6BFC65BFAB6D`
  - Definition: Maximum positive first-precursor MS2 charge.

- **`ChargeRatio_3over2`** — extended; scalar; applies to ms2
  - CV: `QCPRIDE:3263CF93B552C37EE7C0`
  - Definition: Count of MS2 charge 3 divided by charge 2; unavailable unless both are present (rawQC convention).

- **`ChargeRatio_4over2`** — extended; scalar; applies to ms2
  - CV: `QCPRIDE:D916EC99E4BD5371D2AE`
  - Definition: Count of MS2 charge 4 divided by charge 2; unavailable unless both are present (rawQC convention).

:::

### Scan timing (12)

| Registry metric | Tier | Applies to |
| --- | --- | --- |
| `RtRange_MS1` | diagnostic | ms1 |
| `ScanRate_MS1` | core | ms1 |
| `FastestFrequency_MS1` | extended | ms1 |
| `RT_MS1_Quantiles` | extended | ms1 |
| `RT_MS1_IQR` | extended | ms1 |
| `RT_MS1_IQRRate` | extended | ms1 |
| `RtRange_MS2` | diagnostic | ms2 |
| `ScanRate_MS2` | core | ms2 |
| `FastestFrequency_MS2` | extended | ms2 |
| `RT_MS2_Quantiles` | extended | ms2 |
| `RT_MS2_IQR` | extended | ms2 |
| `RT_MS2_IQRRate` | extended | ms2 |

:::{dropdown} Full definitions and CV accessions — Scan timing

- **`RtRange_MS1`** — diagnostic; tuple; applies to ms1
  - CV: `QCPRIDE:77A5ED9E8CB8C7F8CD4F`
  - Definition: Minimum and maximum finite RT in seconds, MS1.

- **`ScanRate_MS1`** — core; scalar; applies to ms1
  - CV: `QCPRIDE:BC47FFC025E342B4AD8D`
  - Definition: Count of finite-RT scans divided by that level's RT span, in scans per minute, MS1.

- **`FastestFrequency_MS1`** — extended; scalar; applies to ms1
  - CV: `MS:4000065`
  - Definition: Maximum scan count in a closed 60-second window divided by 60, including short runs.

- **`RT_MS1_Quantiles`** — extended; tuple; applies to ms1
  - CV: `MS:4000184`
  - Definition: Four consecutive scan-RT quartile widths divided by the MS1 time span.

- **`RT_MS1_IQR`** — extended; scalar; applies to ms1
  - CV: `QCPRIDE:1FBA04CFFBA1AF991EC1`
  - Definition: 75th minus 25th percentile of finite retention times, in seconds, MS1.

- **`RT_MS1_IQRRate`** — extended; scalar; applies to ms1
  - CV: `QCPRIDE:502F90C949DF08AF0C6A`
  - Definition: Count of finite-RT scans in the closed RT interquartile interval divided by its width in seconds, MS1; ID-free proxy.

- **`RtRange_MS2`** — diagnostic; tuple; applies to ms2
  - CV: `QCPRIDE:55D96D3835BFEE5AD7AB`
  - Definition: Minimum and maximum finite RT in seconds, MS2.

- **`ScanRate_MS2`** — core; scalar; applies to ms2
  - CV: `QCPRIDE:D3718560FE9B29568175`
  - Definition: Count of finite-RT scans divided by that level's RT span, in scans per minute, MS2.

- **`FastestFrequency_MS2`** — extended; scalar; applies to ms2
  - CV: `MS:4000066`
  - Definition: Maximum scan count in a closed 60-second window divided by 60, including short runs.

- **`RT_MS2_Quantiles`** — extended; tuple; applies to ms2
  - CV: `MS:4000185`
  - Definition: Four consecutive scan-RT quartile widths divided by the MS2 time span.

- **`RT_MS2_IQR`** — extended; scalar; applies to ms2
  - CV: `QCPRIDE:56E9A52714B2A63AB701`
  - Definition: 75th minus 25th percentile of finite retention times, in seconds, MS2.

- **`RT_MS2_IQRRate`** — extended; scalar; applies to ms2
  - CV: `QCPRIDE:DC206A3283881BDD1E18`
  - Definition: Count of finite-RT scans in the closed RT interquartile interval divided by its width in seconds, MS2; ID-free proxy.

:::

### Spectra and peak summaries (15)

| Registry metric | Tier | Applies to |
| --- | --- | --- |
| `NumberOfMSLevels` | core | all |
| `NumberOfSpectra` | core | all |
| `NumberOfSpectralPeaks` | diagnostic | all |
| `NumberOfSpectra_MS1` | core | ms1 |
| `PeakDensity_MS1_Quantiles` | core | ms1 |
| `PeakTypes_MS1` | diagnostic | ms1 |
| `EstimatedPeakTypes_MS1` | extended | ms1 |
| `Polarity_MS1` | diagnostic | ms1 |
| `InvalidPeakCount_MS1` | diagnostic | ms1 |
| `NumberOfSpectra_MS2` | core | ms2 |
| `PeakDensity_MS2_Quantiles` | core | ms2 |
| `PeakTypes_MS2` | diagnostic | ms2 |
| `EstimatedPeakTypes_MS2` | extended | ms2 |
| `Polarity_MS2` | diagnostic | ms2 |
| `InvalidPeakCount_MS2` | diagnostic | ms2 |

:::{dropdown} Full definitions and CV accessions — Spectra and peak summaries

- **`NumberOfMSLevels`** — core; scalar; applies to all
  - CV: `QCPRIDE:C782CFB84B18715F8D67`
  - Definition: Number of distinct observed MS levels.

- **`NumberOfSpectra`** — core; scalar; applies to all
  - CV: `QCPRIDE:8418755F449095444376`
  - Definition: Total spectrum count over all levels.

- **`NumberOfSpectralPeaks`** — diagnostic; scalar; applies to all
  - CV: `QCPRIDE:EED0065AF36234685A04`
  - Definition: Total valid stored peak data points over all levels; not chromatographic features.

- **`NumberOfSpectra_MS1`** — core; scalar; applies to ms1
  - CV: `MS:4000059`
  - Definition: Count of all MS1 spectra, including empty scans.

- **`PeakDensity_MS1_Quantiles`** — core; tuple; applies to ms1
  - CV: `MS:4000061`
  - Definition: 25th, 50th, 75th percentiles of valid peaks per MS1 spectrum; empty scans included.

- **`PeakTypes_MS1`** — diagnostic; table; applies to ms1
  - CV: `QCPRIDE:658F0AF60FACA5FDD19E`
  - Definition: Counts of centroid, profile and unknown scan representations, MS1.

- **`EstimatedPeakTypes_MS1`** — extended; table; applies to ms1
  - CV: `QCPRIDE:499B085DD486DCBECCC0`
  - Definition: Optional OpenMS PeakTypeEstimator counts for scans with more than ten peaks, MS1; empty if estimation was not requested or eligible.

- **`Polarity_MS1`** — diagnostic; table; applies to ms1
  - CV: `QCPRIDE:5D7B6DFCAC3E42EAF25F`
  - Definition: Counts by recorded polarity, including unknown, MS1.

- **`InvalidPeakCount_MS1`** — diagnostic; scalar; applies to ms1
  - CV: `QCPRIDE:4ADF3C521FC5FF049772`
  - Definition: Count of excluded peaks with non-finite/nonpositive m/z or non-finite/negative intensity, MS1.

- **`NumberOfSpectra_MS2`** — core; scalar; applies to ms2
  - CV: `MS:4000060`
  - Definition: Count of all MS2 spectra, including empty scans.

- **`PeakDensity_MS2_Quantiles`** — core; tuple; applies to ms2
  - CV: `MS:4000062`
  - Definition: 25th, 50th, 75th percentiles of valid peaks per MS2 spectrum; empty scans included.

- **`PeakTypes_MS2`** — diagnostic; table; applies to ms2
  - CV: `QCPRIDE:02BD471D53BE3411913F`
  - Definition: Counts of centroid, profile and unknown scan representations, MS2.

- **`EstimatedPeakTypes_MS2`** — extended; table; applies to ms2
  - CV: `QCPRIDE:201C9F534D91EC37FC1D`
  - Definition: Optional OpenMS PeakTypeEstimator counts for scans with more than ten peaks, MS2; empty if estimation was not requested or eligible.

- **`Polarity_MS2`** — diagnostic; table; applies to ms2
  - CV: `QCPRIDE:6A1C64BAB5CCF09A5DE3`
  - Definition: Counts by recorded polarity, including unknown, MS2.

- **`InvalidPeakCount_MS2`** — diagnostic; scalar; applies to ms2
  - CV: `QCPRIDE:0F8B65A410A572FEE4D7`
  - Definition: Count of excluded peaks with non-finite/nonpositive m/z or non-finite/negative intensity, MS2.

:::

## Provenance and caveats

The complete set of calculated QC metrics remains in mzQC per run
and in losslessly consolidated accession mzQC. A website or
MongoDB projection may select **stored values only**, never generate
report-only scientific measurements. Grouped spectrum counts, MS2/MS3
representation, exploratory reporter-ion evidence and denominators
must retain their original definitions. A missing measurement is
not a zero value, low score, or proof that an experimental method
was absent.

`NumberOfSpectra_MS1`, `NumberOfSpectra_MS2`,
`MS2_PrecursorCharge_Fractions`, `ChromatographyDuration`,
and isolation-width statistics remain in the scientific registry.
Ion-source CV metadata are observed only if the reader supplies it.
RT/duration do not reveal LC separation chemistry. Label/channel
annotations require verified submitter/SDRF evidence for writeback.

See also [QC workflow](qc-workflow.md) and [measurement evidence](measurement-evidence.md).
