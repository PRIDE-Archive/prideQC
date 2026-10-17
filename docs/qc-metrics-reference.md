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

| Display | Registry key | Shape | Applies when | Interpretation / limitation |
| --- | --- | --- | --- | --- |
| Run / MS levels | `NumberOfMSLevels` | scalar | all | Distinct MS levels recorded; not a quality score. |
| Total spectra | `NumberOfSpectra` | scalar | all | Counts every MS level including empty scans. |
| Run time (s) | `ChromatographyDuration` | scalar | all | RT span from spectra, not LC separation chemistry. |
| MS1 spectra | `NumberOfSpectra_MS1` | scalar | ms1 | Includes scans with no detected peaks. |
| Empty MS1 scans | `EmptyScans_MS1` | scalar | ms1 | Zero-TIC scans; interpret alongside the MS1 count. |
| MS1 peaks per scan (quartiles) | `PeakDensity_MS1_Quantiles` | tuple | ms1 | Valid stored peaks per scan, not identified molecular features. |
| Median MS1 TIC | `TIC_MS1_Median` | scalar | ms1 | Arbitrary instrument units; do not compare without calibration. |
| MS1 TIC variability (CV) | `TIC_MS1_CV` | scalar | ms1 | Undefined if the mean TIC is zero or too few scans are available. |
| MS1 scan rate | `ScanRate_MS1` | scalar | ms1 | Scans per minute over the MS1 RT span. |
| MS2 spectra | `NumberOfSpectra_MS2` | scalar | ms2 | Includes scans with no detected peaks. |
| Empty MS2 scans | `EmptyScans_MS2` | scalar | ms2 | Zero-TIC scans; interpret alongside the MS2 count. |
| MS2 peaks per scan (quartiles) | `PeakDensity_MS2_Quantiles` | tuple | ms2 | Valid stored peaks per scan, not identified molecular features. |
| Median MS2 TIC | `TIC_MS2_Median` | scalar | ms2 | Arbitrary instrument units; do not compare without calibration. |
| MS2 scan rate | `ScanRate_MS2` | scalar | ms2 | Scans per minute over the MS2 RT span. |
| Median MS2 isolation width | `IsolationWidth_MS2_Median` | scalar | ms2 | Requires positive recorded isolation offsets; absent is not zero. |
| MS1/MS2 scan ratio | `MS1_to_MS2_Ratio` | scalar | ms1_ms2 | Undefined when there are no MS2 scans; not a performance score. |
| MS2 scans per acquisition cycle | `AcquisitionCycle_MS2Count_Quantiles` | tuple | cycles | MS1-delimited cycles only; not applicable to MS1-free targeted runs. |
| MS2 precursor charge fractions | `MS2_PrecursorCharge_Fractions` | table | ms2 | Denominator is all MS2 scans; unknown charge is retained. |

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

### Acquisition cycles (12)

| Key | CV accession | Shape | Tier | Applicability | Definition |
| --- | --- | --- | --- | --- | --- |
| `AcquisitionCycle_Count` | `QCPRIDE:E5D44EF02DD0C112FCCD` | scalar | diagnostic | cycles | Count of MS1-delimited cycles containing at least one MS2 scan. MS2 scans before the first MS1 are ignored; a trailing cycle is included when it contains MS2 scans. |
| `AcquisitionCycle_MS2Count_Quantiles` | `QCPRIDE:CB0E144250E9279F1725` | tuple | core | cycles | 25th, 50th, 75th percentiles of MS2 scans per observed MS1-delimited cycle. |
| `AcquisitionCycle_MS2Count_Mode` | `QCPRIDE:38064FE0BC760341F5CE` | scalar | diagnostic | cycles | Most frequent MS2 scan count per observed MS1-delimited cycle; lower count wins ties. |
| `AcquisitionCycle_MS2Count_ModalFraction` | `QCPRIDE:F6EA9ED8D480A916D62F` | scalar | extended | cycles | Fraction of observed MS1-delimited cycles having the modal MS2 scan count. |
| `AcquisitionCycle_TargetCoverageFraction` | `QCPRIDE:1E7E0C9468F1B503778C` | scalar | diagnostic | cycles | Fraction of MS2 scans inside observed MS1-delimited cycles with a finite positive first-precursor m/z paired with a positive isolation width. |
| `AcquisitionCycle_TargetEligibleCount` | `QCPRIDE:172D73B0C593B7E3F376` | scalar | extended | cycles | Count of observed MS1-delimited cycles where every MS2 scan has target metadata eligible for target-set/order comparison. |
| `AcquisitionCycle_TargetEligibleFraction` | `QCPRIDE:0EA9DB6C9708D5F09586` | scalar | extended | cycles | Fraction of observed MS1-delimited cycles eligible for target-set/order comparison. |
| `AcquisitionCycle_UniqueTargetCount_Quantiles` | `QCPRIDE:EAE7AB20278D155810C1` | tuple | extended | cycles | 25th, 50th, 75th percentiles of unique first-precursor isolation targets per eligible cycle after 0.1 Th quantization. |
| `AcquisitionCycle_TargetSetDistinctCount` | `QCPRIDE:3DAD05F0A9B6418D738F` | scalar | extended | cycles | Number of distinct per-cycle isolation-target sets among eligible cycles after 0.1 Th quantization. |
| `AcquisitionCycle_TargetSetModalFraction` | `QCPRIDE:C1EE0988D7C3648932F1` | scalar | diagnostic | cycles | Fraction of eligible cycles whose quantized isolation-target set equals the modal set. |
| `AcquisitionCycle_TargetOrderDistinctCount` | `QCPRIDE:971C73C2545A365A4D7E` | scalar | extended | cycles | Number of distinct ordered per-cycle isolation-target sequences among eligible cycles after 0.1 Th quantization. |
| `AcquisitionCycle_TargetOrderModalFraction` | `QCPRIDE:4A5FB670C1ACFB1FBDB2` | scalar | extended | cycles | Fraction of eligible cycles whose quantized ordered isolation-target sequence equals the modal sequence. |

### Chromatography (5)

| Key | CV accession | Shape | Tier | Applicability | Definition |
| --- | --- | --- | --- | --- | --- |
| `ChromatographyDuration` | `MS:4000053` | scalar | core | all | Last minus first finite spectrum RT, over all MS levels. |
| `NumberOfChromatograms` | `QCPRIDE:3A03550E7F73CC8BC2E7` | scalar | diagnostic | chromatograms | Number of stored chromatograms. |
| `NumberOfChromatogramDataPoints` | `QCPRIDE:166BE0366D9F74CDD72A` | scalar | diagnostic | chromatograms | Total stored chromatogram array lengths, not chromatographic peaks. |
| `ChromatogramTypes` | `QCPRIDE:F6028B02FCB099FDD16A` | table | diagnostic | chromatograms | Counts grouped by OpenMS ChromatogramType enum value; engine version is recorded. |
| `ChromatogramRTRange` | `QCPRIDE:DA479CE356AAFB988099` | tuple | extended | chromatograms | Minimum and maximum finite chromatogram time, in seconds. |

### Ion intensity and TIC (22)

| Key | CV accession | Shape | Tier | Applicability | Definition |
| --- | --- | --- | --- | --- | --- |
| `BasePeak_All_Max` | `QCPRIDE:0DFA7AA5174921CAC620` | scalar | extended | all | Largest valid intensity across all spectra, including MS3 and higher; arbitrary intensity units. |
| `BasePeak_MS1_Mean` | `QCPRIDE:D452DD102903397E2819` | scalar | extended | ms1 | Mean of per-spectrum maximum valid intensity, excluding scans with no valid peaks, MS1. |
| `TIC_MS1_Median` | `QCPRIDE:C9C317E0A64A03006F7D` | scalar | core | ms1 | Median sum of valid peak intensities per scan, MS1. |
| `TIC_MS1_CV` | `QCPRIDE:00C4C4D556C5B163B0FE` | scalar | core | ms1 | Sample standard deviation (ddof=1) divided by mean TIC, MS1. |
| `TIC_MS1_Area` | `QCPRIDE:01859CE82042A7057DDF` | scalar | diagnostic | ms1 | Trapezoidal TIC integral against RT, stable time order, intensity times seconds, MS1. |
| `TIC_MS1_Area_RTQuantiles` | `QCPRIDE:59BA94235C553AF4FDC6` | tuple | extended | ms1 | Four exact integrals of piecewise-linear TIC split at scan-RT quartiles, intensity times seconds, MS1. |
| `TIC_MS1_QuartileLogRatios` | `QCPRIDE:604EE74E7854759E05FC` | tuple | extended | ms1 | Natural logs of TIC Q2/Q1, Q3/Q2, Q4/Q3 (Q4=max), MS1. |
| `TIC_MS1_ChangeQuartileLogRatios` | `QCPRIDE:B4CB91D742E1BE7A0B5E` | tuple | extended | ms1 | Natural logs of absolute adjacent TIC change Q2/Q1, Q3/Q2, Q4/Q3, MS1. |
| `TIC_MS1_SignalJump10x_Count` | `QCPRIDE:CC18B4EE12E3DE4D10CC` | scalar | diagnostic | ms1 | Adjacent sorted-TIC ratios >=10, with strictly positive previous TIC, MS1. |
| `TIC_MS1_SignalFall10x_Count` | `QCPRIDE:968FA11A94C7458E4F1C` | scalar | diagnostic | ms1 | Adjacent sorted-TIC ratios <=0.1, with strictly positive previous TIC; zero next TIC counts, MS1. |
| `MedianTIC_in_RT_MS1_IQR` | `QCPRIDE:384CE930503A31FD6181` | scalar | extended | ms1 | Median TIC in middle two scan-index quartile groups (recycle-and-sort assignment), MS1. ID-free proxy; not an identified-peptide metric. |
| `TIC_MS1_MedianInHalfRange` | `QCPRIDE:866968003763CDDE9CF5` | tuple | diagnostic | ms1 | Median TIC in the shortest window containing ceil(n/2) finite-RT MS1 scans; first window wins ties. ID-free proxy. |
| `RT_TIC_Quantiles` | `MS:4000183` | tuple | extended | ms1 | Four intervals at first scan-wise cumulative MS1 TIC quartile crossings, divided by the MS1 time span; this uses accumulation, not integration. |
| `BasePeak_MS2_Mean` | `QCPRIDE:25C2A0E5AF8883957F5E` | scalar | extended | ms2 | Mean of per-spectrum maximum valid intensity, excluding scans with no valid peaks, MS2. |
| `TIC_MS2_Median` | `QCPRIDE:D249949F6A3F6A3EBF3E` | scalar | core | ms2 | Median sum of valid peak intensities per scan, MS2. |
| `TIC_MS2_CV` | `QCPRIDE:8AF7F9131561A60EDB62` | scalar | diagnostic | ms2 | Sample standard deviation (ddof=1) divided by mean TIC, MS2. |
| `TIC_MS2_Area` | `QCPRIDE:9FD64B0DA7A59E8514E2` | scalar | diagnostic | ms2 | Trapezoidal TIC integral against RT, stable time order, intensity times seconds, MS2. |
| `TIC_MS2_Area_RTQuantiles` | `QCPRIDE:E191EB22593B10C83302` | tuple | extended | ms2 | Four exact integrals of piecewise-linear TIC split at scan-RT quartiles, intensity times seconds, MS2. |
| `TIC_MS2_QuartileLogRatios` | `QCPRIDE:80B75AE92157AF8E1090` | tuple | extended | ms2 | Natural logs of TIC Q2/Q1, Q3/Q2, Q4/Q3 (Q4=max), MS2. |
| `TIC_MS2_ChangeQuartileLogRatios` | `QCPRIDE:760CEDD06B184EE701FF` | tuple | extended | ms2 | Natural logs of absolute adjacent TIC change Q2/Q1, Q3/Q2, Q4/Q3, MS2. |
| `TIC_MS2_SignalJump10x_Count` | `QCPRIDE:4BE7DF39FCE21F9E3583` | scalar | diagnostic | ms2 | Adjacent sorted-TIC ratios >=10, with strictly positive previous TIC, MS2. |
| `TIC_MS2_SignalFall10x_Count` | `QCPRIDE:DCF0DBC8D6F83087ED6E` | scalar | diagnostic | ms2 | Adjacent sorted-TIC ratios <=0.1, with strictly positive previous TIC; zero next TIC counts, MS2. |

### Ion mobility / FAIMS (3)

| Key | CV accession | Shape | Tier | Applicability | Definition |
| --- | --- | --- | --- | --- | --- |
| `FAIMS_CV_Values` | `QCPRIDE:E69917C56F3DCCD6D0EE` | tuple | diagnostic | faims | Distinct finite FAIMS compensation voltages, sorted, in volts; not drift times. |
| `FAIMS_CV_Count` | `QCPRIDE:0FF9AAA5ADA010E8801C` | scalar | diagnostic | faims | Number of distinct finite FAIMS compensation voltages. |
| `FAIMS_CV_Range` | `QCPRIDE:3C721D04BCC72F40C1A0` | tuple | diagnostic | faims | Minimum and maximum finite FAIMS compensation voltages, in volts. |

### Other run-level QC (9)

| Key | CV accession | Shape | Tier | Applicability | Definition |
| --- | --- | --- | --- | --- | --- |
| `EmptyScans_MS1` | `QCPRIDE:6F59F45176612D8C6992` | scalar | core | ms1 | MS1 scans with no valid positive intensity (zero TIC). |
| `ObservedMzRange_MS1` | `QCPRIDE:2A7C646C00FBCB34F4F8` | tuple | diagnostic | ms1 | Min/max valid peak m/z; this is an observed range, not the instrument scan window, MS1. |
| `ScanWindow_MS1` | `QCPRIDE:05EB1DDF8CCEE51BE48C` | tuple | diagnostic | ms1 | Envelope of recorded scan-window bounds in Th (m/z), MS1. |
| `AvgCycleTime_MS1` | `QCPRIDE:2C6256D0D438E08FF913` | scalar | diagnostic | ms1 | Mean strictly positive adjacent sorted MS1 RT differences, in seconds. |
| `EmptyScans_MS2` | `QCPRIDE:3BE75D0FF4AD59C42832` | scalar | core | ms2 | MS2 scans with no valid positive intensity (zero TIC). |
| `ObservedMzRange_MS2` | `QCPRIDE:47E4D2D8B1F1F4D61C19` | tuple | diagnostic | ms2 | Min/max valid peak m/z; this is an observed range, not the instrument scan window, MS2. |
| `ScanWindow_MS2` | `QCPRIDE:BBD580CB1663521C4976` | tuple | diagnostic | ms2 | Envelope of recorded scan-window bounds in Th (m/z), MS2. |
| `MzRange_MS2` | `QCPRIDE:3432120DEC9FAE1739E7` | tuple | extended | ms2 | Range of positive finite first-precursor m/z, not fragment peak m/z, MS2. |
| `MS1_to_MS2_Ratio` | `QCPRIDE:3E6B36C0810EEDA24D7A` | scalar | core | ms1_ms2 | MS1 scan count divided by MS2 scan count; undefined when denominator is zero. |

### Precursors and isolation (26)

| Key | CV accession | Shape | Tier | Applicability | Definition |
| --- | --- | --- | --- | --- | --- |
| `PrecursorMz_MS2_Median` | `QCPRIDE:0DB185A646DB806298FC` | scalar | diagnostic | ms2 | Median positive finite first-precursor m/z, MS2. |
| `PrecursorIntensity_MS2_Quantiles` | `QCPRIDE:8FB2EFF7D173D0D4C3A2` | tuple | extended | ms2 | Quartiles of recorded positive first-precursor intensity, MS2; no fragment-TIC substitution. |
| `PrecursorIntensity_MS2_Mean` | `QCPRIDE:36D2085F9869CE667625` | scalar | extended | ms2 | Mean recorded positive first-precursor intensity, MS2. |
| `PrecursorIntensity_MS2_Sd` | `QCPRIDE:AA281BF4225DD99B5942` | scalar | extended | ms2 | Sample SD of recorded positive first-precursor intensity (ddof=1), MS2. |
| `PrecursorIntensity_MS2_MissingCount` | `QCPRIDE:37ABD88F90B17EAA841D` | scalar | diagnostic | ms2 | Scans lacking a positive finite first-precursor intensity, including the OpenMS unset default zero, MS2. |
| `IsolationWidth_MS2_Median` | `QCPRIDE:652589C44AB88BC96F40` | scalar | core | ms2 | Median positive first-precursor lower+upper isolation offsets in Th (m/z), MS2. |
| `IsolationWidth_MS2_Quantiles` | `QCPRIDE:482962258BA982FED6BB` | tuple | diagnostic | ms2 | 25th, 50th, 75th percentiles of positive first-precursor isolation widths in Th (m/z), MS2. |
| `IsolationWidth_MS2_Min` | `QCPRIDE:228CD009EF66CB6342FD` | scalar | extended | ms2 | Minimum positive first-precursor isolation width in Th (m/z), MS2. |
| `IsolationWidth_MS2_Max` | `QCPRIDE:E0ADF25D530854F3543D` | scalar | extended | ms2 | Maximum positive first-precursor isolation width in Th (m/z), MS2. |
| `IsolationWidth_MS2_Count` | `QCPRIDE:0F163ECB26874F0C65C7` | scalar | diagnostic | ms2 | Count of first-precursor isolation widths with a finite positive value, MS2. |
| `IsolationWidth_MS2_FractionLe15` | `QCPRIDE:C093A5D8209394150C2B` | scalar | extended | ms2 | Fraction of recorded positive first-precursor isolation widths <=15 Th (m/z), MS2. |
| `IsolationWidth_MS2_FractionGe15` | `QCPRIDE:DE0B3E119A726F58C6C1` | scalar | extended | ms2 | Fraction of recorded positive first-precursor isolation widths >=15 Th (m/z), MS2. |
| `IsolationPrecursorMz_MS2_Count` | `QCPRIDE:5C5E15B79F1D0F5E5C1F` | scalar | extended | ms2 | Count of finite positive first-precursor m/z values paired with positive isolation widths, MS2. |
| `IsolationPrecursorMz_MS2_UniqueCount` | `QCPRIDE:C98D3B2F6D4A5E4D3E1B` | scalar | extended | ms2 | Count of unique first-precursor m/z values after rounding to 0.1 Th, MS2. |
| `IsolationPrecursorMz_MS2_UniqueFraction` | `QCPRIDE:B7B2C4A7F8F6D1A4C6E9` | scalar | diagnostic | ms2 | Unique first-precursor m/z target count divided by paired isolation-width observations, after 0.1 Th rounding, MS2. |
| `IsolationPrecursorMz_MS2_MaxRepeatFraction` | `QCPRIDE:AAE8B64F7E3D2B6C4D9A` | scalar | extended | ms2 | Largest fraction of paired isolation-width observations sharing one rounded (0.1 Th) first-precursor m/z target, MS2. |
| `IsolationPrecursorMz_MS2_RepeatedFraction` | `QCPRIDE:8D2E4B9C7A1F6E5B3C0D` | scalar | extended | ms2 | Fraction of paired isolation-width observations whose rounded (0.1 Th) first-precursor m/z target occurs at least twice, MS2. |
| `MultiplePrecursors_MS2_Count` | `QCPRIDE:82C9F41E70F7A9F93FCC` | scalar | diagnostic | ms2 | Scans with multiple precursors; first precursor supplies precursor-distribution metrics, MS2. |
| `ExtentPrecursorIntensity_95over5_MS2` | `QCPRIDE:BD6FAF3B9FAECC3FF579` | scalar | extended | ms2 | 95th/5th percentile ratio of positive recorded first-precursor MS2 intensities; ID-free proxy, no fragment-TIC fallback. |
| `MS2_PrecursorCharge_Fractions` | `QCPRIDE:EF3BCB4B9E83C77AAC5B` | table | core | ms2 | Counts and fractions of first-precursor charges, denominator all MS2 scans; missing or nonpositive charge is unknown. |
| `ChargeMean` | `QCPRIDE:E303A6C9571AC1DD9AD6` | scalar | extended | ms2 | Mean positive first-precursor MS2 charge; unknown charges excluded. |
| `ChargeMedian` | `QCPRIDE:3E06B5777AC79EC2AAC6` | scalar | diagnostic | ms2 | Median positive first-precursor MS2 charge; unknown charges excluded. |
| `ChargeMin` | `QCPRIDE:87C7AAEB9CFD9A6651E1` | scalar | extended | ms2 | Minimum positive first-precursor MS2 charge. |
| `ChargeMax` | `QCPRIDE:7BEF7C9E6BFC65BFAB6D` | scalar | extended | ms2 | Maximum positive first-precursor MS2 charge. |
| `ChargeRatio_3over2` | `QCPRIDE:3263CF93B552C37EE7C0` | scalar | extended | ms2 | Count of MS2 charge 3 divided by charge 2; unavailable unless both are present (rawQC convention). |
| `ChargeRatio_4over2` | `QCPRIDE:D916EC99E4BD5371D2AE` | scalar | extended | ms2 | Count of MS2 charge 4 divided by charge 2; unavailable unless both are present (rawQC convention). |

### Scan timing (12)

| Key | CV accession | Shape | Tier | Applicability | Definition |
| --- | --- | --- | --- | --- | --- |
| `RtRange_MS1` | `QCPRIDE:77A5ED9E8CB8C7F8CD4F` | tuple | diagnostic | ms1 | Minimum and maximum finite RT in seconds, MS1. |
| `ScanRate_MS1` | `QCPRIDE:BC47FFC025E342B4AD8D` | scalar | core | ms1 | Count of finite-RT scans divided by that level's RT span, in scans per minute, MS1. |
| `FastestFrequency_MS1` | `MS:4000065` | scalar | extended | ms1 | Maximum scan count in a closed 60-second window divided by 60, including short runs. |
| `RT_MS1_Quantiles` | `MS:4000184` | tuple | extended | ms1 | Four consecutive scan-RT quartile widths divided by the MS1 time span. |
| `RT_MS1_IQR` | `QCPRIDE:1FBA04CFFBA1AF991EC1` | scalar | extended | ms1 | 75th minus 25th percentile of finite retention times, in seconds, MS1. |
| `RT_MS1_IQRRate` | `QCPRIDE:502F90C949DF08AF0C6A` | scalar | extended | ms1 | Count of finite-RT scans in the closed RT interquartile interval divided by its width in seconds, MS1; ID-free proxy. |
| `RtRange_MS2` | `QCPRIDE:55D96D3835BFEE5AD7AB` | tuple | diagnostic | ms2 | Minimum and maximum finite RT in seconds, MS2. |
| `ScanRate_MS2` | `QCPRIDE:D3718560FE9B29568175` | scalar | core | ms2 | Count of finite-RT scans divided by that level's RT span, in scans per minute, MS2. |
| `FastestFrequency_MS2` | `MS:4000066` | scalar | extended | ms2 | Maximum scan count in a closed 60-second window divided by 60, including short runs. |
| `RT_MS2_Quantiles` | `MS:4000185` | tuple | extended | ms2 | Four consecutive scan-RT quartile widths divided by the MS2 time span. |
| `RT_MS2_IQR` | `QCPRIDE:56E9A52714B2A63AB701` | scalar | extended | ms2 | 75th minus 25th percentile of finite retention times, in seconds, MS2. |
| `RT_MS2_IQRRate` | `QCPRIDE:DC206A3283881BDD1E18` | scalar | extended | ms2 | Count of finite-RT scans in the closed RT interquartile interval divided by its width in seconds, MS2; ID-free proxy. |

### Spectra and peak summaries (15)

| Key | CV accession | Shape | Tier | Applicability | Definition |
| --- | --- | --- | --- | --- | --- |
| `NumberOfMSLevels` | `QCPRIDE:C782CFB84B18715F8D67` | scalar | core | all | Number of distinct observed MS levels. |
| `NumberOfSpectra` | `QCPRIDE:8418755F449095444376` | scalar | core | all | Total spectrum count over all levels. |
| `NumberOfSpectralPeaks` | `QCPRIDE:EED0065AF36234685A04` | scalar | diagnostic | all | Total valid stored peak data points over all levels; not chromatographic features. |
| `NumberOfSpectra_MS1` | `MS:4000059` | scalar | core | ms1 | Count of all MS1 spectra, including empty scans. |
| `PeakDensity_MS1_Quantiles` | `MS:4000061` | tuple | core | ms1 | 25th, 50th, 75th percentiles of valid peaks per MS1 spectrum; empty scans included. |
| `PeakTypes_MS1` | `QCPRIDE:658F0AF60FACA5FDD19E` | table | diagnostic | ms1 | Counts of centroid, profile and unknown scan representations, MS1. |
| `EstimatedPeakTypes_MS1` | `QCPRIDE:499B085DD486DCBECCC0` | table | extended | ms1 | Optional OpenMS PeakTypeEstimator counts for scans with more than ten peaks, MS1; empty if estimation was not requested or eligible. |
| `Polarity_MS1` | `QCPRIDE:5D7B6DFCAC3E42EAF25F` | table | diagnostic | ms1 | Counts by recorded polarity, including unknown, MS1. |
| `InvalidPeakCount_MS1` | `QCPRIDE:4ADF3C521FC5FF049772` | scalar | diagnostic | ms1 | Count of excluded peaks with non-finite/nonpositive m/z or non-finite/negative intensity, MS1. |
| `NumberOfSpectra_MS2` | `MS:4000060` | scalar | core | ms2 | Count of all MS2 spectra, including empty scans. |
| `PeakDensity_MS2_Quantiles` | `MS:4000062` | tuple | core | ms2 | 25th, 50th, 75th percentiles of valid peaks per MS2 spectrum; empty scans included. |
| `PeakTypes_MS2` | `QCPRIDE:02BD471D53BE3411913F` | table | diagnostic | ms2 | Counts of centroid, profile and unknown scan representations, MS2. |
| `EstimatedPeakTypes_MS2` | `QCPRIDE:201C9F534D91EC37FC1D` | table | extended | ms2 | Optional OpenMS PeakTypeEstimator counts for scans with more than ten peaks, MS2; empty if estimation was not requested or eligible. |
| `Polarity_MS2` | `QCPRIDE:6A1C64BAB5CCF09A5DE3` | table | diagnostic | ms2 | Counts by recorded polarity, including unknown, MS2. |
| `InvalidPeakCount_MS2` | `QCPRIDE:0F8B65A410A572FEE4D7` | scalar | diagnostic | ms2 | Count of excluded peaks with non-finite/nonpositive m/z or non-finite/negative intensity, MS2. |

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
