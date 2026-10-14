# Experimental TMT evidence from profile spectra

This is an **independent, explicitly opt-in** profile-spectrum supplement to the validated centroid-only `--tmt-reporter-evidence` feature. It uses the OpenMS `PeakPickerHiRes` class on an ephemeral m/z window spanning runtime pyOpenMS TMT reporter channels with 1 Da margin. Native Thermo RAW spectra are not modified, written back, or converted to new files.

```bash
uv run --active --no-project python -m prideqc analyze \
  /path/to/input.raw -o /path/to/new-output \
  --tmt-reporter-evidence --tmt-profile-reporter-evidence
```

The two switches produce **two distinct annotations**: `tmt_reporter_evidence` for *native centroid* data and `tmt_profile_reporter_evidence` for *ephemerally peak-picked, explicitly profile* MS2/MS3 data. They have separate scan counts and denominator definitions and must **not be summed**. The new annotation is exported to summary JSON and a local `QCPRIDE` mzQC metric. Both are exploratory, do not infer an exact multiplexing kit or biological sample/channel occupancy, and cannot populate SDRF `comment[label]`.

Profile processing requires valid positive, sorted m/z values and nonnegative intensities. Failed OpenMS peak picking, invalid sample grids, and absent reporter-region sampling are recorded as separate **inconclusive** categories, not as negative TMT calls. A scan is considered evaluable only if at least one runtime catalog channel has three local raw profile samples within ±0.02 Da. Reporter matches use the same default 10 ppm, 0.1%-of-**full-original-scan** base peak and three-channel screen as the centroid collector. The peak picker is applied to reporter-region m/z samples only and results stay ephemeral.

Per-channel `profile_spectra_with_local_mz_sampling` is reported separately from `spectra_with_hit`; channel-conditional hit fractions are null when the m/z region was not sampled. The screening values are neither isotopic impurity-corrected nor calibrated FDR results. Profile-derived and native centroid signals may have different quantitative biases, and the implementation is not yet experimentally qualified for low-resolution data, coalescent peaks, or all vendor profile representations.

## First real-data gate

Existing inputs (no downloads):
- `~/Documents/datasets/prideqc_tmt_pilot_v1/PXD000001/TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.raw` (profile MS2 pilot previously inconclusive)
- `~/Documents/datasets/prideqc_tmt_pilot_v1/PXD005180/Fusion_KMTMT10RPpH10_P02X.raw` (profile MS3 pilot previously inconclusive)

Begin with **PXD000001** in a fresh output directory, with the native RAW-capable `pyopenms==3.6.0.dev20260910` environment. Check summary JSON and mzQC evidence parity, inspect MS2 `profile_spectra_evaluable`, `profile_peak_picker_failures`, and method/channel coverage. Then test PXD005180 MS3. Do not expect the earlier profile-skipping standalone probe to provide a comparable reporter-positive baseline. Evaluate profile peak-picking against independent protocol labels and negatives before treating results as anything beyond experimental evidence. A crash, zero support or high failure rate does not validate a negative TMT call.

## Reproducibility

The module has offline mock-OpenMS tests for profile spectra, centroid/profile isolation, spectrum-pass accounting, out-of-range and invalid scans, picker exceptions, flat profiles, full-scan base-peak thresholds, and mzQC serialization. Real Thermo RAW profile validation must be run by the user; do not claim it has passed yet. No changes to the validated centroid collector are part of this patch.
