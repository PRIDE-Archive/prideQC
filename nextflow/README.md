# PRIDEQC Nextflow DSL2 — archive-aware pilot

Entry point: `main.nf`; version 2 adapters under `bin/` and modules under
`modules/local/`. Requires a pinned SIF with prideqc + pridepy, Nextflow DSL2,
and a frozen **file-task** and **accession-report** manifest.

## Pipeline

1. `PRIDEQC_MANIFEST`: fail-closed validation of frozen source-archive tasks.
2. `PRIDEQC_PREPARE`: one pridepy download for each listed source file, then
   safely expand a ZIP or TAR.GZ to one *logical MS run* per unit directory.
   Supported inputs: `.raw`, `.raw.zip`, `.raw.gz`, `.mzML`, `.mzML.gz`,
   `.d.zip`, `.d`, `.tar.gz`, `.tgz` and generic bundles ending in `.zip`.
   ZIPs may contain several runs; `.d` folders are never split apart.
   Nested ZIP/TAR.GZ supported to bounded depth. Unsupported/corrupt archives
   produce failure sentinels and cannot lead to a complete study.
3. `PRIDEQC_ANALYZE`: one isolated QC job per logical run. Every reached run
   emits `result.json` with success or failure, and a per-run mzQC only on success.
4. `PRIDEQC_FINALIZE`: dynamically wait for the discovered number of unit
   receipts across the expected source-archive tasks. Write `_audit` JSON/TSV;
   **only** complete success permits a lossless study `_study/*.mzQC` merge.

`expected_file_tasks` counts original PRIDE archive entries, whereas
`expected_logical_runs` is calculated after validated extraction. These must
never be equated when a source ZIP contains multiple RAW files.

## Current limitations and production gates

* Requires `nextflow lint .`, local/remote profile resolution, Nextflow graph
  preview, a small real-archive end-to-end test, and Codon Slurm execution.
* Python offline tests do **not** prove FIRE availability, vendor `.d` reading,
  SIF completeness, resource limits, or Nextflow output/publish semantics.
* Supported v2 source types are metadata-only; SDRF-aligned per-member mapping,
  WIFF/Waters readers and arbitrary archive formats require separate work.
* Archive extraction is bounded to 500 GiB uncompressed per source by default;
  large submissions need a per-cluster capacity policy, not a silent override.
* A scheduler/node crash may prevent a per-unit receipt; Nextflow will report
  task failure and the study must be treated as incomplete, not as a clean run.
* Nextflow workflow completion must **not** be interpreted as accession QC
  success. Check `_audit/accession-status.json` (`status == complete`).
* MongoDB write is intentionally disabled. Export only validated complete mzQC
  through the separately reviewed `export_mongo_ndjson.py` artifact.
* PRIDE team-specific Codon instructions are maintained privately; do not
  publish internal filesystem paths, passwords or FIRE credentials.
