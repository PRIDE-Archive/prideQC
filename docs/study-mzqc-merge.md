# Accession-level mzQC consolidation

`prideqc analyze` writes one mzQC and companion OBO vocabulary for each successful raw file. The accession reporter can also publish a **packaged**, multi-run mzQC 1.0 document. It does **not** average, sum, or otherwise recompute run quality metrics.

## Output contract

```text
PXD041271/
  repository-metadata/             # Original per-run mzQC and OBO files
  _study/
    PXD041271.mzQC                 # runQualities: all eligible per-run entries
    PXD041271.obo                  # Union of compatible QCPRIDE terms
    PXD041271.manifest.json        # Members, checksums, expected and actual counts
  _multiqc/
    multiqc_report.html
    report-info.txt                # Contains study_merge_status
```

The original `runQualities` entries, source filenames, metadata, measurement values, annotations and units are preserved. The original per-run files are never modified. The study-level `controlledVocabularies` includes a union of consistent standard CVs and a **single combined local CV**. `setQualities` is deliberately absent; any cohort statistics require separately defined denominators and grouping semantics.

The merged `.obo` is a required sidecar for `QCPRIDE:` metrics. By default its URI is a local `file://` URI pointing at the output location. On Codon the report wrapper passes a host-visible URI (not the temporary Singularity mount path). Such URIs are **not portable off the host**: republishing the mzQC needs republishing the OBO and setting `--cv-uri` to the new retrievable URI. Existing input-file location URIs are retained as original acquisition provenance and are not guaranteed to be externally resolvable.

## Command-line use

```bash
prideqc mzqc merge \
  --input-dir results/PXD041271/repository-metadata \
  --output results/PXD041271/_study/PXD041271.mzQC \
  --accession PXD041271 \
  --expected-runs 67
```

The `--expected-runs` value must be derived from the file-task manifest, **not** the observed number of successful outputs. Without this flag the manifest status is `count_unverified`, never `complete`. The reporter passes `--expected-runs` from its accession task manifest, and runs the merger only after a successful MultiQC invocation and exact expected-file count match.

The merger rejects invalid input shape, missing OBO sidecars, non-1.0 mzQC, zero/multiple per-run `runQualities`, `setQualities` in sources, duplicate run labels or input-file basenames, mismatch of project accession when declared, unsupported QCPRIDE version, and conflicting standard CV or local OBO term definitions. It does not repair these inconsistencies automatically. Successful output is written atomically per file and the manifest is published last, so consumers can verify SHA-256 digests; repeated runs with the same inputs and intact outputs avoid rewriting.

**Consumer rule:** MultiQC discovery is limited to `repository-metadata/`; do not additionally ingest `_study/<PXD>.mzQC` in the same report. MultiQC therefore continues to parse the originals, while the study mzQC is available for standalone, downstream exchange. No pmultiqc multi-run parser compatibility is assumed until tested with the deployment-pinned version.

## Deployment caveats

- Do not create `_study/<PXD>.mzQC` under an accession while the **old recursive report wrapper** is still deployed; its next MultiQC run would rediscover the merged artifact. For an initial merge pilot, write to an external output directory, or deploy the narrowed-discovery reporter first.
- Updating the reporter script without rebuilding its prideQC SIF will cause `prideqc mzqc merge` to fail because the old image does not contain the new command. Rebuild and checksum-pin the SIF before scheduling the updated report jobs.
- This package has been unit-tested using synthetic 1/67-run mzQC examples and the bundled mzQC v1 JSON schema; a real `PXD041271` merge and Codon pmultiqc rerun have **not** been performed.
- An old `_study` artifact may remain after an incomplete later rerun; consumers must verify the manifest digests and current `report-info.txt` (`study_merge_status=complete`) rather than trusting the presence of the file alone.
- Runtime validation uses a fail-closed set of structural and semantic checks, not the full jsonschema library. The bundled schema is exercised in tests. Full schema and ontology validation of real source files is a recommended separate deployment gate.
