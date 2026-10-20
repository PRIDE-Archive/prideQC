#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { PRIDEQC_MANIFEST } from '../modules/local/manifest/main'
include { PRIDEQC_SOURCE_QC } from '../modules/local/source_qc/main'
include { PRIDEQC_FINALIZE } from '../modules/local/finalize/main'

workflow PRIDEQC_WORKFLOW {
    if (!params.input || !params.report_manifest || !params.outdir || !params.sif) {
        error 'Require --input (file manifest), --report_manifest, --outdir and --sif'
    }
    // The current pilot cannot independently prove that an entire accession
    // inventory is represented, particularly with vendor/mzML alternatives.
    // A CLI/operator Boolean assertion must never enable study publication.
    if (params.allow_complete_study) {
        error 'Full-study merging is disabled until a verified inventory attestation gate is implemented'
    }
    if (!(params.protocol in ['ftp','aspera','s3','globus','fire'])) {
        error "Invalid --protocol: ${params.protocol}"
    }
    if (params.accession && !(params.accession ==~ /PXD[0-9]{6,}/)) {
        error "Invalid --accession: ${params.accession}"
    }
    if (!(params.outdir as String).startsWith('/')) {
        error '--outdir must be an absolute persistent path for published provenance'
    }
    if (!file(params.sif).exists()) {
        error "SIF does not exist: ${params.sif}"
    }

    PRIDEQC_MANIFEST(
        channel.fromPath(params.input, checkIfExists: true),
        channel.fromPath(params.report_manifest, checkIfExists: true),
        file("${projectDir}/bin/plan_manifest.py")
    )

    // The frozen manifest counts SOURCE archives, not extracted logical MS runs.
    def archive_tasks = PRIDEQC_MANIFEST.out.tasks.flatten().map { task_file ->
        def record = new groovy.json.JsonSlurper().parse(task_file)
        tuple(record.accession, record.expected_runs as Integer,
              record.task_id.toString(), task_file)
    }
    // A single Slurm scratch task downloads, expands and analyzes one logical
    // source run. Multi-run archives emit a failure receipt in this pilot lane.
    // Only QC outputs/receipts, never extracted RAW, travel between processes.
    PRIDEQC_SOURCE_QC(archive_tasks,
        file("${projectDir}/bin/source_qc.py"),
        file("${projectDir}/bin/prepare_archive.py"),
        file("${projectDir}/bin/archive_units.py"),
        file("${projectDir}/bin/analyze_unit.py"))

    def units = PRIDEQC_SOURCE_QC.out.completed.flatMap { accession, _expected_archives, _task_id, unit_dirs, _status ->
        def results = unit_dirs instanceof List ? unit_dirs : [unit_dirs]
        results.collect { result_dir -> tuple(accession, result_dir) }
    }
    def expected_counts = PRIDEQC_SOURCE_QC.out.completed
        .map { accession, expected_archives, _task_id, unit_dirs, _status ->
            def results = unit_dirs instanceof List ? unit_dirs : [unit_dirs]
            tuple(groupKey(accession, expected_archives), results.size())
        }
        .groupTuple()
        .map { key, counts -> tuple(key.getGroupTarget(), counts.sum() as Integer) }
    def grouped = units
        .combine(expected_counts, by: 0)
        .map { accession, result_dir, expected_units ->
            tuple(groupKey(accession, expected_units), result_dir)
        }
        .groupTuple()
        .map { key, run_dirs ->
            tuple(key.getGroupTarget(), run_dirs.size() as Integer, run_dirs)
        }
    PRIDEQC_FINALIZE(grouped, file("${projectDir}/bin/finalize_accession.py"))
}
