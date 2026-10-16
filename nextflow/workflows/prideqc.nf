#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { PRIDEQC_MANIFEST } from '../modules/local/manifest/main'
include { PRIDEQC_PREPARE }  from '../modules/local/prepare/main'
include { PRIDEQC_ANALYZE }  from '../modules/local/analyze/main'
include { PRIDEQC_FINALIZE } from '../modules/local/finalize/main'

workflow PRIDEQC_WORKFLOW {
    if (!params.input || !params.report_manifest || !params.outdir || !params.sif) {
        error 'Require --input (file manifest), --report_manifest, --outdir and --sif'
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
    PRIDEQC_PREPARE(archive_tasks,
        file("${projectDir}/bin/prepare_archive.py"),
        file("${projectDir}/bin/archive_units.py"))

    // Split each archive into one independent Slurm QC task per discovered run.
    def logical_units = PRIDEQC_PREPARE.out.units.flatMap { accession, _expected_archives, _task_id, unit_dirs, _prepare_status ->
        def units = unit_dirs instanceof List ? unit_dirs : [unit_dirs]
        units.collect { unit -> tuple(accession, unit) }
    }
    PRIDEQC_ANALYZE(logical_units, file("${projectDir}/bin/analyze_unit.py"))

    // Dynamic accession barrier: sum discovered units only after all expected
    // archive tasks finished. Failure sentinels still count, so nothing vanishes.
    def expected_counts = PRIDEQC_PREPARE.out.units
        .map { accession, expected_archives, _task_id, unit_dirs, _prepare_status ->
            def units = unit_dirs instanceof List ? unit_dirs : [unit_dirs]
            tuple(groupKey(accession, expected_archives), units.size())
        }
        .groupTuple()
        .map { key, counts -> tuple(key.getGroupTarget(), counts.sum() as Integer) }

    // Join each run receipt to the accession's actual discovered-run count.
    // Nextflow groups a study as soon as that many receipts have arrived.
    def grouped = PRIDEQC_ANALYZE.out.unit_result
        .combine(expected_counts, by: 0)
        .map { accession, unit_result, expected_units ->
            tuple(groupKey(accession, expected_units), unit_result)
        }
        .groupTuple()
        .map { key, run_dirs ->
            tuple(key.getGroupTarget(), run_dirs.size() as Integer, run_dirs)
        }
    PRIDEQC_FINALIZE(grouped, file("${projectDir}/bin/finalize_accession.py"))
}
