#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { PRIDEQC_MANIFEST } from '../modules/local/manifest/main'
include { PRIDEQC_ANALYZE }  from '../modules/local/analyze/main'
include { PRIDEQC_MERGE }    from '../modules/local/merge/main'

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
        error '--outdir must be an absolute persistent path for published mzQC OBO provenance'
    }
    if (!file(params.sif).exists()) {
        error "SIF does not exist: ${params.sif}"
    }

    PRIDEQC_MANIFEST(
        Channel.fromPath(params.input, checkIfExists: true),
        Channel.fromPath(params.report_manifest, checkIfExists: true)
    )

    def tasks = PRIDEQC_MANIFEST.out.tasks.flatten().map { task_file ->
        def record = new groovy.json.JsonSlurper().parse(task_file)
        def slug = record.archive_file.replaceAll(/[^A-Za-z0-9._-]/, '_')
        tuple(record.accession, record.expected_runs as Integer,
              record.task_id.toString(), slug, task_file)
    }
    PRIDEQC_ANALYZE(tasks)

    def per_accession = PRIDEQC_ANALYZE.out.run
        .map { accession, expected, mzqc, obo ->
            tuple(groupKey("${accession}:${expected}", expected), mzqc, obo)
        }
        .groupTuple()
        .map { key, mzqcs, obos ->
            def parts = key.getGroupTarget().split(':')
            tuple(parts[0], parts[1].toInteger(), mzqcs, obos)
        }
    PRIDEQC_MERGE(per_accession)
}
