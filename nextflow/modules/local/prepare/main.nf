process PRIDEQC_PREPARE {
    tag "${accession}:archive-${task_id}"
    label 'process_medium'
    container params.sif
    scratch true
    publishDir { "${params.outdir}/${accession}/_archive_status" },
        mode: 'copy', pattern: 'prepared/prepare-status.json', overwrite: true,
        saveAs: { _filename -> "${task_id}.prepare-status.json" }
    input:
    tuple val(accession), val(expected_archives), val(task_id), path(task_json)
    path prepare_script
    path archive_script
    output:
    tuple val(accession), val(expected_archives), val(task_id),
        path('prepared/units/*'), path('prepared/prepare-status.json'), emit: units
    script:
    """
    export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
    python3 ${prepare_script} --task ${task_json} --protocol ${params.protocol} \
      --outdir prepared ${params.checksum_check ? '' : '--no-checksum-check'}
    """
}
