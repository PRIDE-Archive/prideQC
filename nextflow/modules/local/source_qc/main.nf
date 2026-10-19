// First-flight lane only: one source archive MUST contain exactly one logical run.
// No RAW-containing path is a declared output or publishDir pattern.
process PRIDEQC_SOURCE_QC {
    tag "${accession}:archive-${task_id}"
    label 'process_medium'
    container params.sif
    scratch true
    publishDir { "${params.outdir}/${accession}/_archive_status" },
        mode: 'copy', pattern: 'ready/prepare-status.json', overwrite: true,
        saveAs: { _filename -> "${task_id}.prepare-status.json" }
    publishDir { "${params.outdir}/${accession}/_run_outputs" },
        mode: 'copy', pattern: 'ready/units/*', overwrite: true
    input:
    tuple val(accession), val(expected_archives), val(task_id), path(task_json)
    path runner
    path prepare_script
    path archive_script
    path analyzer
    output:
    tuple val(accession), val(expected_archives), val(task_id),
        path('ready/units/*'), path('ready/prepare-status.json'), emit: completed
    script:
    """
    export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
    export TMPDIR="\$PWD/task-temp" XDG_CACHE_HOME="\$PWD/task-cache"
    mkdir -p "\$TMPDIR" "\$XDG_CACHE_HOME"
    python3 ${runner} --task ${task_json} --protocol ${params.protocol} \
        --scratch-dir node-local-raw --outdir ready \
        ${params.checksum_check ? '' : '--no-checksum-check'} \
        ${params.estimate_mass_error ? '--estimate-mass-error' : ''} \
        ${params.estimate_mass_shifts ? '--estimate-mass-shifts' : ''}
    """
}
