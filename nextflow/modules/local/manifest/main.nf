process PRIDEQC_MANIFEST {
    tag 'manifest-validation'
    label 'process_low'
    container params.sif
    input:
    path file_manifest
    path report_manifest
    output:
    path 'plan/tasks/*.json', emit: tasks
    path 'plan/plan.json', emit: plan
    shell:
    '''
    plan_manifest.py \
      --manifest !{file_manifest} \
      --report-manifest !{report_manifest} \
      --outdir plan \
      !{params.accession ? '--accession '+params.accession : ''}
    '''
}
