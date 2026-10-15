#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { PRIDEQC_WORKFLOW } from './workflows/prideqc'

workflow {
    PRIDEQC_WORKFLOW()
}
