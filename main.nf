#!/usr/bin/env nextflow
/*
 * seq-triage: entry point
 *
 * main.nf only checks the parameters and hands over. The analysis wiring will
 * live in workflows/seq_triage.nf (nf-core convention: the entry point
 * handles input, the workflow file handles logic).
 */

include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'

workflow {
    // Stop before any work if a parameter is missing, mistyped or points at
    // a file that does not exist. Rules come from nextflow_schema.json.
    validateParameters()

    // Print the parameters that differ from their defaults, so every run log
    // records exactly how the pipeline was invoked.
    log.info paramsSummaryLog(workflow)

    log.info "Parameters valid. (Skeleton only: no analysis steps yet.)"
}
