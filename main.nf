#!/usr/bin/env nextflow
/*
 * seq-triage: entry point
 *
 * main.nf checks the inputs and builds the sample channel. The analysis
 * wiring will live in workflows/seq_triage.nf (nf-core convention: the entry
 * point handles input, the workflow file handles logic).
 */

include { validateParameters; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'
include { SEQ_TRIAGE } from './workflows/seq_triage'

workflow {
    // Stop before any work if a parameter is missing, mistyped or points at
    // a file that does not exist. Rules come from nextflow_schema.json, which
    // also validates the samplesheet against assets/schema_input.json.
    validateParameters()

    // Print the parameters that differ from their defaults, so every run log
    // records exactly how the pipeline was invoked.
    log.info paramsSummaryLog(workflow)

    // One [meta, fastq_1, fastq_2] list per row. Schema fields marked "meta"
    // become keys of the meta map: id, declared_platform, declared_type.
    def rows = samplesheetToList(params.input, "${projectDir}/assets/schema_input.json")

    // R1 and R2 the same file: JSON Schema cannot compare two fields, so this
    // is checked here. It must be a hard error: seqkit would see identical
    // read counts and call the sample "paired" with high confidence.
    def same_mates = rows.findAll { _meta, r1, r2 -> r2 && r1 == r2 }
    if (same_mates) {
        error("fastq_1 and fastq_2 are the same file for sample(s): " +
              same_mates.collect { meta, _r1, _r2 -> meta.id }.join(', '))
    }

    ch_samples = channel.fromList(rows)
        .map { meta, r1, r2 ->
            // single_end comes from the samplesheet LAYOUT only. Whether the
            // reads agree with it is for the detection step to judge.
            def single_end = !r2
            // nf-schema returns a blank cell as [] (empty list). Convert to
            // null here, once, so downstream code has a single meaning of
            // "not declared". Non-blank values are kept exactly as written:
            // normalising "Nanopore" or "N/A" is the R screen's job, and it
            // reports that it did so (CLAUDE.md D7).
            def clean = meta.collectEntries { k, v -> [ k, v == [] ? null : v ] }
            [ clean + [ single_end: single_end ], single_end ? [ r1 ] : [ r1, r2 ] ]
        }

    SEQ_TRIAGE(ch_samples, file(params.input), file(params.metric_guide), file(params.species_groups))
}
