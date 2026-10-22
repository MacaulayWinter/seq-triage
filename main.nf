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

// Empty cells arrive as [] (D22); case and spaces do not make two
// declarations differ.
def declaredValue(v) {
    return (v == null || v == []) ? null : v.toString().trim().toLowerCase()
}

// FASTQ files in a fastq_dir folder (D75), sorted by name.
def dirFastqs(dir) {
    def all   = dir.listFiles()
    def found = all.findAll { f -> f.name ==~ /.+\.(f(ast)?q|fa|fasta|fna)(\.(gz|bz2|xz|zst))?$/ }.sort { f -> f.name }
    if (!found) {
        // Say what the folder IS when it is not reads yet (D79).
        def names = all.collect { f -> f.name }
        def what  = names.any { n -> n ==~ /.+\.(pod5|fast5)$/ } ? ' It holds raw nanopore signal (POD5/FAST5): basecall it first (e.g. dorado), then give the FASTQ or BAM.'
                  : names.contains('RunInfo.xml') || names.contains('Data') ? ' It looks like an Illumina run folder (BCL): demultiplex it first (bcl-convert or bcl2fastq), then give the FASTQ files.'
                  : names.any { n -> n.endsWith('.sra') } ? ' It holds SRA files: convert them first (fasterq-dump).'
                  : ''
        error("fastq_dir ${dir} holds no FASTQ or FASTA read files (.fastq/.fq, .fa/.fasta/.fna, plain or compressed).${what}")
    }
    return found
}

// One sample from its samplesheet rows: [ meta, [ files ] ]. Files are in
// lane order, mates together (L1 R1, L1 R2, L2 R1, ...), as nf-core's
// CAT_FASTQ expects. meta gains n_units only when there is more than one
// lane or file, so a one-file sample's meta (and its cache) is unchanged.
def sampleFromRows(List rs) {
    def meta  = rs[0][0].collectEntries { k, v -> [ k, v == [] ? null : v ] }
    def first = rs[0]
    if (first[3]) {
        return [ meta + [ single_end: false, input_format: 'bam' ], [ first[3] ] ]
    }
    def single_end = !first[2]
    def units = rs.collectMany { _m, r1, r2, _bam, dir, _target, _gsize -> dir ? dirFastqs(dir).collect { f -> [ f ] } : [ r2 ? [ r1, r2 ] : [ r1 ] ] }
    def extra = units.size() > 1 ? [ n_units: units.size() ] : [:]
    return [ meta + [ single_end: single_end ] + extra, units.flatten() ]
}

workflow {
    // --help / --helpFull: nf-schema has already printed the help (generated
    // from nextflow_schema.json, D53). Stop here, so a help request is not
    // also reported as a run with no --input.
    if (params.containsKey('help') || params.containsKey('helpFull')) {
        exit(0)
    }

    // Stop before any work if a parameter is missing, mistyped or points at
    // a file that does not exist. Rules come from nextflow_schema.json, which
    // also validates the samplesheet against assets/schema_input.json.
    validateParameters()

    // Print the parameters that differ from their defaults, so every run log
    // records exactly how the pipeline was invoked.
    log.info paramsSummaryLog(workflow)

    // One [meta, fastq_1, fastq_2, bam] list per row. Schema fields marked
    // "meta" become keys of the meta map: id, declared_platform, declared_type.
    // A row gives FASTQ files or one BAM, never both (the schema checks it).
    def rows = samplesheetToList(params.input, "${projectDir}/assets/schema_input.json")
    def same_mates = rows.findAll { _meta, r1, r2, _bam, _dir, _target, _gsize -> r2 && r1 == r2 }
    if (same_mates) {
        error("fastq_1 and fastq_2 are the same file for sample(s): " +
              same_mates.collect { meta, _r1, _r2, _bam, _dir, _target, _gsize -> meta.id }.join(', '))
    }

    // Several rows with the same sample ID are one sample whose reads are in
    // several files: lanes, top-up runs (D75). They must agree on everything
    // but the files, or the run stops here, naming the sample.
    def problems = []
    rows.groupBy { row -> row[0].id }.each { id, rs ->
        def files = rs.collectMany { _meta, r1, r2, bam, dir, _target, _gsize -> [ r1, r2, bam, dir ].findAll { f -> f } }
        def twice = files.findAll { f -> files.count(f) > 1 }.unique()
        if (twice) {
            problems << "${id}: the same file is listed twice (${twice.collect { f -> f.name }.join(', ')})"
        }
        if (rs.size() == 1) {
            return
        }
        if (rs.any { row -> row[3] }) {
            problems << "${id}: a BAM, SAM or CRAM file must be the only row of its sample"
        } else if (rs.collect { row -> row[2] ? 'paired' : 'single' }.unique().size() > 1) {
            problems << "${id}: some rows are paired (fastq_2 given) and some are not"
        }
        [ declared_platform: 'platform', declared_type: 'sample_type', declared_molecule: 'molecule' ].each { key, column ->
            def values = rs.collect { row -> declaredValue(row[0][key]) }.unique()
            if (values.size() > 1) {
                problems << "${id}: rows disagree on ${column} (${values.collect { v -> v ?: 'blank' }.join(', ')})"
            }
        }
        def targets = rs.collect { row -> declaredValue(row[5]) }.unique()
        if (targets.size() > 1) {
            problems << "${id}: rows disagree on target_taxid (${targets.collect { v -> v ?: 'blank' }.join(', ')})"
        }
        def sizes = rs.collect { row -> declaredValue(row[6]) }.unique()
        if (sizes.size() > 1) {
            problems << "${id}: rows disagree on genome_size (${sizes.collect { v -> v ?: 'blank' }.join(', ')})"
        }
    }
    if (problems) {
        error("Samplesheet rows of the same sample must agree (D75):\n  " + problems.join('\n  '))
    }

    ch_samples = channel.fromList(rows.groupBy { row -> row[0].id }.values() as List)
        .map { rs -> sampleFromRows(rs) }

    SEQ_TRIAGE(ch_samples, file(params.input), file(params.metric_guide), file(params.species_groups), file(params.genome_sizes))
}
