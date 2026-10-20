/*
 * DETECT_PLATFORM: infer platform and layout (bin/detect_platform.R)
 *
 * Takes seqkit's whole-file statistics for the core call, plus the reads
 * themselves: the script reads only the first 2,000 records of each file,
 * for header format and quality values (CLAUDE.md D33).
 *
 * Per sample, so it runs in parallel and resumes per sample. Emits a one-row
 * TSV that the workflow joins back onto meta.
 */
process DETECT_PLATFORM {
    tag "${meta.id}"
    label 'process_single'

    container 'docker.io/rocker/r-ver:4.5.2@sha256:fd4ccdd3a4a6f7ef805e2daeee2a0fe3bf126bc231f36351223baecf5a595a4c'

    input:
    tuple val(meta), path(stats), path(reads)

    output:
    tuple val(meta), path("${meta.id}.detection.tsv"), emit: tsv
    tuple val("${task.process}"), val('R'), eval("Rscript --version 2>&1 | sed 's/^.*version //; s/ .*//'"), topic: versions

    script:
    """
    detect_platform.R ${meta.id} ${stats} ${meta.id}.detection.tsv ${reads}
    """
}
