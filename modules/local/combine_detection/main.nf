/*
 * COMBINE_DETECTION: one detection row for a sample in several files
 * (bin/combine_detection.R, CLAUDE.md D75)
 *
 * Each lane (a mate pair, or one file of single-end reads) was detected on
 * its own by DETECT_PLATFORM; this joins the lane rows, in lane order, into
 * one row in the same format, so the rest of the workflow cannot tell a
 * several-file sample from a one-file one.
 */
process COMBINE_DETECTION {
    tag "${id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'docker.io/rocker/r-ver:4.5.2@sha256:fd4ccdd3a4a6f7ef805e2daeee2a0fe3bf126bc231f36351223baecf5a595a4c'

    input:
    tuple val(id), path(lane_tsvs)

    output:
    tuple val(id), path("${id}.detection.tsv"), emit: tsv
    tuple val("${task.process}"), val('R'), eval("Rscript --version 2>&1 | sed 's/^.*version //; s/ .*//'"), topic: versions

    script:
    """
    combine_detection.R ${id} ${id}.detection.tsv ${lane_tsvs}
    """
}
