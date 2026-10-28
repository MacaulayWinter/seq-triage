/*
 * STUDY_METADATA: the study's own variables, checked and described
 * (bin/study_metadata.R, CLAUDE.md D95)
 *
 * One task for the whole sheet. Study metadata is never used for QC: it is
 * carried to the analysis pipelines, checked against the optional data
 * column guide, and the design as read is written in plain language.
 */
process STUDY_METADATA {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'docker.io/rocker/r-ver:4.5.2@sha256:fd4ccdd3a4a6f7ef805e2daeee2a0fe3bf126bc231f36351223baecf5a595a4c'

    input:
    path samplesheet
    path schema
    path metadata      // [] when not given
    path guide        // [] when not given

    output:
    path 'study_metadata.tsv',  emit: tsv
    path 'metadata_check.tsv',  emit: check
    path 'metadata_check.md',   emit: md
    path 'column_guide_draft.csv', optional: true, emit: draft
    tuple val("${task.process}"), val('R'), eval("Rscript --version 2>&1 | sed 's/^.*version //; s/ .*//'"), topic: versions

    script:
    def m = metadata   ? "--metadata ${metadata}"     : ''
    def gd = guide ? "--column_guide ${guide}" : ''
    """
    study_metadata.R ${samplesheet} ${schema} ${m} ${gd}
    """
}
