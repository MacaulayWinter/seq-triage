/*
 * RECOMMEND: the next nf-core pipeline for each ready sample (bin/recommend.R)
 *
 * One task after the report: reads the triage summary and metrics, the
 * user's samplesheet (for the original file paths) and the pinned pipeline
 * table, and writes recommendations/ (CLAUDE.md D66).
 */
process RECOMMEND {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'docker.io/rocker/r-ver:4.5.2@sha256:fd4ccdd3a4a6f7ef805e2daeee2a0fe3bf126bc231f36351223baecf5a595a4c'

    input:
    path summary      // triage_summary.tsv
    path metrics      // triage_metrics.tsv
    path samplesheet  // the user's samplesheet, as given
    path pipelines    // assets/downstream_pipelines.tsv or the user's own

    output:
    path 'recommendations', emit: dir
    tuple val("${task.process}"), val('R'), eval("Rscript --version 2>&1 | sed 's/^.*version //; s/ .*//'"), topic: versions

    script:
    """
    recommend.R ${summary} ${metrics} ${samplesheet} ${pipelines} recommendations
    """
}
