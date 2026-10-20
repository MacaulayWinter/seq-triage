/*
 * TRIAGE_REPORT: one summary for all samples (bin/triage_report.R)
 *
 * A single task over every sample: it must see the whole batch to check
 * that no sample has gone missing on the way.
 */
process TRIAGE_REPORT {
    label 'process_single'

    container 'docker.io/rocker/r-ver:4.5.2@sha256:fd4ccdd3a4a6f7ef805e2daeee2a0fe3bf126bc231f36351223baecf5a595a4c'

    input:
    path screen       // metadata_screen.tsv: every row of the samplesheet
    path samples      // samples.tsv: per-sample state after detection and routing
    path metric_guide
    path qc_values    // qc_values.tsv: every QC number, one row per sample/tool/field

    output:
    path 'triage_summary.tsv', emit: tsv
    path 'triage_summary.md',  emit: md
    path 'triage_metrics.tsv', emit: metrics
    tuple val("${task.process}"), val('R'), eval("Rscript --version 2>&1 | sed 's/^.*version //; s/ .*//'"), topic: versions

    script:
    """
    triage_report.R ${screen} ${samples} ${metric_guide} ${qc_values} triage
    """
}
