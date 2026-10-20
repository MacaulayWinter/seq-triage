/*
 * SCREEN_METADATA: warn about samplesheet inconsistencies (bin/screen_metadata.R)
 *
 * One task for the whole sheet, not one per sample: collisions and shared
 * files only show up when rows are compared (CLAUDE.md D5).
 */
process SCREEN_METADATA {
    label 'process_single'

    // Tag + digest: the tag says which R, the digest guarantees the exact
    // image. Multi-arch index (amd64 + arm64). The registry is written out
    // because nextflow.config sets docker.registry = 'quay.io' for nf-core.
    container 'docker.io/rocker/r-ver:4.5.2@sha256:fd4ccdd3a4a6f7ef805e2daeee2a0fe3bf126bc231f36351223baecf5a595a4c'

    input:
    path samplesheet

    output:
    path 'metadata_screen.tsv', emit: tsv
    tuple val("${task.process}"), val('R'), eval("Rscript --version 2>&1 | sed 's/^.*version //; s/ .*//'"), topic: versions

    script:
    """
    screen_metadata.R ${samplesheet} metadata_screen.tsv
    """
}
