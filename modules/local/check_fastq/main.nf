/*
 * CHECK_FASTQ: is every file readable and well formed? (bin/check_fastq.sh)
 *
 * Runs before anything else reads the files. A damaged file becomes a FAIL
 * row in the report instead of stopping the whole batch (CLAUDE.md D45).
 */
process CHECK_FASTQ {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'docker.io/rocker/r-ver:4.5.2@sha256:fd4ccdd3a4a6f7ef805e2daeee2a0fe3bf126bc231f36351223baecf5a595a4c'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.check.tsv"), emit: tsv
    tuple val("${task.process}"), val('gzip'), eval("gzip --version 2>&1 | head -1 | sed 's/^gzip //'"), topic: versions

    script:
    """
    check_fastq.sh ${meta.id} ${reads}
    """
}
