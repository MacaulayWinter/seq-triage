/*
 * CHECK_BAM: is a BAM intact from header to end? (bin/check_bam.sh)
 *
 * Runs before the BAM is converted to FASTQ. A damaged BAM becomes a FAIL
 * row in the report instead of stopping the whole batch (CLAUDE.md D45, D55).
 * Same samtools build as SAMTOOLS_COLLATEFASTQ (nf-core), so the version
 * recorded is the version that read the file.
 */
process CHECK_BAM {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/htslib_samtools_xz:63148285647cb8e3'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.id}.check.tsv"), emit: tsv
    tuple val(meta), path("${meta.id}.header.tsv"), emit: header
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), topic: versions

    script:
    """
    check_bam.sh ${meta.id} ${bam}
    """
}
