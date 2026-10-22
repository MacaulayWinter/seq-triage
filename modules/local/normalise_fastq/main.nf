/*
 * NORMALISE_FASTQ: every read file becomes gzipped FASTQ (bin/normalise_fastq.sh)
 *
 * Runs only for samples with a file that is not already .gz: plain FASTQ,
 * .bz2, .xz or .zst (CLAUDE.md D69). A file that cannot be decompressed
 * becomes a FAIL row (D45). Same pinned samtools image as CHECK_BAM, which
 * carries gzip, bzip2, xz, zstd and bgzip; nothing is published.
 */
process NORMALISE_FASTQ {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/htslib_samtools_xz:63148285647cb8e3'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.normalise.tsv"), path("${meta.id}_*.fastq.gz", arity: '0..*'), emit: reads
    tuple val("${task.process}"), val('bgzip'), eval("bgzip --version | sed '1!d;s/.* //'"), topic: versions

    script:
    """
    normalise_fastq.sh ${meta.id} ${task.cpus} ${reads}
    """
}
