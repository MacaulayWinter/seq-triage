/*
 * DEINTERLEAVE: one interleaved FASTQ -> R1 and R2 (bin/deinterleave_fastq.sh, D78)
 *
 * Only for samples whose detection found both mates in one file
 * (detected_layout 'interleaved'), just before QC, so fastp and FastQC see a
 * pair. Same pinned samtools image as NORMALISE_FASTQ (bash, awk, bgzip);
 * nothing is published.
 */
process DEINTERLEAVE {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/htslib_samtools_xz:63148285647cb8e3'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}_{1,2}.fastq.gz"), emit: reads
    tuple val("${task.process}"), val('bgzip'), eval("bgzip --version | sed '1!d;s/.* //'"), topic: versions

    script:
    """
    deinterleave_fastq.sh ${meta.id} ${task.cpus} ${reads}
    """
}
