/*
 * KRAKEN2_BATCH: every sample in one task, one database load (D50)
 *
 * Used instead of the per-sample nf-core KRAKEN2 when --kraken2_batch is
 * set: one read of the database for the whole run, at the cost of running
 * samples one after another. Best for big databases on an HPC.
 *
 * Several samples may share a file (the test profile does), so each file
 * is staged into its own numbered folder (in1/, in2/, ...) and the
 * manifest says which files belong to which sample.
 */
process KRAKEN2_BATCH {
    tag "${metas.size()} samples"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0f/0f827dcea51be6b5c32255167caa2dfb65607caecdc8b067abd6b71c267e2e82/data' :
        'community.wave.seqera.io/library/kraken2_coreutils_pigz:920ecc6b96e2ba71' }"

    input:
    tuple val(metas), val(manifest), path(reads, stageAs: 'in*/*')
    path db

    output:
    tuple val(metas), path('*.kraken2.report.txt'), emit: reports
    tuple val("${task.process}"), val('kraken2'), eval('kraken2 --version 2>&1 | head -1 | sed "s/^.*Kraken version //; s/ .*//"'), topic: versions

    script:
    """
    cat > manifest.tsv <<'MANIFEST'
${manifest}
MANIFEST
    kraken2_batch.sh ${db} ${task.cpus} manifest.tsv
    """
}
