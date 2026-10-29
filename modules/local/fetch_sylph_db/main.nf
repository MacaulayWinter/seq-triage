/*
 * FETCH_SYLPH_DB: the sylph database (D99), downloaded once
 *
 * Runs only when the database named by --sylph_db is a URL not yet in the
 * cache. The task downloads it itself (a download Nextflow stages would fail
 * before the task exists and stop the whole run, D96), into
 * --sylph_cache/, where every later run finds it. A failure never stops the
 * run: off-target is then estimated by Kraken2/Bracken, as before D99.
 */
process FETCH_SYLPH_DB {
    label 'process_single'
    storeDir "${params.sylph_cache}"
    errorStrategy 'ignore'

    conda "${moduleDir}/environment.yml"
    container 'docker.io/rocker/r-ver:4.5.2@sha256:fd4ccdd3a4a6f7ef805e2daeee2a0fe3bf126bc231f36351223baecf5a595a4c'

    input:
    val url

    output:
    path "${file(url).name}", emit: db

    script:
    def name = file(url).name
    """
    Rscript -e 'options(timeout = 7200)
                download.file("${url}", "${name}.part", mode = "wb", quiet = TRUE)'
    mv ${name}.part ${name}
    """
}
