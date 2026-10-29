/*
 * FETCH_PATHOGEN_LIST: the list of bacterial human pathogens, resolved to
 * NCBI taxids (bin/pathogen_list.R, D96)
 *
 * Runs only when the list has a version not yet in the cache. The task itself
 * downloads the list and NCBI's taxonomy dump: a download Nextflow stages
 * fails before the task exists, which errorStrategy cannot catch, and would
 * stop the whole run (found in testing, D96). The result is stored in
 * --pathogen_cache/<version>/ and reused by every later run. A failure never
 * stops the run: the report then says the list was not checked.
 */
process FETCH_PATHOGEN_LIST {
    label 'process_medium'          // reading NCBI's names.dmp takes about 4 GB
    storeDir "${params.pathogen_cache}/${version}"
    errorStrategy 'ignore'

    conda "${moduleDir}/environment.yml"
    container 'docker.io/rocker/r-ver:4.5.2@sha256:fd4ccdd3a4a6f7ef805e2daeee2a0fe3bf126bc231f36351223baecf5a595a4c'

    input:
    val record_url      // the list (an .xlsx, or an archive holding it)
    val taxdump_url     // NCBI taxdump.tar.gz
    val version
    val doi

    // No 'versions' output: storeDir allows only path and val outputs. It is
    // the same R image as TRIAGE_REPORT, which records the R version.
    output:
    path 'pathogen_list', emit: dir

    script:
    """
    Rscript -e 'options(timeout = 900)
                download.file("${record_url}", "${record_url.endsWith('.xlsx') ? 'list.xlsx' : 'list.zip'}", mode = "wb", quiet = TRUE)
                download.file("${taxdump_url}", "taxdump.tar.gz", mode = "wb", quiet = TRUE)'
    mkdir pathogen_list
    pathogen_list.R ${record_url.endsWith('.xlsx') ? 'list.xlsx' : 'list.zip'} taxdump.tar.gz ${version} ${doi} pathogen_list/pathogen_list
    """
}
