/*
 * FETCH_PATHOGEN_LIST: the list of bacterial human pathogens, resolved to
 * NCBI taxids (bin/pathogen_list.R, D96)
 *
 * Runs only when Zenodo has a version not yet in the cache. Nextflow
 * downloads the Zenodo archive and NCBI's taxonomy dump; the result is stored
 * in --pathogen_cache/<version>/ and reused by every later run. A failure
 * never stops the run: the report then says the list was not checked.
 */
process FETCH_PATHOGEN_LIST {
    label 'process_medium'          // reading NCBI's names.dmp takes about 4 GB
    storeDir "${params.pathogen_cache}/${version}"
    errorStrategy 'ignore'

    conda "${moduleDir}/environment.yml"
    container 'docker.io/rocker/r-ver:4.5.2@sha256:fd4ccdd3a4a6f7ef805e2daeee2a0fe3bf126bc231f36351223baecf5a595a4c'

    input:
    path record, stageAs: 'zenodo_record.zip'   // the Zenodo archive holding Data S1
    path taxdump                                // NCBI taxdump.tar.gz
    val version
    val doi

    output:
    path 'pathogen_list', emit: dir
    tuple val("${task.process}"), val('R'), eval("Rscript --version 2>&1 | sed 's/^.*version //; s/ .*//'"), topic: versions

    script:
    """
    mkdir pathogen_list
    pathogen_list.R zenodo_record.zip ${taxdump} ${version} ${doi} pathogen_list/pathogen_list
    """
}
