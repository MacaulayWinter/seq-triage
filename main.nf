#!/usr/bin/env nextflow
/*
 * seq-triage: entry point
 *
 * main.nf checks the inputs and builds the sample channel. The analysis
 * wiring will live in workflows/seq_triage.nf (nf-core convention: the entry
 * point handles input, the workflow file handles logic).
 */

include { validateParameters; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'
include { SEQ_TRIAGE } from './workflows/seq_triage'
include { FETCH_PATHOGEN_LIST } from './modules/local/fetch_pathogen_list/main'

// Empty cells arrive as [] (D22); case and spaces do not make two
// declarations differ.
def declaredValue(v) {
    return (v == null || v == []) ? null : v.toString().trim().toLowerCase()
}

// D95: study metadata columns are accepted (nf-schema's own check of unknown
// headers is off), but a header that looks like a pipeline column stops the
// run, so a typo can never be silently treated as something else (D19: a
// missed fastq_2 would turn paired reads single-end). "Looks like": the
// same letters ignoring case and punctuation, one keystroke away, or a
// common synonym (R1, read1, sample_id ...). A column named in the column
// guide (--column_guide) is confirmed as study metadata and is
// not stopped: listing it is how the user says the name is intended.
def levenshtein(String a, String b) {
    def d = (0..b.size()).collect { j -> j }
    (1..a.size()).each { i ->
        def prev = d[0]
        d[0] = i
        (1..b.size()).each { j ->
            def tmp = d[j]
            d[j] = [ d[j] + 1, d[j - 1] + 1, prev + (a[i - 1] == b[j - 1] ? 0 : 1) ].min()
            prev = tmp
        }
    }
    return d[b.size()]
}
def normCol(String x) {
    return x.toLowerCase().replaceAll(/[^a-z0-9]/, '')
}
def headerOf(path) {
    def line = path.withReader { r -> r.readLine() } ?: ''
    def sep = path.name.toLowerCase().endsWith('.tsv') ? '\t' : ','
    return line.replaceFirst(/^\uFEFF/, '').split(sep, -1).collect { c -> c.trim().replaceAll(/^"|"$/, '') }.findAll { c -> c }
}
def checkHeaders(path, List known, List confirmed) {
    def synonyms = [ sampleid: 'sample', samplename: 'sample', name: 'sample', id: 'sample',
                     r1: 'fastq_1', read1: 'fastq_1', reads1: 'fastq_1', fq1: 'fastq_1', forward: 'fastq_1',
                     r2: 'fastq_2', read2: 'fastq_2', reads2: 'fastq_2', fq2: 'fastq_2', reverse: 'fastq_2',
                     cram: 'bam', sam: 'bam', fastqdirectory: 'fastq_dir', folder: 'fastq_dir' ]
    def study = headerOf(path).findAll { c -> !(c in known) }
    def problems = []
    study.findAll { c -> !(c in confirmed) }.each { c ->
        def n = normCol(c)
        // The closest pipeline column: the same letters first, then a synonym, then one keystroke away.
        def hit = known.find { k -> normCol(k) == n } ?: synonyms[n] ?: known.find { k -> levenshtein(normCol(k), n) <= 1 }
        if (hit) {
            problems << "'${c}' looks like the pipeline column '${hit}'"
        }
    }
    if (problems) {
        error("Samplesheet header: " + problems.join('; ') + ". Nothing has been assumed: if it is that column, rename it exactly " +
              "(as in the README); if it is your own study variable, list it in --column_guide to confirm the name is intended.")
    }
    return study
}

// FASTQ files in a fastq_dir folder (D75), sorted by name.
def dirFastqs(dir) {
    def all   = dir.listFiles()
    def found = all.findAll { f -> f.name ==~ /.+\.(f(ast)?q|fa|fasta|fna)(\.(gz|bz2|xz|zst))?$/ }.sort { f -> f.name }
    if (!found) {
        // Say what the folder IS when it is not reads yet (D79).
        def names = all.collect { f -> f.name }
        def what  = names.any { n -> n ==~ /.+\.(pod5|fast5)$/ } ? ' It holds raw nanopore signal (POD5/FAST5): basecall it first (e.g. dorado), then give the FASTQ or BAM.'
                  : names.contains('RunInfo.xml') || names.contains('Data') ? ' It looks like an Illumina run folder (BCL): demultiplex it first (bcl-convert or bcl2fastq), then give the FASTQ files.'
                  : names.any { n -> n.endsWith('.sra') } ? ' It holds SRA files: convert them first (fasterq-dump).'
                  : ''
        error("fastq_dir ${dir} holds no FASTQ or FASTA read files (.fastq/.fq, .fa/.fasta/.fna, plain or compressed).${what}")
    }
    return found
}

// One sample from its samplesheet rows: [ meta, [ files ] ]. Files are in
// lane order, mates together (L1 R1, L1 R2, L2 R1, ...), as nf-core's
// CAT_FASTQ expects. meta gains n_units only when there is more than one
// lane or file, so a one-file sample's meta (and its cache) is unchanged.
def sampleFromRows(List rs) {
    def meta  = rs[0][0].collectEntries { k, v -> [ k, v == [] ? null : v ] }
    def first = rs[0]
    if (first[3]) {
        return [ meta + [ single_end: false, input_format: 'bam' ], [ first[3] ] ]
    }
    def single_end = !first[2]
    def units = rs.collectMany { _m, r1, r2, _bam, dir, _target, _gsize, _group -> dir ? dirFastqs(dir).collect { f -> [ f ] } : [ r2 ? [ r1, r2 ] : [ r1 ] ] }
    def extra = units.size() > 1 ? [ n_units: units.size() ] : [:]
    return [ meta + [ single_end: single_end ] + extra, units.flatten() ]
}

// D96: Bartlett et al. 2022's list of bacterial human pathogens, from the
// authors' GitHub repository: the newest commit that changed the list file is
// looked up at each run (the user's choice, 2026-10-29: the authors corrected
// names there in 2025 without a new Zenodo version). Cited with the paper and
// the Zenodo DOI. NCBI's taxonomy turns names into taxids once per version.
def PATHOGEN_REPO() { 'padpadpadpad/bartlett_et_al_2022_human_pathogens' }
def PATHOGEN_FILE() { 'data/bacteria_human_pathogens.xlsx' }
def PATHOGEN_DOI() { '10.5281/zenodo.7337534' }
def NCBI_TAXDUMP_URL() { 'https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz' }

// WHAT IT DOES  Decides where the pathogen list comes from (D96): skipped;
//               a prepared folder (--pathogen_list); the cached copy of the
//               newest commit; a fresh build of it; or, offline, the newest
//               cached copy, labelled as such. Never stops the run.
// RETURNS       [mode: 'none' | 'dir' | 'fetch', note: <for the report>,
//               dir: <folder> | version, doi, zip: <download URL>]
def pathogenListSource() {
    if (params.skip_pathogen_flags) {
        return [mode: 'none', note: 'turned off with --skip_pathogen_flags']
    }
    if (params.pathogen_list) {
        def d = file(params.pathogen_list, type: 'dir', checkIfExists: true)
        if (!d.resolve('pathogen_list.tsv').exists() || !d.resolve('pathogen_list_families.tsv').exists()) {
            error("--pathogen_list must be a folder written by bin/pathogen_list.R (pathogen_list.tsv, pathogen_list_families.tsv, pathogen_list_unresolved.tsv): ${d}")
        }
        return [mode: 'dir', dir: d, note: 'the copy given with --pathogen_list (latest not checked)']
    }
    def cache = file(params.pathogen_cache, type: 'dir')
    def today = new Date().format('yyyy-MM-dd')
    try {
        def conn = new java.net.URI("https://api.github.com/repos/${PATHOGEN_REPO()}/commits?path=${PATHOGEN_FILE()}&per_page=1").toURL().openConnection()
        conn.setConnectTimeout(20000)
        conn.setReadTimeout(20000)
        def commit = new groovy.json.JsonSlurper().parse(conn.getInputStream())[0]
        def sha = commit.sha.toString()
        // e.g. github-2380c3a265-2025-07-28: the commit and the day it was made
        def version = "github-${sha.take(10)}-${commit.commit.committer.date.toString().take(10)}"
        def hit = cache.resolve(version).resolve('pathogen_list')
        if (hit.resolve('pathogen_list.tsv').exists()) {
            return [mode: 'dir', dir: hit, note: "the newest on the authors' GitHub on ${today}, from the cache"]
        }
        return [mode: 'fetch', version: version, doi: PATHOGEN_DOI(),
                zip: "https://raw.githubusercontent.com/${PATHOGEN_REPO()}/${sha}/${PATHOGEN_FILE()}",
                note: "the newest on the authors' GitHub on ${today}, downloaded for this run"]
    } catch (Exception e) {
        def cached = cache.exists() ? cache.listFiles().findAll { v -> v.resolve('pathogen_list/pathogen_list.tsv').exists() }.sort { v -> v.lastModified() } : []
        if (cached) {
            return [mode: 'dir', dir: cached[-1].resolve('pathogen_list'),
                    note: "from the cache; latest not checked, as GitHub could not be reached (${e.getClass().getSimpleName()})"]
        }
        return [mode: 'none', note: "the pathogen list is unavailable: GitHub could not be reached and nothing is cached (${e.getClass().getSimpleName()})"]
    }
}

workflow {
    // --help / --helpFull: nf-schema has already printed the help (generated
    // from nextflow_schema.json, D53). Stop here, so a help request is not
    // also reported as a run with no --input.
    if (params.containsKey('help') || params.containsKey('helpFull')) {
        exit(0)
    }

    // Stop before any work if a parameter is missing, mistyped or points at
    // a file that does not exist. Rules come from nextflow_schema.json, which
    // also validates the samplesheet against assets/schema_input.json.
    validateParameters()

    // Print the parameters that differ from their defaults, so every run log
    // records exactly how the pipeline was invoked.
    log.info paramsSummaryLog(workflow)

    // One [meta, fastq_1, fastq_2, bam] list per row. Schema fields marked
    // "meta" become keys of the meta map: id, declared_platform, declared_type.
    // A row gives FASTQ files or one BAM, never both (the schema checks it).
    // D95: any other column is study metadata, carried to metadata/study_metadata.tsv.
    def known_cols = new groovy.json.JsonSlurper().parse(file("${projectDir}/assets/schema_input.json")).items.properties.keySet() as List
    def confirmed = params.column_guide ? file(params.column_guide, checkIfExists: true).readLines().drop(1).collect { l ->
        def c = l.split(params.column_guide.toString().toLowerCase().endsWith('.tsv') ? '\t' : ',', -1)[0]
        c.trim().replaceAll(/^"|"$/, '') } : []
    def study_cols = checkHeaders(file(params.input), known_cols, confirmed)
    if (study_cols) {
        log.info "Study metadata columns (carried through, not used for QC): ${study_cols.join(', ')}"
    }
    def rows = samplesheetToList(params.input, "${projectDir}/assets/schema_input.json")
    def same_mates = rows.findAll { _meta, r1, r2, _bam, _dir, _target, _gsize, _group -> r2 && r1 == r2 }
    if (same_mates) {
        error("fastq_1 and fastq_2 are the same file for sample(s): " +
              same_mates.collect { meta, _r1, _r2, _bam, _dir, _target, _gsize, _group -> meta.id }.join(', '))
    }

    // Several rows with the same sample ID are one sample whose reads are in
    // several files: lanes, top-up runs (D75). They must agree on everything
    // but the files, or the run stops here, naming the sample.
    def problems = []
    rows.groupBy { row -> row[0].id }.each { id, rs ->
        def files = rs.collectMany { _meta, r1, r2, bam, dir, _target, _gsize, _group -> [ r1, r2, bam, dir ].findAll { f -> f } }
        def twice = files.findAll { f -> files.count(f) > 1 }.unique()
        if (twice) {
            problems << "${id}: the same file is listed twice (${twice.collect { f -> f.name }.join(', ')})"
        }
        if (rs.size() == 1) {
            return
        }
        if (rs.any { row -> row[3] }) {
            problems << "${id}: a BAM, SAM or CRAM file must be the only row of its sample"
        } else if (rs.collect { row -> row[2] ? 'paired' : 'single' }.unique().size() > 1) {
            problems << "${id}: some rows are paired (fastq_2 given) and some are not"
        }
        [ declared_platform: 'platform', declared_type: 'sample_type', declared_molecule: 'molecule', declared_library: 'library' ].each { key, column ->
            def values = rs.collect { row -> declaredValue(row[0][key]) }.unique()
            if (values.size() > 1) {
                problems << "${id}: rows disagree on ${column} (${values.collect { v -> v ?: 'blank' }.join(', ')})"
            }
        }
        def targets = rs.collect { row -> declaredValue(row[5]) }.unique()
        if (targets.size() > 1) {
            problems << "${id}: rows disagree on target_taxid (${targets.collect { v -> v ?: 'blank' }.join(', ')})"
        }
        def sizes = rs.collect { row -> declaredValue(row[6]) }.unique()
        if (sizes.size() > 1) {
            problems << "${id}: rows disagree on genome_size (${sizes.collect { v -> v ?: 'blank' }.join(', ')})"
        }
        def groups = rs.collect { row -> declaredValue(row[7]) }.unique()
        if (groups.size() > 1) {
            problems << "${id}: rows disagree on group (${groups.collect { v -> v ?: 'blank' }.join(', ')})"
        }
    }
    if (problems) {
        error("Samplesheet rows of the same sample must agree (D75):\n  " + problems.join('\n  '))
    }

    ch_samples = channel.fromList(rows.groupBy { row -> row[0].id }.values() as List)
        .map { rs -> sampleFromRows(rs) }

    // D96: the known-pathogen list. Looked up here, on the machine running
    //      Nextflow (an HPC login node has the internet; compute nodes may
    //      not); built by FETCH_PATHOGEN_LIST only for a version not cached.
    def pl = pathogenListSource()
    def no_list = file("${projectDir}/assets/no_pathogen_list", type: 'dir')
    def ch_pathogen_dir
    def ch_pathogen_note
    if (pl.mode == 'fetch') {
        log.info "Pathogen list: building version ${pl.version} (first use; cached in ${params.pathogen_cache})"
        FETCH_PATHOGEN_LIST(pl.zip, NCBI_TAXDUMP_URL(), pl.version, pl.doi)
        ch_pathogen_dir  = FETCH_PATHOGEN_LIST.out.dir.ifEmpty(no_list)
        ch_pathogen_note = FETCH_PATHOGEN_LIST.out.dir.map { _d -> pl.note }
            .ifEmpty("the pathogen list could not be built (download or processing failed, see the run log)")
    } else {
        ch_pathogen_dir  = channel.value(pl.mode == 'dir' ? pl.dir : no_list)
        ch_pathogen_note = channel.value(pl.note)
    }

    SEQ_TRIAGE(ch_samples, file(params.input), file(params.metric_guide), file(params.species_groups), file(params.genome_sizes),
               ch_pathogen_dir, ch_pathogen_note)
}
