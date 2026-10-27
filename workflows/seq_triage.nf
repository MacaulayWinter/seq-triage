/*
 * SEQ_TRIAGE: screen -> detect -> join -> route -> QC -> report
 *
 * The channel diagram in CLAUDE.md ("Architecture") is the map for this
 * file. Our own logic is the joins, routeSample() and the branch; every
 * standard tool is an nf-core module.
 */

include { SCREEN_METADATA } from '../modules/local/screen_metadata/main'
include { CHECK_FASTQ     } from '../modules/local/check_fastq/main'
include { DETECT_PLATFORM } from '../modules/local/detect_platform/main'
include { TRIAGE_REPORT   } from '../modules/local/triage_report/main'
include { SEQKIT_STATS    } from '../modules/nf-core/seqkit/stats/main'
include { FASTP           } from '../modules/nf-core/fastp/main'
include { FASTQC          } from '../modules/nf-core/fastqc/main'
include { NANOPLOT        } from '../modules/nf-core/nanoplot/main'
include { MULTIQC         } from '../modules/nf-core/multiqc/main'
include { UNTAR           } from '../modules/nf-core/untar/main'
include { KRAKEN2_KRAKEN2 as KRAKEN2 } from '../modules/nf-core/kraken2/kraken2/main'
include { BRACKEN_BRACKEN as BRACKEN } from '../modules/nf-core/bracken/bracken/main'
include { KRAKEN2_BATCH   } from '../modules/local/kraken2_batch/main'
include { RECOMMEND       } from '../modules/local/recommend/main'
include { CHECK_BAM       } from '../modules/local/check_bam/main'
include { NORMALISE_FASTQ } from '../modules/local/normalise_fastq/main'
include { SAMTOOLS_COLLATEFASTQ } from '../modules/nf-core/samtools/collatefastq/main'
include { CAT_FASTQ       } from '../modules/nf-core/cat/fastq/main'
include { COMBINE_DETECTION } from '../modules/local/combine_detection/main'
include { DEINTERLEAVE    } from '../modules/local/deinterleave/main'
include { SORTMERNA       } from '../modules/nf-core/sortmerna/main'
include { SORTMERNA as SORTMERNA_INDEX } from '../modules/nf-core/sortmerna/main'
include { GUNZIP as GUNZIP_RRNA        } from '../modules/nf-core/gunzip/main'

// Empty TSV cells arrive as "", so convert them to null: one meaning of
// "not present" throughout the pipeline (CLAUDE.md D22).
// Does a gzipped FASTQ hold at least one read? Reads only the first byte.
def hasReads(f) {
    new java.util.zip.GZIPInputStream(f.newInputStream()).withCloseable { s -> s.read() != -1 }
}

// Reads converted from a BAM (D55): the layout comes from the reads' own
// flags. READ1 and READ2 both present -> paired; otherwise the reads with
// neither flag (unaligned long-read BAMs from dorado or PacBio, single-end
// Illumina) -> single. Reads left out are said in the note, never dropped
// silently.
def bamReads(Map meta, pair, other, singleton) {
    def (r1, r2) = [ pair ].flatten().sort { f -> f.name }
    def note
    def out
    if (hasReads(r1) && hasReads(r2)) {
        out  = [ meta + [ single_end: false ], [ r1, r2 ] ]
        note = 'Reads extracted from the BAM/SAM/CRAM file: paired by their READ1/READ2 flags.'
        if (hasReads(other) || hasReads(singleton)) {
            note += ' Some reads had no mate or no READ1/READ2 flag and are not used.'
        }
    } else if (hasReads(other)) {
        out  = [ meta + [ single_end: true ], [ other ] ]
        note = 'Reads extracted from the BAM/SAM/CRAM file: no READ1/READ2 flags, so treated as one file of reads.'
        if (hasReads(r1) || hasReads(r2)) {
            note += ' Reads flagged as only one mate are not used.'
        }
    } else {
        def first = hasReads(r1) ? r1 : (hasReads(r2) ? r2 : other)
        out  = [ meta + [ single_end: true ], [ first ] ]
        note = hasReads(first) ? 'Reads extracted from the BAM/SAM/CRAM file: only one mate present, treated as single.'
                               : 'Reads extracted from the BAM/SAM/CRAM file: the BAM holds no primary reads.'
    }
    out[0] = out[0] + [ bam_note: note ]
    return out
}

// The BAM header's @RG PL as a platform (D55). PL is typed by whoever made
// the file, so it is treated like a declared platform, never as evidence.
// Only values that name a platform seq-triage routes are used; several
// different ones mean nothing can be taken from it.
def rgPlatform(String pl) {
    if (!pl) {
        return null
    }
    def known = [ illumina: 'illumina', ont: 'ont', nanopore: 'ont', oxford_nanopore: 'ont', pacbio: 'pacbio' ]
    def hits = pl.split(',').collect { p -> known[p.trim().toLowerCase()] }.unique()
    return hits.size() == 1 ? hits[0] : null
}

// The lanes of a sample (D75): [ meta, files ] per lane (a mate pair, or
// one file of single-end reads), each with id "<sample>.u<k>", unit_of and
// unit. A one-file sample comes back unchanged.
def lanes(Map meta, reads) {
    def files = [ reads ].flatten()
    if (!meta.n_units) {
        return [ [ meta, reads ] ]
    }
    def groups = meta.single_end ? files.collect { f -> [ f ] } : files.collate(2)
    return groups.withIndex().collect { g, k ->
        [ meta + [ id: "${meta.id}.u${k + 1}".toString(), unit_of: meta.id, unit: k + 1 ], g ]
    }
}

def blankToNull(String s) {
    return (s == null || s.trim().isEmpty()) ? null : s
}

/*
 * QC values -> one long table (CLAUDE.md D34)
 *
 * Each helper turns one tool's output file into rows of
 *     sample <TAB> tool <TAB> file <TAB> field <TAB> value
 * where `field` uses the same names as metric_guide.tsv's source_field.
 * This is format conversion only: no thresholds, no judgement. All of that
 * stays in bin/triage_report.R. Done here because Groovy reads JSON and zip
 * files natively, and the pinned R container cannot read JSON.
 */
def qcRow(String sample, String tool, String file, String field, Object value) {
    return [ sample, tool, file, field, value ].join('\t')
}

// Nested JSON -> [ "summary.after_filtering.q30_rate", 0.88 ] pairs. Only
// numbers are kept; lists (per-cycle curves) and text are not needed.
def flattenJson(Object node, String prefix) {
    if (node instanceof Map) {
        return node.collectMany { k, v -> flattenJson(v, prefix ? "${prefix}.${k}" : k.toString()) }
    }
    if (node instanceof Number) {
        return [ [ prefix, node ] ]
    }
    return []
}

def fastpRows(String sample, json) {
    def parsed = new groovy.json.JsonSlurper().parseText(json.text)
    return flattenJson(parsed, '').collect { kv -> qcRow(sample, 'fastp', json.name, kv[0], kv[1]) }
}

// NanoStats.txt written with --tsv_stats: a header line, then "field<TAB>value".
def nanoplotRows(String sample, stats) {
    return stats.readLines().drop(1)
        .collect { line -> line.split('\t') }
        .findAll { parts -> parts.size() == 2 }
        .collect { parts -> qcRow(sample, 'NanoPlot', stats.name, parts[0].trim(), parts[1].trim()) }
}

// FastQC's zip holds summary.txt: "STATUS<TAB>Module name<TAB>file" per module.
// withCloseable closes the zip afterwards (Nextflow's strict syntax has no
// try/finally).
def fastqcRows(String sample, zip) {
    return new java.util.zip.ZipFile(zip.toFile()).withCloseable { zf ->
        def entry = zf.entries().find { e -> e.name.endsWith('/summary.txt') }
        if (!entry) {
            return []
        }
        zf.getInputStream(entry).text.readLines()
            .collect { line -> line.split('\t') }
            .findAll { parts -> parts.size() >= 2 }
            .collect { parts -> qcRow(sample, 'FastQC', zip.name, parts[1], parts[0]) }
    }
}

// Kraken2 report: "pct<TAB>clade reads<TAB>direct reads<TAB>rank<TAB>taxid<TAB>name"
// (name indented by depth). Field = "rank|taxid|name", value = clade reads.
def kraken2Rows(String sample, report) {
    return report.readLines()
        .collect { line -> line.split('\t') }
        .findAll { p -> p.size() >= 6 }
        .collect { p -> qcRow(sample, 'Kraken2', report.name, "${p[3].trim()}|${p[4].trim()}|${p[5].trim()}", p[1].trim()) }
}

// Bracken species table: name, taxonomy_id, taxonomy_lvl, kraken_assigned_reads,
// added_reads, new_est_reads, fraction_total_reads. Value = re-estimated reads.
def brackenRows(String sample, tsv) {
    return tsv.splitCsv(header: true, sep: '\t').collect { r ->
        qcRow(sample, 'Bracken', tsv.name, "${r.taxonomy_lvl}|${r.taxonomy_id}|${r.name}", r.new_est_reads)
    }
}

// Largest species-level clade count in a Kraken2 report. Bracken stops with
// "no reads found" (and writes nothing, which would stop the run) when no
// species has at least 10 reads, so it only runs above that (D41).
def maxSpeciesReads(report) {
    def counts = report.readLines()
        .collect { line -> line.split('\t') }
        .findAll { p -> p.size() >= 6 && p[3].trim() == 'S' }
        .collect { p -> p[1].trim() as Long }
    return counts ? counts.max() : 0L
}

// SortMeRNA log (D81): "Total reads passing E-value threshold = 1234 (12.34)".
// Field rrna_fraction = that percentage / 100.
def sortmernaRows(String sample, log) {
    def line = log.readLines().find { l -> l.contains('Total reads passing E-value threshold') }
    if (!line) {
        return []
    }
    // One find() only: Groovy's truth test of a Matcher already calls
    // find(), so "m && m.find()" looks for a SECOND match and drops the row
    // (found on real SortMeRNA logs, D81).
    def m = (line =~ /\(\s*([0-9.]+)\s*%?\s*\)/)
    return m.find() ? [ qcRow(sample, 'SortMeRNA', log.name, 'rrna_fraction', (m.group(1) as BigDecimal) / 100) ] : []
}

// seqkit stats: one row per FASTQ file; the report sums over a sample's files.
def seqkitRows(String sample, tsv) {
    return tsv.splitCsv(header: true, sep: '\t').collectMany { row ->
        [ 'num_seqs', 'sum_len', 'sum_n' ].collect { k -> qcRow(sample, 'seqkit', row.file, k, row[k]) }
    }
}

/*
 * routeSample: decide which QC a sample gets, and on what basis (CLAUDE.md D2).
 *   damaged file              -> route 'unreadable' (FAIL; nothing else runs, D45)
 *   no reads                  -> route 'empty'    (FAIL; no QC can run)
 *   high-confidence detection -> route by DETECTED platform
 *   declared PacBio without high-confidence PacBio HiFi detection
 *                             -> 'unrouted' (routed_by 'pacbio_not_hifi'; only HiFi
 *                                is QC'd, D64)
 *   medium/low detection that CONTRADICTS the declared platform
 *                             -> 'unrouted' (routed_by 'contradicted'; UNCERTAIN,
 *                                review): neither source can be trusted enough to
 *                                pick the QC (D60)
 *   otherwise, if declared    -> route by DECLARED platform (verdict capped at WARN)
 *   otherwise                 -> 'unrouted'      (UNCERTAIN; no QC run)
 * Then one safety rule (D33): a sample headed for ONT QC with two files is
 * stopped ('unrouted', routed_by 'blocked'). Nanopore reads are never
 * paired, and NanoPlot fails on two files, which would stop the whole run.
 * A sample confidently identified as a platform with no QC here (PacBio,
 * MGI/DNBSEQ, Element) is also 'unrouted' (routed_by 'unsupported', D40).
 * Returns meta plus route ('illumina' | 'ont' | 'pacbio' | 'empty' | 'unreadable' | 'unrouted') and
 * routed_by ('detected' | 'declared' | 'contradicted' | 'pacbio_not_hifi' | 'blocked' |
 * 'unsupported' | 'none'). PacBio is routed only as HiFi (D64): detection names
 * PacBio without PLATFORM_NOT_HANDLED only when HiFi markers are present.
 */
// Platforms seq-triage QCs, by read class (D69): short reads get fastp +
// FastQC, long reads NanoPlot.
def readClass(String platform) {
    return [ illumina: 'short', dnbseq: 'short', element: 'short', ont: 'long', pacbio: 'long' ][platform]
}

def routeSample(Map meta) {
    def r
    if (meta.detection_flags == 'FILE_UNREADABLE') {
        return meta + [ route: 'unreadable', routed_by: 'none' ]
    }
    // Only index or barcode reads were given (D78): nothing to QC.
    if ((meta.detection_flags ?: '').split(',').contains('INDEX_READS')) {
        return meta + [ route: 'unrouted', routed_by: 'index_reads' ]
    }
    if (meta.n_reads == 0) {
        r = [ route: 'empty', routed_by: 'none' ]
    } else if (meta.confidence == 'high') {
        r = [ route: meta.detected_platform, routed_by: 'detected' ]
    } else if (meta.declared_platform && readClass(meta.detected_platform as String)
               && readClass(meta.detected_platform as String) != readClass(meta.declared_platform as String)) {
        // D60, refined in D69: a contradiction is short reads vs long reads.
        // A medium "illumina" call is a SHORT-READ call (without names the
        // reads cannot tell Illumina from MGI or Element), so a declared
        // dnbseq or element is not contradicted by it.
        r = [ route: 'unrouted', routed_by: 'contradicted' ]
    } else if (meta.declared_platform == 'pacbio') {
        r = [ route: 'unrouted', routed_by: 'pacbio_not_hifi' ]
    } else if (meta.declared_platform) {
        r = [ route: meta.declared_platform, routed_by: 'declared' ]
    } else {
        r = [ route: 'unrouted', routed_by: 'none' ]
    }
    // Several single files of long reads (lanes, a barcode folder; D75)
    // are merged, so only a mate PAIR of long reads is blocked.
    if (readClass(r.route as String) == 'long' && !meta.single_end) {
        r = [ route: 'unrouted', routed_by: 'blocked' ]
    }
    // A named platform seq-triage has no QC for (PacBio, MGI/DNBSEQ,
    // Element; D40): reported by name, no QC.
    def not_handled = (meta.detection_flags ?: '').split(',').contains('PLATFORM_NOT_HANDLED')
    if (r.routed_by == 'detected' && (!readClass(r.route as String) || not_handled)) {
        r = [ route: 'unrouted', routed_by: 'unsupported' ]
    }
    return meta + r
}

workflow SEQ_TRIAGE {

    take:
    ch_samples      // [ meta, [ reads ] ] from the samplesheet
    samplesheet     // the raw samplesheet file, for the metadata screen
    metric_guide    // assets/metric_guide.tsv or the user's own
    species_groups  // assets/species_groups.tsv or the user's own
    genome_sizes    // assets/species_genome_size.tsv (NCBI), for depth of coverage (D80)

    main:

    // 1. Metadata screen: one task for the whole sheet, then one
    //    [ id, row ] element per sample so it can be joined.
    SCREEN_METADATA(samplesheet)
    ch_screen = SCREEN_METADATA.out.tsv
        .splitCsv(header: true, sep: '\t')
        .map { row -> [ row.sample, row ] }

    // 1b. BAM input (D55): check each BAM is intact from header to end
    //     (CHECK_BAM), then convert it to FASTQ with samtools collate |
    //     samtools fastq (nf-core SAMTOOLS_COLLATEFASTQ). Secondary and supplementary alignments are left out
    //     (conf/modules.config), so each read is counted once. A damaged BAM
    //     becomes a FAIL row, like a damaged FASTQ (D45). Everything after
    //     this point sees FASTQ only.
    ch_input = ch_samples.branch { meta, _files ->
        bam:   meta.input_format == 'bam'
        fastq: true
    }
    // CRAM needs the reference it was written against (D69): --cram_fasta, or
    // samtools fetches it from EBI's reference registry by checksum.
    ch_cram_ref = params.cram_fasta ? file(params.cram_fasta, checkIfExists: true) : []
    ch_bam = ch_input.bam.map { meta, files -> [ meta, files[0] ] }
    CHECK_BAM(ch_bam, ch_cram_ref)
    ch_bam_check = ch_bam
        .join(CHECK_BAM.out.tsv, failOnMismatch: true, failOnDuplicate: true)
        .join(CHECK_BAM.out.header, failOnMismatch: true, failOnDuplicate: true)
        .branch { meta, bam, check, header ->
            def row = check.splitCsv(header: true, sep: '\t')[0]
            def hdr = header.splitCsv(header: true, sep: '\t')
            ok: row.status == 'ok'
                return [ meta + [ bam_rg_platform: hdr ? blankToNull(hdr[0].rg_platform) : null,
                                  bam_programs   : hdr ? blankToNull(hdr[0].programs)    : null ], bam ]
            damaged: true
                return [ meta, bam, "${row.file}: ${row.problem}" ]
        }
    // collate first: mates are written as pairs only when they sit next to
    // each other, and a coordinate-sorted BAM scatters them (the first test
    // run kept 60 of 200 reads without it). The reference is used for CRAM
    // only; samtools builds its .fai itself.
    SAMTOOLS_COLLATEFASTQ(ch_bam_check.ok,
        params.cram_fasta ? [ [ id: 'cram_reference' ], file(params.cram_fasta, checkIfExists: true), [] ] : [ [], [], [] ],
        false)
    ch_from_bam = SAMTOOLS_COLLATEFASTQ.out.fastq
        .join(SAMTOOLS_COLLATEFASTQ.out.fastq_other, failOnMismatch: true, failOnDuplicate: true)
        .join(SAMTOOLS_COLLATEFASTQ.out.fastq_singleton, failOnMismatch: true, failOnDuplicate: true)
        .map { meta, pair, other, singleton -> bamReads(meta, pair, other, singleton) }

    // 1c. Other compression (D69): any FASTQ file that is not .gz (plain,
    //     .bz2, .xz, .zst) is rewritten as .fastq.gz first, so every later
    //     step sees one format. Samples already all .gz skip this, so their
    //     cached results stay valid. A file that will not decompress becomes
    //     a FAIL row, like a damaged BAM or gzip (D45).
    ch_fastq = ch_input.fastq.branch { _meta, files ->
        // D79: FASTA reads are converted to FASTQ too, compressed or not.
        gz:    files.every { f -> f.name ==~ /.+\.f(ast)?q\.gz$/ }
        other: true
    }
    NORMALISE_FASTQ(ch_fastq.other)
    ch_normalised = NORMALISE_FASTQ.out.reads.branch { meta, report, files ->
        def rows = report.splitCsv(header: true, sep: '\t')
        ok: rows.every { r -> r.status == 'ok' }
            // Back in input order (NORMALISE_FASTQ numbers them 1, 2, ...):
            // lane order and mate order matter for several-file samples (D75).
            // FASTA reads (D79) are marked, so the report can say why there
            // is no quality to judge.
            def fasta = rows.any { r -> r.file ==~ /.+\.(fa|fasta|fna)(\.(gz|bz2|xz|zst))?$/ }
            return [ meta + (fasta ? [ from_fasta: true ] : [:]), [ files ].flatten().sort { f -> (f.name =~ /_(\d+)\.fastq\.gz$/)[0][1] as Integer } ]
        damaged: true
            return [ meta, [ files ].flatten(), rows.findAll { r -> r.status != 'ok' }.collect { r -> "${r.file}: ${r.problem}" }.join('; ') ]
    }
    // Every input that could not be read before the integrity check, in one
    // place: [ meta, files, problem ].
    ch_damaged_input = ch_bam_check.damaged.map { meta, bam, problem -> [ meta, [ bam ], problem ] }
        .mix(ch_normalised.damaged)

    ch_reads = ch_fastq.gz.mix(ch_normalised.ok).mix(ch_from_bam)

    // 2. Integrity check (D45): is every file readable and well-formed
    //    FASTQ? A damaged sample skips every later step and becomes a FAIL
    //    row; the rest of the batch carries on. Readable samples keep their
    //    meta unchanged, so cached results stay valid.
    CHECK_FASTQ(ch_reads)
    ch_integrity = ch_reads
        .join(CHECK_FASTQ.out.tsv, failOnMismatch: true, failOnDuplicate: true)
        .branch { meta, reads, check ->
            def rows = check.splitCsv(header: true, sep: '\t')
            ok: rows.every { r -> r.status == 'ok' }
                return [ meta, reads ]
            unreadable: true
                return [ meta, rows.findAll { r -> r.status != 'ok' }.collect { r -> "${r.file}: ${r.problem}" }.join('; '),
                         [ reads ].flatten().size() ]
        }
    ch_readable = ch_integrity.ok

    // 3. Read statistics and platform detection, per readable sample.
    //    Detection needs the stats AND the reads (headers, quality values),
    //    so the stats are joined back onto the reads by meta.
    //    A sample in several files (D75) is detected per lane (or per file
    //    for single-end reads), then the lanes are combined: a top-up run at
    //    another read length would blur the read-length spread if the files
    //    were pooled first. Each lane is its own task ("<sample>.u<k>");
    //    one-file samples are untouched.
    ch_units = ch_readable.flatMap { meta, reads -> lanes(meta, reads) }
    SEQKIT_STATS(ch_units)
    DETECT_PLATFORM(
        SEQKIT_STATS.out.stats.join(ch_units, failOnMismatch: true, failOnDuplicate: true)
    )
    ch_unit_tsv = DETECT_PLATFORM.out.tsv.branch { meta, _tsv ->
        lane:  meta.unit_of
        whole: true
    }
    COMBINE_DETECTION(
        ch_unit_tsv.lane
            .map { meta, tsv -> [ groupKey(meta.unit_of, meta.n_units as Integer), meta.unit, tsv ] }
            .groupTuple()
            .map { id, units, tsvs -> [ id.toString(), [ units, tsvs ].transpose().sort { u -> u[0] }.collect { u -> u[1] } ] }
    )
    // Unreadable samples get a detection row saying so, so every sample
    // still reaches the joins below (and the report) exactly once.
    ch_detect = ch_unit_tsv.whole
        .map { meta, tsv -> [ meta.id, tsv ] }
        .mix(COMBINE_DETECTION.out.tsv)
        .map { id, tsv -> [ id, tsv.splitCsv(header: true, sep: '\t')[0] ] }
        .mix(ch_integrity.unreadable.map { meta, problem, n_files -> [ meta, problem, n_files ] }
            .mix(ch_damaged_input.map { meta, files, problem -> [ meta, problem, files.size() ] })
            .map { meta, problem, n_files ->
                [ meta.id, [ n_files: "${n_files}", n_reads: '0', detected_platform: 'uncertain',
                             detected_layout: 'uncertain', confidence: 'uncertain',
                             detection_flags: 'FILE_UNREADABLE', reasons: problem ] ]
            })

    // 4. Join everything onto the sample by ID. failOnMismatch: a sample
    //    missing from either side stops the run instead of silently
    //    dropping out. The screen's NORMALISED declared_platform replaces
    //    the raw value from nf-schema (the report reads the raw value and the
    //    flags straight from the screen output).
    //    meta carries ONLY what downstream tasks need: meta is part of every
    //    task's cache key, so anything extra here (such as the flag count)
    //    makes a wording change in the screen re-run QC. The first -resume
    //    showed exactly that.
    ch_triaged = ch_reads
        .mix(ch_damaged_input.map { meta, files, _problem -> [ meta, files ] })
        .map { meta, reads -> [ meta.id, meta, reads ] }
        .join(ch_screen, failOnMismatch: true, failOnDuplicate: true)
        .join(ch_detect, failOnMismatch: true, failOnDuplicate: true)
        .map { _id, meta, reads, scr, det ->
            // BAM header platform (D55): stands in for a blank samplesheet
            // platform; a samplesheet that disagrees with it is flagged.
            def declared      = blankToNull(scr.declared_platform)
            def declared_from = declared ? 'samplesheet' : null
            def rg            = rgPlatform(meta.bam_rg_platform as String)
            def header_note   = null
            def header_flag   = null
            if (!declared && rg) {
                declared      = rg
                declared_from = 'bam_header'
                header_note   = "No platform in the samplesheet; the BAM header names ${rg} (@RG PL:${meta.bam_rg_platform}), so that is used as the declared platform."
            } else if (declared && rg && rg != declared) {
                header_flag   = 'BAM_HEADER_DISAGREES'
                header_note   = "The samplesheet says ${declared}, but the BAM header says ${rg} (@RG PL:${meta.bam_rg_platform})."
            }
            def merged = meta + [
                declared_platform    : declared,
                declared_from        : declared_from,
                declared_type        : blankToNull(scr.declared_type),
                n_files              : det.n_files as Integer,
                n_reads              : det.n_reads as Long,
                detected_platform    : blankToNull(det.detected_platform),
                detected_layout      : blankToNull(det.detected_layout),
                confidence           : blankToNull(det.confidence),
                detection_flags      : blankToNull([ det.detection_flags, header_flag ].findAll { x -> x }.join(',')),
                reasons              : blankToNull([ meta.bam_note, header_note,
                                                     meta.from_fasta ? 'Reads were given as FASTA, which has no quality values; seq-triage wrote a placeholder quality, so quality is not judged and the sample cannot PASS.' : null,
                                                     det.reasons ].findAll { x -> x }.join(' ')),
                // Read-start concentration for the sample-type call (D65).
                start_share          : blankToNull(det.start_share),
                probe_reads          : blankToNull(det.probe_reads),
                // One genome or many, from the reads alone (D86).
                gc_excess            : blankToNull(det.gc_excess),
                kmer_repeated        : blankToNull(det.kmer_repeated),
                kmer_bases           : blankToNull(det.kmer_bases),
                kmer_genome_mb       : blankToNull(det.kmer_genome_mb),
            ]
            [ routeSample(merged), reads ]
        }

    // 5. Route. Only samples with a platform get QC; 'empty' and
    //    'unrouted' go straight to the report. A sample in several files
    //    headed for QC is merged into one file (or one pair) first, with
    //    nf-core's CAT_FASTQ (D75): QC, the screen and the type call judge
    //    the whole sample. Its files were each checked above.
    ch_to_merge = ch_triaged.branch { meta, _reads ->
        merge: meta.n_units && readClass(meta.route as String)
        as_is: true
    }
    CAT_FASTQ(ch_to_merge.merge)
    //    Library structure (D78), for samples headed for QC: an interleaved
    //    file is split into R1 and R2; a cell-barcode R1 (10x-style) or an
    //    index read given as R2 is set aside, so QC judges the read of the
    //    sample's molecules, as single-end. meta changes only for these.
    ch_shaped = ch_to_merge.as_is.mix(CAT_FASTQ.out.reads).branch { meta, _reads ->
        def f = (meta.detection_flags ?: '').split(',')
        interleaved: readClass(meta.route as String) == 'short' && meta.detected_layout == 'interleaved'
        one_mate:    readClass(meta.route as String) == 'short' && (f.contains('BARCODE_READ') || f.contains('INDEX_READ_AS_MATE'))
        as_is:       true
    }
    DEINTERLEAVE(ch_shaped.interleaved)
    ch_qc_reads = ch_shaped.as_is
        .mix(DEINTERLEAVE.out.reads.map { meta, reads -> [ meta + [ single_end: false ], reads ] })
        .mix(ch_shaped.one_mate.map { meta, reads ->
            def pair = [ reads ].flatten()
            def keep = (meta.detection_flags ?: '').split(',').contains('BARCODE_READ') ? pair[1] : pair[0]
            [ meta + [ single_end: true, qc_mate: keep == pair[1] ? 'R2' : 'R1' ], [ keep ] ]
        })
    ch_routed = ch_qc_reads.branch { meta, _reads ->
        short_reads: readClass(meta.route as String) == 'short'    // illumina, dnbseq, element
        ont:      meta.route == 'ont'
        pacbio:   meta.route == 'pacbio'
        no_qc:    true
    }

    // 6. QC. fastp's three flags: don't discard trimmed reads, don't save
    //    failed reads, don't merge pairs. This is triage, not preprocessing.
    // MGI/DNBSEQ samples get MGI's adapter sequences, so trimming does not
    // depend on fastp's auto-detection (D69); Illumina and Element (Adept,
    // Illumina-compatible adapters) keep auto-detection.
    def mgi_adapters = file("${projectDir}/assets/adapters/mgi.fasta", checkIfExists: true)
    FASTP(ch_routed.short_reads.map { meta, reads -> [ meta, reads, meta.route == 'dnbseq' ? mgi_adapters : [] ] }, false, false, false)
    FASTQC(ch_routed.short_reads)
    // PacBio HiFi gets the same long-read QC as nanopore (D64); its limits
    // differ only in the metric guide (pacbio_reads rows).
    NANOPLOT(ch_routed.ont.mix(ch_routed.pacbio))

    // 6b. Contamination screen (CLAUDE.md D41), only when a Kraken2 database
    //     is given. Every sample sent to QC is classified by Kraken2 (host
    //     and unclassified fractions, species for ONT). Illumina samples
    //     also get Bracken's species-level re-estimate; it is calibrated for
    //     a fixed read length, so it is not used on nanopore reads.
    ch_kraken_reports = channel.empty()
    ch_bracken_tables = channel.empty()
    if (params.kraken2_db) {
        if (params.kraken2_db.endsWith('.tar.gz') || params.kraken2_db.endsWith('.tgz')) {
            UNTAR(channel.of([ [ id: 'kraken2_db' ], file(params.kraken2_db, checkIfExists: true) ]))
            ch_db = UNTAR.out.untar.map { _meta, dir -> dir }.first()
        } else {
            ch_db = channel.value(file(params.kraken2_db, checkIfExists: true))
        }
        ch_to_classify = ch_routed.short_reads.mix(ch_routed.ont).mix(ch_routed.pacbio)
        if (params.kraken2_batch) {
            // One task for every sample: the database is read once (D50).
            // Files are staged as in1/, in2/, ... (several samples may share
            // a file), and the manifest maps each sample to its files.
            ch_batch = ch_to_classify.toList()
                .filter { items -> items.size() > 0 }
                .map { items ->
                    def files = []
                    def lines = items.collect { meta, reads ->
                        def staged = [ reads ].flatten().collect { r -> files << r; "in${files.size()}/${r.name}" }
                        [ meta.id, meta.single_end ? 'true' : 'false', staged.join(',') ].join('\t')
                    }
                    [ items.collect { item -> item[0] }, lines.join('\n'), files ]
                }
            KRAKEN2_BATCH(ch_batch, ch_db)
            ch_kraken_reports = KRAKEN2_BATCH.out.reports.flatMap { metas, reports ->
                [ reports ].flatten().collect { r -> [ metas.find { m -> r.name == "${m.id}.kraken2.report.txt" }, r ] }
            }
        } else {
            KRAKEN2(ch_to_classify, ch_db, false, false)
            ch_kraken_reports = KRAKEN2.out.report
        }
        BRACKEN(ch_kraken_reports.filter { meta, report -> readClass(meta.route as String) == 'short' && maxSpeciesReads(report) >= 10 }, ch_db)
        ch_bracken_tables = BRACKEN.out.reports
    }

    // 6c. rRNA screen (D81), only with --rrna_db: SortMeRNA on samples
    //     declared molecule = rna on short-read routes. References may be
    //     one FASTA or a .txt manifest of several; .gz ones are unzipped
    //     first. The index is built once and reused by every sample.
    ch_rrna_logs = channel.empty()
    if (params.rrna_db) {
        def refs = params.rrna_db.endsWith('.txt')
            ? file(params.rrna_db, checkIfExists: true).readLines().findAll { l -> l.trim() && !l.startsWith('#') }.collect { l -> file(l.trim(), checkIfExists: true) }
            : [ file(params.rrna_db, checkIfExists: true) ]
        ch_refs_in = channel.fromList(refs).branch { f ->
            gz:    f.name.endsWith('.gz')
            plain: true
        }
        GUNZIP_RRNA(ch_refs_in.gz.map { f -> [ [ id: f.baseName ], f ] })
        ch_refs = ch_refs_in.plain.mix(GUNZIP_RRNA.out.gunzip.map { _m, f -> f })
            .collect().map { fs -> [ [ id: 'rrna_db' ], fs.sort { f -> f.name } ] }
        SORTMERNA_INDEX(channel.of([ [ id: 'rrna_index', single_end: true ], [] ]), ch_refs, [ [], [] ])
        ch_rna = ch_routed.short_reads.filter { meta, _reads -> (meta.declared_molecule ?: '').toString().toLowerCase() == 'rna' }
        SORTMERNA(ch_rna, ch_refs.first(), SORTMERNA_INDEX.out.index.first())
        ch_rrna_logs = SORTMERNA.out.log
    }

    ch_multiqc_files = channel.empty()
        .mix(FASTP.out.json.map    { _meta, f -> f })
        .mix(FASTQC.out.zip.map    { _meta, f -> f })
        .mix(NANOPLOT.out.txt.map  { _meta, f -> f })
        .mix(ch_kraken_reports.map { _meta, f -> f })
        .collect()
    MULTIQC(ch_multiqc_files.map { files -> [ [ id: 'seq-triage' ], files, [], [], [], [] ] })

    // 7. Report over ALL samples, routed or not. samples.tsv is built from
    //    meta, one line per sample, sorted so the file is reproducible.
    def cols = [ 'id', 'single_end', 'n_files', 'n_reads', 'detected_platform', 'detected_layout',
                 'confidence', 'detection_flags', 'reasons', 'route', 'routed_by',
                 'declared_platform', 'declared_from', 'start_share', 'probe_reads',
                 'gc_excess', 'kmer_repeated', 'kmer_bases', 'kmer_genome_mb' ]
    ch_samples_tsv = ch_triaged
        .map { meta, _reads -> cols.collect { c -> meta[c] == null ? '' : meta[c].toString() }.join('\t') }
        .collectFile(name: 'samples.tsv', seed: (['sample'] + cols.drop(1)).join('\t'), newLine: true, sort: true)

    //    QC values from every tool, as one long table (helpers above).
    //    seqkit runs on every sample, so this channel is never empty and the
    //    report always runs, even if no sample was routed to QC.
    ch_qc_values = channel.empty()
        .mix(SEQKIT_STATS.out.stats.flatMap { meta, f -> seqkitRows((meta.unit_of ?: meta.id) as String, f) })
        .mix(FASTP.out.json.flatMap         { meta, f -> fastpRows(meta.id, f) })
        .mix(NANOPLOT.out.txt.flatMap       { meta, f -> [ f ].flatten().findAll { x -> x.name.endsWith('NanoStats.txt') }.collectMany { x -> nanoplotRows(meta.id, x) } })
        .mix(FASTQC.out.zip.flatMap         { meta, f -> [ f ].flatten().collectMany { x -> fastqcRows(meta.id, x) } })
        .mix(ch_kraken_reports.flatMap      { meta, f -> kraken2Rows(meta.id, f) })
        .mix(ch_bracken_tables.flatMap      { meta, f -> brackenRows(meta.id, f) })
        .mix(ch_rrna_logs.flatMap           { meta, f -> sortmernaRows(meta.id, f) })
        .collectFile(name: 'qc_values.tsv', seed: "sample\ttool\tfile\tfield\tvalue", newLine: true, sort: true)

    TRIAGE_REPORT(SCREEN_METADATA.out.tsv, ch_samples_tsv, metric_guide, ch_qc_values, species_groups, genome_sizes)

    // 7b. Next pipeline for each ready sample (D66), with a ready-to-run
    //     samplesheet and params file in that pipeline's own format.
    RECOMMEND(TRIAGE_REPORT.out.tsv, TRIAGE_REPORT.out.metrics, samplesheet, file(params.downstream_pipelines))

    // 8. Software versions actually used (reproducibility layer 4): every
    //    module reports [process, tool, version] to the 'versions' topic.
    //    MultiQC reports separately: it consumes outputs, so publishing to the
    //    topic itself could make the pipeline wait forever.
    channel.topic('versions')
        .mix(MULTIQC.out.versions)
        .map { process, tool, version -> "${process.tokenize(':')[-1]}\t${tool}\t${version}" }
        .unique()
        .collectFile(name: 'software_versions.tsv', storeDir: "${params.outdir}/pipeline_info",
                     seed: "process\ttool\tversion", newLine: true, sort: true)

    emit:
    report = TRIAGE_REPORT.out.md
}
