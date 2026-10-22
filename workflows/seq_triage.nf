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
include { CHECK_BAM       } from '../modules/local/check_bam/main'
include { SAMTOOLS_COLLATEFASTQ } from '../modules/nf-core/samtools/collatefastq/main'

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
        note = 'Reads extracted from BAM: paired by their READ1/READ2 flags.'
        if (hasReads(other) || hasReads(singleton)) {
            note += ' Some reads had no mate or no READ1/READ2 flag and are not used.'
        }
    } else if (hasReads(other)) {
        out  = [ meta + [ single_end: true ], [ other ] ]
        note = 'Reads extracted from BAM: no READ1/READ2 flags, so treated as one file of reads.'
        if (hasReads(r1) || hasReads(r2)) {
            note += ' Reads flagged as only one mate are not used.'
        }
    } else {
        def first = hasReads(r1) ? r1 : (hasReads(r2) ? r2 : other)
        out  = [ meta + [ single_end: true ], [ first ] ]
        note = hasReads(first) ? 'Reads extracted from BAM: only one mate present, treated as single.'
                               : 'Reads extracted from BAM: the BAM holds no primary reads.'
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
    def known = [ illumina: 'illumina', ont: 'ont', nanopore: 'ont', oxford_nanopore: 'ont' ]
    def hits = pl.split(',').collect { p -> known[p.trim().toLowerCase()] }.unique()
    return hits.size() == 1 ? hits[0] : null
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
 *   otherwise, if declared    -> route by DECLARED platform (verdict capped at WARN)
 *   otherwise                 -> 'unrouted'      (UNCERTAIN; no QC run)
 * Then one safety rule (D33): a sample headed for ONT QC with two files is
 * stopped ('unrouted', routed_by 'blocked'). Nanopore reads are never
 * paired, and NanoPlot fails on two files, which would stop the whole run.
 * A sample confidently identified as a platform with no QC here (PacBio,
 * MGI/DNBSEQ, Element) is also 'unrouted' (routed_by 'unsupported', D40).
 * Returns meta plus route ('illumina' | 'ont' | 'empty' | 'unreadable' | 'unrouted') and
 * routed_by ('detected' | 'declared' | 'blocked' | 'unsupported' | 'none').
 */
def routeSample(Map meta) {
    def r
    if (meta.detection_flags == 'FILE_UNREADABLE') {
        return meta + [ route: 'unreadable', routed_by: 'none' ]
    }
    if (meta.n_reads == 0) {
        r = [ route: 'empty', routed_by: 'none' ]
    } else if (meta.confidence == 'high') {
        r = [ route: meta.detected_platform, routed_by: 'detected' ]
    } else if (meta.declared_platform) {
        r = [ route: meta.declared_platform, routed_by: 'declared' ]
    } else {
        r = [ route: 'unrouted', routed_by: 'none' ]
    }
    if (r.route == 'ont' && meta.n_files > 1) {
        r = [ route: 'unrouted', routed_by: 'blocked' ]
    }
    // A named platform seq-triage has no QC for (PacBio, MGI/DNBSEQ,
    // Element; D40): reported by name, no QC.
    if (r.routed_by == 'detected' && !(r.route in [ 'illumina', 'ont' ])) {
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
    ch_bam = ch_input.bam.map { meta, files -> [ meta, files[0] ] }
    CHECK_BAM(ch_bam)
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
    // run kept 60 of 200 reads without it). No reference: BAM, not CRAM.
    SAMTOOLS_COLLATEFASTQ(ch_bam_check.ok, [ [], [], [] ], false)
    ch_from_bam = SAMTOOLS_COLLATEFASTQ.out.fastq
        .join(SAMTOOLS_COLLATEFASTQ.out.fastq_other, failOnMismatch: true, failOnDuplicate: true)
        .join(SAMTOOLS_COLLATEFASTQ.out.fastq_singleton, failOnMismatch: true, failOnDuplicate: true)
        .map { meta, pair, other, singleton -> bamReads(meta, pair, other, singleton) }
    ch_reads = ch_input.fastq.mix(ch_from_bam)

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
    SEQKIT_STATS(ch_readable)
    DETECT_PLATFORM(
        SEQKIT_STATS.out.stats.join(ch_readable, failOnMismatch: true, failOnDuplicate: true)
    )
    // Unreadable samples get a detection row saying so, so every sample
    // still reaches the joins below (and the report) exactly once.
    ch_detect = DETECT_PLATFORM.out.tsv
        .map { meta, tsv -> [ meta.id, tsv.splitCsv(header: true, sep: '\t')[0] ] }
        .mix(ch_integrity.unreadable.map { meta, problem, n_files -> [ meta, problem, n_files ] }
            .mix(ch_bam_check.damaged.map { meta, _bam, problem -> [ meta, problem, 1 ] })
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
        .mix(ch_bam_check.damaged.map { meta, bam, _problem -> [ meta, [ bam ] ] })
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
                reasons              : blankToNull([ meta.bam_note, header_note, det.reasons ].findAll { x -> x }.join(' ')),
            ]
            [ routeSample(merged), reads ]
        }

    // 5. Route. Only samples with a platform get QC; 'empty' and
    //    'unrouted' go straight to the report.
    ch_routed = ch_triaged.branch { meta, _reads ->
        illumina: meta.route == 'illumina'
        ont:      meta.route == 'ont'
        no_qc:    true
    }

    // 6. QC. fastp's three flags: don't discard trimmed reads, don't save
    //    failed reads, don't merge pairs. This is triage, not preprocessing.
    FASTP(ch_routed.illumina.map { meta, reads -> [ meta, reads, [] ] }, false, false, false)
    FASTQC(ch_routed.illumina)
    NANOPLOT(ch_routed.ont)

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
        ch_to_classify = ch_routed.illumina.mix(ch_routed.ont)
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
        BRACKEN(ch_kraken_reports.filter { meta, report -> meta.route == 'illumina' && maxSpeciesReads(report) >= 10 }, ch_db)
        ch_bracken_tables = BRACKEN.out.reports
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
                 'declared_platform', 'declared_from' ]
    ch_samples_tsv = ch_triaged
        .map { meta, _reads -> cols.collect { c -> meta[c] == null ? '' : meta[c].toString() }.join('\t') }
        .collectFile(name: 'samples.tsv', seed: (['sample'] + cols.drop(1)).join('\t'), newLine: true, sort: true)

    //    QC values from every tool, as one long table (helpers above).
    //    seqkit runs on every sample, so this channel is never empty and the
    //    report always runs, even if no sample was routed to QC.
    ch_qc_values = channel.empty()
        .mix(SEQKIT_STATS.out.stats.flatMap { meta, f -> seqkitRows(meta.id, f) })
        .mix(FASTP.out.json.flatMap         { meta, f -> fastpRows(meta.id, f) })
        .mix(NANOPLOT.out.txt.flatMap       { meta, f -> [ f ].flatten().findAll { x -> x.name.endsWith('NanoStats.txt') }.collectMany { x -> nanoplotRows(meta.id, x) } })
        .mix(FASTQC.out.zip.flatMap         { meta, f -> [ f ].flatten().collectMany { x -> fastqcRows(meta.id, x) } })
        .mix(ch_kraken_reports.flatMap      { meta, f -> kraken2Rows(meta.id, f) })
        .mix(ch_bracken_tables.flatMap      { meta, f -> brackenRows(meta.id, f) })
        .collectFile(name: 'qc_values.tsv', seed: "sample\ttool\tfile\tfield\tvalue", newLine: true, sort: true)

    TRIAGE_REPORT(SCREEN_METADATA.out.tsv, ch_samples_tsv, metric_guide, ch_qc_values, species_groups)

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
