/*
 * SEQ_TRIAGE: screen -> detect -> join -> route -> QC -> report
 *
 * The channel diagram in CLAUDE.md ("Architecture") is the map for this
 * file. Our own logic is the joins, routeSample() and the branch; every
 * standard tool is an nf-core module.
 */

include { SCREEN_METADATA } from '../modules/local/screen_metadata/main'
include { DETECT_PLATFORM } from '../modules/local/detect_platform/main'
include { TRIAGE_REPORT   } from '../modules/local/triage_report/main'
include { SEQKIT_STATS    } from '../modules/nf-core/seqkit/stats/main'
include { FASTP           } from '../modules/nf-core/fastp/main'
include { FASTQC          } from '../modules/nf-core/fastqc/main'
include { NANOPLOT        } from '../modules/nf-core/nanoplot/main'
include { MULTIQC         } from '../modules/nf-core/multiqc/main'

// Empty TSV cells arrive as "", so convert them to null: one meaning of
// "not present" throughout the pipeline (CLAUDE.md D22).
def blankToNull(String s) {
    return (s == null || s.trim().isEmpty()) ? null : s
}

/*
 * routeSample: decide which QC a sample gets, and on what basis (CLAUDE.md D2).
 *   no reads                  -> route 'empty'    (FAIL; no QC can run)
 *   high-confidence detection -> route by DETECTED platform
 *   otherwise, if declared    -> route by DECLARED platform (verdict capped at WARN)
 *   otherwise                 -> 'unrouted'      (UNCERTAIN; no QC run)
 * Returns meta plus route ('illumina' | 'ont' | 'empty' | 'unrouted') and
 * routed_by ('detected' | 'declared' | 'none').
 */
def routeSample(Map meta) {
    if (meta.n_reads == 0) {
        return meta + [ route: 'empty', routed_by: 'none' ]
    }
    if (meta.confidence == 'high') {
        return meta + [ route: meta.detected_platform, routed_by: 'detected' ]
    }
    if (meta.declared_platform) {
        return meta + [ route: meta.declared_platform, routed_by: 'declared' ]
    }
    return meta + [ route: 'unrouted', routed_by: 'none' ]
}

workflow SEQ_TRIAGE {

    take:
    ch_samples      // [ meta, [ reads ] ] from the samplesheet
    samplesheet     // the raw samplesheet file, for the metadata screen
    metric_guide    // assets/metric_guide.tsv or the user's own

    main:

    // 1. Metadata screen: one task for the whole sheet, then one
    //    [ id, row ] element per sample so it can be joined.
    SCREEN_METADATA(samplesheet)
    ch_screen = SCREEN_METADATA.out.tsv
        .splitCsv(header: true, sep: '\t')
        .map { row -> [ row.sample, row ] }

    // 2. Read statistics and platform detection, per sample.
    SEQKIT_STATS(ch_samples)
    DETECT_PLATFORM(SEQKIT_STATS.out.stats)
    ch_detect = DETECT_PLATFORM.out.tsv
        .map { meta, tsv -> [ meta.id, tsv.splitCsv(header: true, sep: '\t')[0] ] }

    // 3. Join everything onto the sample by ID. failOnMismatch: a sample
    //    missing from either side stops the run instead of silently
    //    dropping out. The screen's NORMALISED declared_platform replaces
    //    the raw value from nf-schema (the report reads the raw value and the
    //    flags straight from the screen output).
    //    meta carries ONLY what downstream tasks need: meta is part of every
    //    task's cache key, so anything extra here (such as the flag count)
    //    makes a wording change in the screen re-run QC. The first -resume
    //    showed exactly that.
    ch_triaged = ch_samples
        .map { meta, reads -> [ meta.id, meta, reads ] }
        .join(ch_screen, failOnMismatch: true, failOnDuplicate: true)
        .join(ch_detect, failOnMismatch: true, failOnDuplicate: true)
        .map { _id, meta, reads, scr, det ->
            def merged = meta + [
                declared_platform    : blankToNull(scr.declared_platform),
                declared_type        : blankToNull(scr.declared_type),
                n_files              : det.n_files as Integer,
                n_reads              : det.n_reads as Long,
                detected_platform    : blankToNull(det.detected_platform),
                detected_layout      : blankToNull(det.detected_layout),
                confidence           : blankToNull(det.confidence),
                detection_flags      : blankToNull(det.detection_flags),
                reasons              : blankToNull(det.reasons),
            ]
            [ routeSample(merged), reads ]
        }

    // 4. Route. Only samples with a platform get QC; 'empty' and
    //    'unrouted' go straight to the report.
    ch_routed = ch_triaged.branch { meta, _reads ->
        illumina: meta.route == 'illumina'
        ont:      meta.route == 'ont'
        no_qc:    true
    }

    // 5. QC. fastp's three flags: don't discard trimmed reads, don't save
    //    failed reads, don't merge pairs. This is triage, not preprocessing.
    FASTP(ch_routed.illumina.map { meta, reads -> [ meta, reads, [] ] }, false, false, false)
    FASTQC(ch_routed.illumina)
    NANOPLOT(ch_routed.ont)

    ch_multiqc_files = channel.empty()
        .mix(FASTP.out.json.map    { _meta, f -> f })
        .mix(FASTQC.out.zip.map    { _meta, f -> f })
        .mix(NANOPLOT.out.txt.map  { _meta, f -> f })
        .collect()
    MULTIQC(ch_multiqc_files.map { files -> [ [ id: 'seq-triage' ], files, [], [], [], [] ] })

    // 6. Report over ALL samples, routed or not. samples.tsv is built from
    //    meta, one line per sample, sorted so the file is reproducible.
    def cols = [ 'id', 'single_end', 'n_files', 'n_reads', 'detected_platform', 'detected_layout',
                 'confidence', 'detection_flags', 'reasons', 'route', 'routed_by' ]
    ch_samples_tsv = ch_triaged
        .map { meta, _reads -> cols.collect { c -> meta[c] == null ? '' : meta[c].toString() }.join('\t') }
        .collectFile(name: 'samples.tsv', seed: (['sample'] + cols.drop(1)).join('\t'), newLine: true, sort: true)
    TRIAGE_REPORT(SCREEN_METADATA.out.tsv, ch_samples_tsv, metric_guide)

    // 7. Software versions actually used (reproducibility layer 4): every
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
