#!/usr/bin/env Rscript
# Test for bin/triage_report.R
#
# Builds SYNTHETIC report inputs (metadata screen, per-sample state, QC
# values) in a temporary folder, one sample per rule, runs the report with
# the real assets/metric_guide.tsv, and checks each sample's verdict and the
# result of the metric the case is about. Also checks that the report
# refuses to run when a sample sent to QC has no QC numbers.
#
# Run from the repository root:
#   Rscript tests/triage_report/test_triage_report.R
# Exits non-zero if any check differs, so CI can use it.

tmp <- tempfile("report_test_")
dir.create(tmp)
write_tsv <- function(x, f) write.table(x, f, sep = "\t", quote = FALSE, row.names = FALSE, na = "")

# Healthy QC values; each case overrides what it is about.
FASTP_OK <- c("summary.after_filtering.total_reads" = 500000, "summary.after_filtering.q30_rate" = 0.9,
              "filtering_result.passed_filter_reads" = 495000, "summary.before_filtering.total_reads" = 500000,
              "duplication.rate" = 0.05, "adapter_cutting.adapter_trimmed_reads" = 1000,
              "summary.after_filtering.gc_content" = 0.5, "insert_size.peak" = 300,
              "summary.after_filtering.read1_mean_length" = 150)
NANO_OK  <- c(number_of_reads = 20000, number_of_bases = 1e8, n50 = 10000, mean_qual = 15, median_read_length = 5000)
SEQKIT_OK <- c(num_seqs = 250000, sum_len = 75000000, sum_n = 1000)

# Kraken2 report lines as the workflow flattens them: "rank|taxid|name" = clade reads.
# A clean E. coli isolate: 1,000,000 reads, 2% unclassified, a little Shigella.
KRAKEN_OK <- c("U|0|unclassified" = 20000, "R|1|root" = 980000, "S|562|Escherichia coli" = 900000,
               "S|623|Shigella flexneri" = 30000, "S|573|Klebsiella pneumoniae" = 5000, "S|9606|Homo sapiens" = 2000)

# D67: a virus on its own (taxonomy says Viruses), and a sample whose reads
# are mostly unclassified with one recognised species.
VIRUS <- c("D|10239|Viruses" = 950000, "S|694009|Betacoronavirus pandemicum" = 950000, "S|562|Escherichia coli" = 0,
           "S|623|Shigella flexneri" = 0, "S|573|Klebsiella pneumoniae" = 0, "S|9606|Homo sapiens" = 0)
# The same virus as the Oct 2025 databases report it: Viruses at rank R1
# (an "acellular root" since NCBI's 2024 change), the virus at rank S1.
VIRUS_R1 <- c("R1|10239|Viruses" = 950000, "S|694009|Betacoronavirus pandemicum" = 950000,
              "S1|2697049|Severe acute respiratory syndrome coronavirus 2" = 940000, "S|562|Escherichia coli" = 0,
              "S|623|Shigella flexneri" = 0, "S|573|Klebsiella pneumoniae" = 0, "S|9606|Homo sapiens" = 0)
MOSTLY_UNKNOWN <- c("U|0|unclassified" = 940000, "S|562|Escherichia coli" = 55000, "S|623|Shigella flexneri" = 0,
                    "S|573|Klebsiella pneumoniae" = 0, "S|9606|Homo sapiens" = 0)

# Six species, none above 25%: a community (D65).
META <- c("S|562|Escherichia coli" = 200000, "S|623|Shigella flexneri" = 0, "S|573|Klebsiella pneumoniae" = 180000,
          "S|1280|Staphylococcus aureus" = 170000, "S|287|Pseudomonas aeruginosa" = 160000,
          "S|1352|Enterococcus faecium" = 150000, "S|1313|Streptococcus pneumoniae" = 120000,
          "G|561|Escherichia" = 200000, "G|570|Klebsiella" = 180000, "G|1279|Staphylococcus" = 170000,
          "G|286|Pseudomonas" = 160000, "G|1350|Enterococcus" = 150000, "G|1301|Streptococcus" = 120000)
# D73: one new phage split over three sister species of one genus (layer E
# round 2): species say "many organisms", the genus says one. Human reads in
# genus Homo must not count.
PHAGE_SPLIT <- c("U|0|unclassified" = 100000, "R1|10239|Viruses" = 900000,
                 "G|2732974|Drulisvirus" = 900000, "G|9605|Homo" = 5000000,
                 "S|2560502|Drulisvirus KpS2" = 300000, "S|2560503|Drulisvirus F19" = 250000,
                 "S|2560504|Drulisvirus KpV71" = 200000, "S|2560505|Drulisvirus KpV41" = 150000,
                 "S|562|Escherichia coli" = 0, "S|623|Shigella flexneri" = 0, "S|573|Klebsiella pneumoniae" = 0,
                 "S|9606|Homo sapiens" = 0)
# D76: a human sample (human reads at species and genus level, a trace of E. coli).
HUMAN <- c("U|0|unclassified" = 20000, "S|9606|Homo sapiens" = 950000, "G|9605|Homo" = 950000,
           "S|562|Escherichia coli" = 1000, "S|623|Shigella flexneri" = 0, "S|573|Klebsiella pneumoniae" = 0)
# Community species without any genus rows: D73 cannot run, so not high.
META_NO_GENUS <- META[!startsWith(names(META), "G|")]

# Cases: platform route, metadata, detection, QC overrides (NA removes a
# field), and what to expect.
cases <- list(
  list(id = "ILL_GOOD",          route = "illumina", verdict = "PASS"),
  list(id = "ILL_LOW_READS",     route = "illumina", fastp = c("summary.after_filtering.total_reads" = 50000),
       verdict = "WARN", check = c(total_reads = "WARN")),
  list(id = "ILL_AT_WARN_LIMIT", route = "illumina", fastp = c("summary.after_filtering.total_reads" = 100000),
       verdict = "PASS", check = c(total_reads = "PASS")),
  list(id = "ILL_Q30_FAIL",      route = "illumina", fastp = c("summary.after_filtering.q30_rate" = 0.6),
       verdict = "FAIL", check = c(q30_rate = "FAIL")),
  list(id = "ILL_DUP_DNA",       route = "illumina", fastp = c("duplication.rate" = 0.6),
       type = "isolate", molecule = "dna",              # fully declared: no hint
       verdict = "WARN", check = c(duplication_rate = "WARN")),
  list(id = "ILL_DUP_UNLABELLED", route = "illumina", fastp = c("duplication.rate" = 0.6),
       verdict = "WARN", check = c(duplication_rate = "WARN"), hint = TRUE),
  list(id = "ILL_DUP_RNA",       route = "illumina", fastp = c("duplication.rate" = 0.6), molecule = "rna",
       verdict = "PASS", check = c(duplication_rate = "exempt")),
  list(id = "ILL_NO_ADAPTER",    route = "illumina", fastp = c("adapter_cutting.adapter_trimmed_reads" = NA, "insert_size.peak" = NA),
       verdict = "PASS", check = c(adapter_trimmed_pct = "not_available", insert_size_peak = "not_available")),
  list(id = "ILL_N_RATE",        route = "illumina", seqkit = c(sum_n = 1500000),
       verdict = "WARN", check = c(n_base_rate = "WARN")),
  list(id = "ILL_FASTQC_INFO",   route = "illumina", fastqc = "FAIL",
       verdict = "PASS", check = c(fastqc_per_base_sequence_content = "info")),
  list(id = "ONT_GOOD",          route = "ont", verdict = "PASS"),
  list(id = "ONT_AMPLICON",      route = "ont", nano = c(n50 = 600), type = "amplicon",
       verdict = "PASS", check = c(read_n50 = "exempt")),
  list(id = "ONT_UNLABELLED",    route = "ont", nano = c(n50 = 600),
       verdict = "FAIL", check = c(read_n50 = "FAIL"), hint = TRUE),
  list(id = "ONT_PLACEHOLDER_Q", route = "ont", nano = c(mean_qual = 40), flags = "QUALITY_UNINFORMATIVE",
       confidence = "medium", routed_by = "declared",
       verdict = "WARN", check = c(mean_read_quality = "not_meaningful")),
  list(id = "CEILING_WARN_QC_FAIL", route = "ont", nano = c(number_of_reads = 500), confidence = "medium",
       routed_by = "declared", verdict = "FAIL", check = c(number_of_reads = "FAIL")),
  list(id = "ONT_U_DECLARED_DNA", route = "ont", flags = "U_BASES", molecule = "dna", verdict = "WARN"),
  list(id = "UNROUTED",          route = "unrouted", routed_by = "none", confidence = "uncertain",
       declared = NA, verdict = "UNCERTAIN"),
  # D45: a damaged file is a FAIL row, not a stopped run.
  list(id = "UNREADABLE",        route = "unreadable", routed_by = "none", confidence = "uncertain",
       flags = "FILE_UNREADABLE", verdict = "FAIL"),
  # Decision C1 (validation round 1): headers confirm ONT, but quality is a
  # placeholder, so it cannot be assessed: WARN, never PASS.
  list(id = "ONT_PLACEHOLDER_HIGH", route = "ont", nano = c(mean_qual = 40), flags = "QUALITY_UNINFORMATIVE",
       verdict = "WARN", check = c(mean_read_quality = "not_meaningful")),
  # Contamination screen (D41). Shigella counts with E. coli (species_groups.tsv).
  list(id = "SCREEN_CLEAN",       route = "illumina", screen = TRUE, verdict = "PASS",
       check = c(host_fraction = "PASS", off_target_fraction = "PASS", dominant_species = "info")),
  list(id = "SCREEN_HOST_WARN",   route = "ont", screen = TRUE, kraken = c("S|9606|Homo sapiens" = 50000),
       type = "isolate", verdict = "WARN", check = c(host_fraction = "WARN")),          # declared: no hint
  list(id = "SCREEN_HOST_FAIL",   route = "ont", screen = TRUE, kraken = c("S|9606|Homo sapiens" = 200000),
       verdict = "FAIL", check = c(host_fraction = "FAIL"), hint = TRUE),              # blank type: metagenome hint
  list(id = "SCREEN_HOST_META",   route = "ont", screen = TRUE, kraken = c("S|9606|Homo sapiens" = 200000), type = "metagenome",
       verdict = "PASS", check = c(host_fraction = "exempt")),
  list(id = "SCREEN_MIXED",       route = "ont", screen = TRUE, kraken = c("S|573|Klebsiella pneumoniae" = 300000),
       verdict = "FAIL", check = c(off_target_fraction = "FAIL"), hint = TRUE),
  list(id = "SCREEN_SAME_GENUS",  route = "ont", screen = TRUE, kraken = c("S|564|Escherichia fergusonii" = 90000),
       verdict = "WARN", check = c(off_target_fraction = "WARN"), hint = TRUE),
  list(id = "SCREEN_BRACKEN",     route = "illumina", screen = TRUE, bracken = c("S|562|Escherichia coli" = 700000, "S|573|Klebsiella pneumoniae" = 300000),
       verdict = "FAIL", check = c(off_target_fraction = "FAIL"), hint = TRUE),
  list(id = "SCREEN_NO_SPECIES",  route = "ont", screen = TRUE, kraken = c("U|0|unclassified" = 999995, "R|1|root" = 5, "S|562|Escherichia coli" = 5, "S|623|Shigella flexneri" = 0, "S|573|Klebsiella pneumoniae" = 0, "S|9606|Homo sapiens" = 0),
       verdict = "PASS", check = c(off_target_fraction = "not_available", unclassified_fraction = "info")),
  # D47: trace species (< 1% each) do not count as off-target, however many;
  # a real contaminant above 1% still does. 40 trace species of 1,000 reads
  # each (~0.1% each, ~4% in total) on a clean culture:
  list(id = "SCREEN_TRACE_TAIL",   route = "ont", screen = TRUE, type = "isolate",
       kraken = setNames(rep(1000, 40), sprintf("S|%d|Trace species %d", 900001:900040, 1:40)),
       verdict = "PASS", check = c(off_target_fraction = "PASS")),
  list(id = "SCREEN_TAIL_PLUS_REAL", route = "ont", screen = TRUE, type = "isolate",
       kraken = c(setNames(rep(1000, 40), sprintf("S|%d|Trace species %d", 900001:900040, 1:40)),
                  "S|573|Klebsiella pneumoniae" = 60000),
       verdict = "WARN", check = c(off_target_fraction = "WARN")),
  # D40: a platform seq-triage names but has no QC for.
  list(id = "UNSUPPORTED_PACBIO", route = "unrouted", routed_by = "unsupported", detected = "pacbio",
       declared = NA, flags = "PLATFORM_NOT_HANDLED", verdict = "UNCERTAIN"),
  # D55: BAM header @RG PL. Blank samplesheet platform, header names it, reads
  # confirm: the header counts as the declaration, so PASS is possible.
  list(id = "BAM_HEADER_DECLARED", route = "illumina", declared = NA, hdr_declared = "illumina",
       declared_from = "bam_header", verdict = "PASS"),
  # Same reads, no platform anywhere: inferred, at most WARN (control).
  list(id = "BAM_NO_PLATFORM",     route = "illumina", declared = NA, verdict = "WARN"),
  # D64: PacBio HiFi QC. Same NanoPlot fields as nanopore, own quality limit
  # (Q20 / Q10): Q15 passes for nanopore but WARNs for HiFi, which shows the
  # report uses the platform's own guide row, not the first row of that name.
  list(id = "HIFI_GOOD",       route = "pacbio", declared = "pacbio", nano = c(mean_qual = 25),
       verdict = "PASS", check = c(mean_read_quality = "PASS")),
  list(id = "HIFI_Q15",        route = "pacbio", declared = "pacbio", nano = c(mean_qual = 15),
       verdict = "WARN", check = c(mean_read_quality = "WARN")),
  list(id = "ONT_Q15",         route = "ont", nano = c(mean_qual = 15),
       verdict = "PASS", check = c(mean_read_quality = "PASS")),
  list(id = "HIFI_INFERRED",   route = "pacbio", declared = NA, nano = c(mean_qual = 25), verdict = "WARN"),
  list(id = "PACBIO_NOT_HIFI", route = "unrouted", routed_by = "pacbio_not_hifi", declared = "pacbio",
       detected = "uncertain", confidence = "uncertain", verdict = "UNCERTAIN"),
  # D65: sample type. expect_type = "type confidence"; expect_flag = a flag that must be set.
  # (Not ty/flag: R's $ partially matches them to the existing type/flags fields.)
  list(id = "TYPE_AMP_UNDECLARED", route = "illumina", start = 0.30, verdict = "PASS",
       expect_type = "amplicon high", hint = TRUE),
  list(id = "TYPE_AMP_DECLARED",   route = "illumina", start = 0.30, type = "amplicon", verdict = "PASS", expect_type = "amplicon high"),
  list(id = "TYPE_AMP_SAYS_ISOLATE", route = "illumina", start = 0.30, type = "isolate", verdict = "WARN",
       expect_type = "amplicon high", expect_flag = "TYPE_CONTRADICTS"),
  list(id = "TYPE_RNA_LIKE",       route = "illumina", start = 0.13, verdict = "PASS", expect_type = "uncertain uncertain"),
  list(id = "TYPE_FEW_READS",      route = "illumina", start = 0.30, probe = 100, verdict = "PASS", expect_type = "uncertain uncertain"),
  list(id = "TYPE_NO_SCREEN",      route = "illumina", start = 0.01, verdict = "PASS", expect_type = "uncertain uncertain"),
  list(id = "TYPE_ISOLATE_ILM",    route = "illumina", start = 0.01, type = "isolate", screen = TRUE, verdict = "PASS",
       expect_type = "isolate high"),
  list(id = "TYPE_META",           route = "illumina", start = 0.01, type = "metagenome", screen = TRUE, kraken = META,
       verdict = "PASS", expect_type = "metagenome high"),
  list(id = "TYPE_META_SAYS_ISOLATE", route = "illumina", start = 0.01, type = "isolate", screen = TRUE, kraken = META,
       verdict = "FAIL", expect_type = "metagenome high", expect_flag = "TYPE_CONTRADICTS"),
  list(id = "TYPE_ONT_ISOLATE",    route = "ont", type = "amplicon", screen = TRUE, verdict = "PASS", expect_type = "isolate medium"),
  # Nanopore 16S amplicon of a community: composition says metagenome, but
  # long reads cannot rule out an amplicon, so only medium confidence.
  list(id = "TYPE_ONT_16S",        route = "ont", type = "amplicon", screen = TRUE, kraken = META,
       verdict = "PASS", expect_type = "metagenome medium"),
  # Long reads cannot rule out an amplicon, so a community is metagenome at
  # medium confidence only, and a declared isolate is not contradicted.
  list(id = "TYPE_ONT_META_SAYS_ISOLATE", route = "ont", type = "isolate", screen = TRUE, kraken = META,
       verdict = "FAIL", expect_type = "metagenome medium"),
  # D67 (after validation layer E): fragmented viral amplicons start at random,
  # like an isolate, so a virus-dominated sample is isolate at most MEDIUM;
  # and an isolate call resting on few classified reads is medium.
  list(id = "TYPE_VIRUS_FRAGMENTED", route = "illumina", start = 0.03, screen = TRUE, kraken = VIRUS,
       verdict = "PASS", expect_type = "isolate medium", check = c(viral_fraction = "info")),
  list(id = "TYPE_VIRUS_RANK_R1", route = "illumina", start = 0.03, screen = TRUE, kraken = VIRUS_R1,
       verdict = "PASS", expect_type = "isolate medium"),
  list(id = "TYPE_MOSTLY_UNCLASSIFIED", route = "illumina", start = 0.02, screen = TRUE, kraken = MOSTLY_UNKNOWN,
       verdict = "PASS", expect_type = "isolate medium"),
  # D73 (after layer E round 2): metagenome HIGH needs the top genus under 50%.
  list(id = "TYPE_PHAGE_SPLIT",    route = "illumina", start = 0.02, screen = TRUE, kraken = PHAGE_SPLIT,
       verdict = "FAIL", hint = TRUE, expect_type = "metagenome medium"),
  list(id = "TYPE_META_NO_GENUS",  route = "illumina", start = 0.01, screen = TRUE, kraken = META_NO_GENUS,
       verdict = "FAIL", hint = TRUE, expect_type = "metagenome medium"),
  # D76: the target organism. Human declared as the target: host_fraction not
  # judged, and the sample is one organism.
  list(id = "HUMAN_DECLARED",      route = "illumina", start = 0.01, type = "isolate", target = "9606", screen = TRUE, kraken = HUMAN,
       verdict = "PASS", check = c(host_fraction = "exempt", target_fraction = "info"), expect_type = "isolate high"),
  # Human-dominated, nothing declared: human sample or swamped isolate? UNCERTAIN.
  list(id = "HUMAN_UNDECLARED",    route = "illumina", start = 0.01, screen = TRUE, kraken = HUMAN,
       verdict = "UNCERTAIN", check = c(host_fraction = "exempt"), expect_type = "uncertain uncertain", expect_flag = "HOST_DOMINANT"),
  # A yeast the database cannot see: contamination not judged, cannot PASS.
  list(id = "TARGET_NOT_SEEN",     route = "illumina", start = 0.01, type = "isolate", target = "4932", screen = TRUE, kraken = MOSTLY_UNKNOWN,
       verdict = "WARN", check = c(off_target_fraction = "not_available"), expect_type = "uncertain uncertain", expect_flag = "TARGET_NOT_SEEN"),
  # Declared K. pneumoniae, but the reads are E. coli: off-target is measured
  # against the target, so it FAILs, and the flag names the problem.
  list(id = "TARGET_NOT_DOMINANT", route = "illumina", start = 0.01, type = "isolate", target = "573", screen = TRUE,
       verdict = "FAIL", check = c(off_target_fraction = "FAIL"), expect_flag = "TARGET_NOT_DOMINANT"),
  list(id = "TARGET_OK",           route = "illumina", start = 0.01, type = "isolate", target = "562", screen = TRUE,
       verdict = "PASS", check = c(off_target_fraction = "PASS", target_fraction = "info")),
  # D78: library structure from the reads.
  list(id = "INDEX_READS_ONLY",    route = "unrouted", routed_by = "index_reads", detected = "uncertain", confidence = "uncertain",
       flags = "INDEX_READS", verdict = "UNCERTAIN"),
  list(id = "INDEX_AS_MATE",       route = "illumina", flags = "INDEX_READ_AS_MATE", verdict = "WARN"),
  list(id = "BISULFITE_LIB",       route = "illumina", start = 0.01, type = "isolate", screen = TRUE, kraken = META, flags = "BISULFITE",
       verdict = "PASS", check = c(host_fraction = "exempt", off_target_fraction = "exempt"), expect_type = "uncertain uncertain"),
  # A single-cell R2 can start at few places; it is still not an "amplicon".
  list(id = "SINGLE_CELL_LIB",     route = "illumina", start = 0.30, flags = "BARCODE_READ", verdict = "PASS",
       expect_type = "uncertain uncertain"),
  # D80: depth of coverage = all bases (75,000,000 in SEQKIT_OK) / genome size.
  list(id = "DEPTH_FROM_SPECIES",  route = "illumina", start = 0.01, type = "isolate", screen = TRUE,
       verdict = "PASS", expect_depth = "14.6x"),                       # E. coli: NCBI 5,150,000 bp
  list(id = "DEPTH_DECLARED",      route = "illumina", type = "isolate", gsize = "2000000",
       verdict = "PASS", expect_depth = "37.5x"),
  list(id = "DEPTH_EUKARYOTE",     route = "illumina", start = 0.01, type = "isolate", target = "4932", screen = TRUE, kraken = MOSTLY_UNKNOWN,
       verdict = "WARN", expect_depth = "6.8x"),                        # yeast: NCBI 11,000,000 bp; the database cannot see it
  list(id = "DEPTH_METAGENOME",    route = "illumina", start = 0.01, type = "metagenome", screen = TRUE, kraken = META,
       verdict = "PASS", expect_depth = "not computed: not a single genome"),
  list(id = "DEPTH_RNA",           route = "illumina", type = "isolate", molecule = "rna",
       verdict = "PASS", expect_depth = "not computed: RNA (expression, not genome coverage)"),
  list(id = "DEPTH_NO_SIZE",       route = "illumina", type = "isolate",
       verdict = "PASS", expect_depth = "no genome size: declare genome_size or target_taxid (or run the contamination screen)"),
  # D85: RNA reads that start unevenly are not proof of an amplicon.
  list(id = "RNA_UNEVEN_STARTS",    route = "illumina", start = 0.30, molecule = "rna", verdict = "PASS", expect_type = "amplicon medium"),
  # D81: rRNA share of RNA samples.
  list(id = "RRNA_DEPLETED",       route = "illumina", molecule = "rna", rrna = 0.05, verdict = "PASS", check = c(rrna_fraction = "PASS")),
  list(id = "RRNA_FAILED_DEPLETION", route = "illumina", molecule = "rna", rrna = 0.82, verdict = "FAIL", hint = TRUE, check = c(rrna_fraction = "FAIL")),  # +hint: a 16S amplicon looks like this
  list(id = "RRNA_NOT_SCREENED",   route = "illumina", molecule = "rna", verdict = "PASS", check = c(rrna_fraction = "not_available")),
  list(id = "RRNA_16S_AMPLICON",   route = "illumina", molecule = "rna", type = "amplicon", rrna = 0.97, verdict = "PASS", check = c(rrna_fraction = "exempt")),
  # D86: a metagenome call is high only when the reads look like many genomes.
  list(id = "META_KMERS_RECUR",    route = "illumina", start = 0.01, type = "metagenome", screen = TRUE, kraken = META, krep = "0.65",
       verdict = "PASS", expect_type = "metagenome medium"),        # an archaeon the database scatters
  list(id = "META_NARROW_GC",      route = "illumina", start = 0.01, type = "metagenome", screen = TRUE, kraken = META, gcx = "0.029",
       verdict = "PASS", expect_type = "metagenome medium"),        # a phage
  list(id = "META_TOO_LITTLE_SEQ", route = "illumina", start = 0.01, type = "metagenome", screen = TRUE, kraken = META, kbases = "2500000",
       verdict = "PASS", expect_type = "metagenome medium"),
  list(id = "META_ALL_AGREE",      route = "illumina", start = 0.01, type = "metagenome", screen = TRUE, kraken = META, krep = "0.10", gcx = "0.11",
       verdict = "PASS", expect_type = "metagenome high", check = c(kmer_repeat_share = "info", gc_spread = "info")),
  # D69: MGI/DNBSEQ and Element are QC'd like Illumina (fastp, own guide rows).
  list(id = "DNBSEQ_GOOD",       route = "dnbseq", declared = "dnbseq", verdict = "PASS", check = c(q30_rate = "PASS")),
  list(id = "ELEMENT_LOW_READS", route = "element", declared = "element", fastp = c("summary.after_filtering.total_reads" = 50000),
       verdict = "WARN", check = c(total_reads = "WARN")),
  # D60: a medium call contradicting the declaration is not routed.
  list(id = "CONTRADICTED_MEDIUM", route = "unrouted", routed_by = "contradicted", declared = "ont",
       detected = "illumina", confidence = "medium", verdict = "UNCERTAIN"),
  # Samplesheet and header disagree: two metadata sources conflict, no PASS.
  list(id = "BAM_HEADER_DISAGREES", route = "illumina", flags = "BAM_HEADER_DISAGREES", verdict = "WARN")
)

build <- function(cases, dir) {
  screen <- samples <- qc <- list()
  for (cs in cases) {
    route    <- cs$route
    declared <- if ("declared" %in% names(cs)) cs$declared else if (route %in% c("illumina", "dnbseq", "element", "ont", "pacbio")) route else NA
    screen[[cs$id]] <- data.frame(sample = cs$id, declared_platform_raw = declared, declared_platform = declared,
      declared_type = if (is.null(cs$type)) NA else cs$type,
      declared_molecule = if (is.null(cs$molecule)) NA else cs$molecule, n_flags = "0", flags = NA,
      declared_target = if (is.null(cs$target)) NA else cs$target,
      declared_genome_size = if (is.null(cs$gsize)) NA else cs$gsize)
    samples[[cs$id]] <- data.frame(sample = cs$id, single_end = "false", n_files = "1", n_reads = "1000",
      detected_platform = if (!is.null(cs$detected)) cs$detected else if (route %in% c("illumina", "dnbseq", "element", "ont", "pacbio")) route else "uncertain",
      detected_layout = "single", confidence = if (is.null(cs$confidence)) "high" else cs$confidence,
      detection_flags = if (is.null(cs$flags)) NA else cs$flags, reasons = "synthetic",
      route = route, routed_by = if (is.null(cs$routed_by)) "detected" else cs$routed_by,
      declared_platform = if (!is.null(cs$hdr_declared)) cs$hdr_declared else declared,
      declared_from = if (!is.null(cs$declared_from)) cs$declared_from else if (!is.na(declared)) "samplesheet" else NA,
      start_share = if (is.null(cs$start)) NA else cs$start,
      probe_reads = if (is.null(cs$start)) NA else if (is.null(cs$probe)) 2000 else cs$probe,
      # D86: one genome or many; defaults look like a community.
      gc_excess = if (is.null(cs$gcx)) "0.09" else cs$gcx,
      kmer_repeated = if (is.null(cs$krep)) "0.15" else cs$krep,
      kmer_bases = if (is.null(cs$kbases)) "7500000" else cs$kbases,
      kmer_genome_mb = "1.0")
    add <- function(tool, base, over) {
      v <- base; if (length(over)) v[names(over)] <- over
      v <- v[!is.na(v)]
      data.frame(sample = cs$id, tool = tool, file = "f", field = names(v), value = format(v, scientific = FALSE, trim = TRUE))
    }
    qc[[paste0(cs$id, "s")]] <- add("seqkit", SEQKIT_OK, cs$seqkit)
    if (route %in% c("illumina", "dnbseq", "element")) {
      qc[[paste0(cs$id, "f")]] <- add("fastp", FASTP_OK, cs$fastp)
      st <- if (is.null(cs$fastqc)) "PASS" else cs$fastqc
      qc[[paste0(cs$id, "q")]] <- data.frame(sample = cs$id, tool = "FastQC", file = "f",
        field = c("Per base sequence content", "Overrepresented sequences", "Per sequence GC content"), value = st)
    }
    if (route %in% c("ont", "pacbio")) qc[[paste0(cs$id, "n")]] <- add("NanoPlot", NANO_OK, cs$nano)
    # D81: an rRNA share from SortMeRNA (rrna = the fraction).
    if (!is.null(cs$rrna)) qc[[paste0(cs$id, "r")]] <- data.frame(sample = cs$id, tool = "SortMeRNA", file = "f", field = "rrna_fraction", value = format(cs$rrna))
    if (isTRUE(cs$screen)) {
      kr <- add("Kraken2", KRAKEN_OK, cs$kraken)
      if (!is.null(cs$kraken) && any(names(cs$kraken) != "U|0|unclassified")) {
        # keep root consistent: root = classified total
        v <- as.numeric(kr$value); n <- kr$field
        v[n == "R|1|root"] <- sum(v[grepl("^S\\|", n)])
        kr$value <- format(v, scientific = FALSE, trim = TRUE)
        if ("R|1|root" %in% names(cs$kraken)) kr$value[n == "R|1|root"] <- format(cs$kraken[["R|1|root"]], scientific = FALSE)
      }
      qc[[paste0(cs$id, "k")]] <- kr
      if (!is.null(cs$bracken)) qc[[paste0(cs$id, "b")]] <- add("Bracken", cs$bracken, NULL)
    }
  }
  write_tsv(do.call(rbind, screen),  file.path(dir, "screen.tsv"))
  write_tsv(do.call(rbind, samples), file.path(dir, "samples.tsv"))
  write_tsv(do.call(rbind, qc),      file.path(dir, "qc.tsv"))
}

run_report <- function(dir) {
  system2("Rscript", c("--vanilla", "bin/triage_report.R", file.path(dir, "screen.tsv"), file.path(dir, "samples.tsv"),
                       "assets/metric_guide.tsv", file.path(dir, "qc.tsv"), "assets/species_groups.tsv", file.path(dir, "out"),
                       "assets/species_genome_size.tsv"),
          stdout = FALSE, stderr = FALSE)
}

# Two report runs: without the contamination screen, and with it (when the
# screen is on, EVERY sample sent to QC must have Kraken2 results, so the two
# kinds cannot share one report).
screened <- vapply(cases, function(cs) isTRUE(cs$screen), TRUE)
summ <- mets <- NULL
for (part in list(cases[!screened], cases[screened])) {
  d <- tempfile("run_", tmpdir = tmp); dir.create(d)
  build(part, d)
  status <- run_report(d)
  if (status != 0) stop("triage_report.R exited with status ", status)
  summ <- rbind(summ, read.delim(file.path(d, "out_summary.tsv"), colClasses = "character", na.strings = ""))
  mets <- rbind(mets, read.delim(file.path(d, "out_metrics.tsv"), colClasses = "character", na.strings = ""))
}

res <- do.call(rbind, lapply(cases, function(cs) {
  got_v <- summ$verdict[summ$sample == cs$id]
  got_c <- if (length(cs$check)) paste(names(cs$check), mets$result[match(paste(cs$id, names(cs$check)), paste(mets$sample, mets$metric))], sep = "=", collapse = " ") else ""
  exp_c <- if (length(cs$check)) paste(names(cs$check), cs$check, sep = "=", collapse = " ") else ""
  got_h <- !is.na(summ$hint[summ$sample == cs$id])
  r <- summ[summ$sample == cs$id, ]
  got_t <- if (is.null(cs$expect_type)) "" else paste(r$detected_type, r$type_confidence)
  got_f <- if (is.null(cs$expect_flag)) "" else if (grepl(cs$expect_flag, ifelse(is.na(r$detection_flags), "", r$detection_flags))
) cs$expect_flag else "no flag"
  # D80: expect_depth = the depth_of_coverage value, or its note when not computed.
  got_d <- if (is.null(cs$expect_depth)) "" else {
    m <- mets[mets$sample == cs$id & mets$metric == "depth_of_coverage", ]
    if (!nrow(m)) "no depth row" else if (m$result == "info") m$value else m$note }
  got_c <- paste(got_c, got_d)
  exp_c <- paste(exp_c, if (is.null(cs$expect_depth)) "" else cs$expect_depth)
  ok <- identical(got_v, cs$verdict) && got_c == exp_c && got_h == isTRUE(cs$hint) &&
        got_t == (if (is.null(cs$expect_type)) "" else cs$expect_type) && got_f == (if (is.null(cs$expect_flag)) "" else cs$expect_flag)
  data.frame(case = cs$id, got = paste(got_v, got_c, if (got_h) "+hint" else "", got_t, got_f),
             expected = paste(cs$verdict, exp_c, if (isTRUE(cs$hint)) "+hint" else "", if (is.null(cs$expect_type)) "" else cs$expect_type,
                              if (is.null(cs$expect_flag)) "" else cs$expect_flag),
             result = if (ok) "PASS" else "FAIL")
}))

# Gate: a sample routed to QC with no QC numbers must stop the report.
gate <- file.path(tmp, "gate"); dir.create(gate)
build(list(list(id = "NO_QC", route = "illumina", verdict = "")), gate)
qcg <- read.delim(file.path(gate, "qc.tsv"), colClasses = "character")
write_tsv(qcg[qcg$tool == "seqkit", ], file.path(gate, "qc.tsv"))
gate_status <- run_report(gate)
res <- rbind(res, data.frame(case = "GATE_MISSING_QC", got = sprintf("exit %d", gate_status),
                             expected = "non-zero exit", result = if (gate_status != 0) "PASS" else "FAIL"))

# Gate: screen on, but one sample sent to QC has no Kraken2 results.
gate2 <- file.path(tmp, "gate2"); dir.create(gate2)
build(list(list(id = "HAS_SCREEN", route = "ont", screen = TRUE, verdict = ""),
           list(id = "NO_SCREEN",  route = "ont", verdict = "")), gate2)
gate2_status <- run_report(gate2)
res <- rbind(res, data.frame(case = "GATE_MISSING_KRAKEN2", got = sprintf("exit %d", gate2_status),
                             expected = "non-zero exit", result = if (gate2_status != 0) "PASS" else "FAIL"))

print(res, row.names = FALSE, right = FALSE)
n_fail <- sum(res$result == "FAIL")
cat(sprintf("\n%d / %d cases as expected\n", nrow(res) - n_fail, nrow(res)))
unlink(tmp, recursive = TRUE)
if (n_fail > 0) quit(status = 1)
