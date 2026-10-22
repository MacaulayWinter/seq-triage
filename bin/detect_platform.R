#!/usr/bin/env Rscript
#*******************************************************************************
# detect_platform.R: infer sequencing platform and layout from the reads
#
# WHAT IT DOES  Calls each sample Illumina, ONT, a named platform seq-triage
#               does not handle (PacBio, MGI/DNBSEQ, Element) or uncertain,
#               with a confidence level, a layout, warning flags and a
#               plain-language reason.
# HOW           1. CORE CALL from read-length shape (seqkit stats, whole
#                  file). The longest read and the length spread must agree;
#                  one alone is never enough (CLAUDE.md D3).
#               2. SUPPORTING EVIDENCE from the first PROBE_N reads: read-name
#                  format, library adapters and the number of distinct
#                  quality values. These raise or lower confidence but never
#                  make an Illumina/ONT call.
#               2b. INSTRUMENT MARKERS (read names, MGI/PacBio adapters) NAME a
#                  platform seq-triage does not handle, if the read lengths
#                  agree (D40).
#               3. LAYOUT from the number of files, their read counts and
#                  their read names.
# RULE          When the evidence is weak or disagrees, say so. "uncertain"
#               is a result, not an error (D1).
# CALIBRATION   Every limit below was checked against 21 real and test files
#               (CLAUDE.md "Test data survey" and D33). Observed ranges are
#               noted next to each limit.
#
# BAM INPUT     (D55) --bam-programs=<names>: the program names (@PG PN)
#               in a BAM's header. Software that only runs on one platform's
#               data is evidence written by the instrument's own pipeline,
#               not typed by a person: PacBio software names PacBio like a
#               read-name marker; basecaller or Illumina software supports
#               or contradicts the call. @RG PL, which people type, is not
#               used here: the workflow treats it as a declared platform.
#
# Usage:  detect_platform.R <sample_id> <seqkit_stats.tsv> <output.tsv> <reads_1> [reads_2]
#                           [--bam-programs=ccs,samtools,...]
#*******************************************************************************

args <- commandArgs(trailingOnly = TRUE)
opt  <- grepl("^--bam-programs=", args)
bam_programs <- if (any(opt)) strsplit(sub("^--bam-programs=", "", args[opt][1]), ",")[[1]] else character(0)
args <- args[!opt]
if (length(args) < 4 || length(args) > 5) {
  stop("Usage: detect_platform.R <sample_id> <seqkit_stats.tsv> <output.tsv> <reads_1> [reads_2] [--bam-programs=...]")
}
sample_id <- args[1]
stats_tsv <- args[2]
out_tsv   <- args[3]
reads     <- args[4:length(args)]


# Limits -----------------------------------------------------------------------

SHORT_MAX   <- 600    # longest read <= this: short reads. Illumina seen: 75-151 (2x300 MiSeq allowed)
LONG_MIN    <- 1000   # longest read >= this: long reads. ONT seen: 3,523-773,417
NARROW_MAX  <- 0.25   # spread <= this: uniform lengths. Illumina seen: 0.00-0.08; PacBio HiFi 0.16
WIDE_MIN    <- 0.5    # spread >= this: variable lengths. ONT seen: 0.82-1.98
BINNED_MAX  <- 8      # <= this many distinct quality values: binned (NovaSeq/NextSeq: 4-6)
MIN_READS   <- 100    # fewer reads than this: too few to trust the length shape
# Library structure (D78, provisional). Biological reads are longer than 30
# bp; index reads are 6-12 bp, 10x cell barcode + UMI reads 26-28 bp.
INDEX_MAX      <- 30
BIO_MIN        <- 50     # the other mate of a barcode or index read
# Bisulfite / EM-seq convert unmethylated C to T, so one strand's reads keep
# little C (or, for PBAT, little G) while the other base is unchanged. Two
# conditions: the depleted base under 10%, AND under half the other one.
# Checked on real reads: rice WGBS R1 3.3% C, Panax WGBS R1 6.2% C (plants
# keep more C: non-CpG methylation), against 192 ordinary real files whose
# lowest C or G share is 13.4% (C and G within a factor of 2 in every one).
BISULFITE_MAX   <- 0.10
BISULFITE_RATIO <- 0.50
INTERLEAVED_MIN <- 0.9   # share of probed record pairs with the same read name
PROBE_N     <- 2000   # reads read from the start of each file for headers and quality
HEADER_FRAC <- 0.9    # share of probed headers that must match a format to count as evidence

CONF_LEVELS <- c("uncertain", "low", "medium", "high")
cap_conf <- function(conf, max_level) {
  CONF_LEVELS[min(match(conf, CONF_LEVELS), match(max_level, CONF_LEVELS))]
}


# Helpers ----------------------------------------------------------------------

# WHAT IT DOES  Reads the first n records of a FASTQ (gzipped or not) and
#               returns their headers, read names and quality strings.
probe_fastq <- function(path, n) {
  con <- gzfile(path, "r")
  on.exit(close(con))
  x <- readLines(con, n = 4 * n, warn = FALSE)
  k <- length(x) %/% 4
  if (k == 0) return(list(headers = character(0), names = character(0), seq = character(0), qual = character(0)))
  i <- seq_len(k)
  headers <- x[4 * i - 3]
  # Read name = first word of the header, without a trailing /1 or /2, so
  # mates can be compared. Casava 1.8 puts the mate number after the space,
  # so it is already outside the first word.
  first <- sub("^@", "", sub(" .*$", "", headers))
  list(headers = headers, names = sub("/[12]$", "", first), seq = x[4 * i - 2], qual = x[4 * i])
}

# Read-name formats written by each instrument's software (D40). Each is a
# fixed, structured convention, so it does not turn up by chance.
#   pacbio   movie name + well:  m84093_240509_111803_s2/63964829/ccs (HiFi)
#                                m161120_144858_42197_c1010..._s1_p0/9807 (CLR)
#   dnbseq   MGI flowcell + Lane/Column/Row:  V350125880L1C002R03800156256
#   element  instrument ID "AV" + 6 digits:  AV224102:run:flowcell:1:10102:0388:0011
#            Element deliberately mimics the Illumina format, so it is tested
#            FIRST: otherwise Element reads could count as Illumina evidence.
#   illumina instrument:run:flowcell:lane:tile:x:y, or the older
#            instrument:lane:tile:x:y#index
#   ont      a read UUID, or MinKNOW's runid= tag
HEADER_PATTERNS <- c(
  pacbio   = "(^@| )m[0-9]+[eU]?_[0-9]{6}_[0-9]+[^/ ]*/[0-9]+(/| |$)",
  dnbseq   = "(^@| )[A-Z]{1,2}[0-9]{8,}L[0-9]C[0-9]{3}R[0-9]{3}[0-9]+",
  element  = "(^@| )AV[0-9]{6}:",
  illumina = "(^@| )[A-Za-z0-9_-]+:[0-9]+:[A-Za-z0-9_-]+:[0-9]+:[0-9]+:[0-9]+:[0-9]+|(^@| )[A-Za-z0-9_-]+:[0-9]+:[0-9]+:[0-9]+:[0-9]+(#[A-Za-z0-9+]*)?(/[12])?( |$)",
  ont      = "[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}|runid="
)

# WHAT IT DOES  Which platform the read names point to: one of
#               names(HEADER_PATTERNS), or "neutral" (unrecognised or mixed).
#               The first pattern that matches HEADER_FRAC of reads wins.
# NOTE          A trailing /1 or /2 is deliberately ignored: SRA adds it to
#               nanopore reads too (found in the survey), so it is no
#               evidence of Illumina mates.
header_evidence <- function(headers) {
  if (!length(headers)) return("neutral")
  for (p in names(HEADER_PATTERNS)) {
    if (mean(grepl(HEADER_PATTERNS[[p]], headers)) >= HEADER_FRAC) return(p)
  }
  "neutral"
}

# Library adapters: synthetic sequences each platform's kits attach to the
# DNA (D40). They survive in reads when the DNA fragment is shorter than the
# read, or when nobody trimmed them. 20 bases, searched on both strands; an
# exact chance match is about 1 in 10^12 per position.
# CAUTION  Illumina adapters identify the LIBRARY, not the instrument:
#          Element (Adept) and MGI (App-A) can sequence converted Illumina
#          libraries. So Illumina adapters are reported but never raise
#          confidence in an Illumina call (conservative rule, D1).
# Sources: Illumina Adapter Sequences document (TruSeq, Nextera); fastp
# documentation (MGI); Porechop adapter list (ONT); PacBio SMRTbell.
ADAPTERS <- list(
  illumina_library = c(truseq_r1 = "AGATCGGAAGAGCACACGTC", truseq_r2 = "AGATCGGAAGAGCGTCGTGT",
                       nextera = "CTGTCTCTTATACACATCTC"),
  dnbseq           = c(mgi_r1 = "AAGTCGGAGGCCAAGCGGTC", mgi_r2 = "AAGTCGGATCGTAGCCATGT"),
  ont              = c(ligation = "AATGTACTTCGTTCAGTTAC", rapid = "GTTTTCGCATTTATCGTGAA"),
  pacbio           = c(smrtbell = "ATCTCTCTCAACAACAACAA")
)
ADAPTER_MIN_FRAC  <- 0.005   # share of probed reads that must carry the adapter...
ADAPTER_MIN_READS <- 10      # ...and at least this many reads
ADAPTER_END_BP    <- 300     # long reads: adapters sit at the ends, so search only there

revcomp <- function(s) intToUtf8(rev(utf8ToInt(chartr("ACGT", "TGCA", s))))

# WHAT IT DOES  For each adapter group, the share of reads carrying it.
#               Returns a named numeric vector (fraction of probed reads).
adapter_scan <- function(seq) {
  if (!length(seq)) return(setNames(numeric(length(ADAPTERS)), names(ADAPTERS)))
  long <- nchar(seq) > 2 * ADAPTER_END_BP
  seq[long] <- paste0(substr(seq[long], 1, ADAPTER_END_BP), "N",
                      substring(seq[long], nchar(seq[long]) - ADAPTER_END_BP + 1))
  vapply(ADAPTERS, function(group) {
    hit <- Reduce(`|`, lapply(group, function(a) grepl(a, seq, fixed = TRUE) | grepl(revcomp(a), seq, fixed = TRUE)))
    mean(hit)
  }, 0)
}
adapter_platforms <- function(frac, n_probed) {
  names(frac)[frac >= ADAPTER_MIN_FRAC & frac * n_probed >= ADAPTER_MIN_READS]
}

# Platforms seq-triage can name but does not run QC on, and the read type
# each produces (used to check a marker agrees with the read lengths).
NOT_HANDLED <- c(pacbio = "long", dnbseq = "short", element = "short")

# Programs (BAM @PG PN, lower case) that only run on one platform's data, as
# part of the instrument's own processing. Checked against the PacBio test
# BAM (PN:ccs) and the dorado, PacBio and Picard documentation. Aligners and
# samtools (any platform) are deliberately absent: they say nothing.
BAM_SOFTWARE <- c(ccs = "pacbio", lima = "pacbio", baz2bam = "pacbio", bam2bam = "pacbio",
                  skera = "pacbio", isoseq = "pacbio", isoseq3 = "pacbio",
                  dorado = "ont", guppy = "ont", guppy_basecaller = "ont", basecall_client = "ont",
                  illuminabasecallstosam = "illumina")
DISPLAY <- c(illumina = "Illumina", ont = "ONT", pacbio = "PacBio", dnbseq = "MGI/DNBSEQ", element = "Element AVITI")

# WHAT IT DOES  What the set of quality values says: "illumina" (binned, as
#               NovaSeq/NextSeq write), "excluded" (one value only: a
#               placeholder, D26) or "neutral".
quality_evidence <- function(qual) {
  if (!length(qual)) return(list(evidence = "neutral", n_distinct = 0L))
  n <- length(unique(utf8ToInt(paste(qual, collapse = ""))))
  ev <- if (n == 1) "excluded" else if (n <= BINNED_MAX) "illumina" else "neutral"
  list(evidence = ev, n_distinct = n)
}

# WHAT IT DOES  TRUE when the reads are written with U instead of T, as
#               direct-RNA basecalls often are. Not platform evidence: it
#               tells the report the molecule, to compare with what was
#               declared. cDNA (any RNA-seq on Illumina) is written with T,
#               so its absence proves nothing.
written_as_rna <- function(seq) {
  if (!length(seq)) return(FALSE)
  mean(grepl("U", seq, fixed = TRUE) & !grepl("T", seq, fixed = TRUE)) >= HEADER_FRAC
}

# WHAT IT DOES  The core call for one file, from seqkit's longest read and
#               length quartiles. Returns call, max length and spread.
core_call <- function(row) {
  max_len <- as.numeric(row[["max_len"]])
  med     <- as.numeric(row[["Q2"]])
  spread  <- if (med > 0) (as.numeric(row[["Q3"]]) - as.numeric(row[["Q1"]])) / med else NA_real_
  call <- if (is.na(spread))                               "uncertain"
          else if (max_len <= SHORT_MAX && spread <= NARROW_MAX) "illumina"
          else if (max_len >= LONG_MIN  && spread >= WIDE_MIN)   "ont"
          else if (max_len >= LONG_MIN  && spread <= NARROW_MAX) "long_narrow"
          else                                                   "uncertain"
  list(call = call, max_len = max_len, spread = spread)
}


# Read inputs ------------------------------------------------------------------

# seqkit stats --tabular --all: one row per file. Match rows to the read
# files by name rather than trusting row order.
# seqkit can fail SILENTLY: on a truncated gzip it prints "unexpected EOF"
# but exits 0 and writes a table with no row for that file (found in
# validation round 1). So a missing row means seqkit could not read the file.
stats <- read.delim(stats_tsv, check.names = FALSE, stringsAsFactors = FALSE, colClasses = c(file = "character"))
unread <- reads[!basename(reads) %in% basename(stats$file)]
if (length(unread)) {
  stop(sprintf(paste("Sample %s: seqkit could not read %s.",
                     "The file is probably truncated or corrupt (for example, an interrupted download).",
                     "Check it with: gzip -t <file>"),
               sample_id, paste(basename(unread), collapse = " and ")), call. = FALSE)
}
stats <- stats[match(basename(reads), basename(stats$file)), , drop = FALSE]

n_files <- length(reads)
n_reads <- as.numeric(stats$num_seqs)
flags   <- character(0)
why     <- character(0)

write_result <- function(platform, layout, confidence, evidence = list()) {
  out <- data.frame(
    sample            = sample_id,
    n_files           = n_files,
    n_reads           = format(n_reads[1], scientific = FALSE),   # fastq_1 (pairs, if paired)
    detected_platform = platform,
    detected_layout   = layout,
    confidence        = confidence,
    detection_flags   = paste(unique(flags), collapse = ","),
    reasons           = paste(why, collapse = " "),
    max_len           = paste(evidence$max_len, collapse = ";"),
    length_spread     = paste(evidence$spread,  collapse = ";"),
    header_evidence   = paste(evidence$header,  collapse = ";"),
    quality_values    = paste(evidence$n_q,     collapse = ";"),
    u_bases           = paste(evidence$u,       collapse = ";"),
    adapters          = paste(evidence$adapters, collapse = ";"),
    start_share       = paste(evidence$start_share, collapse = ";"),
    probe_reads       = paste(evidence$probe_reads, collapse = ";"),
    stringsAsFactors  = FALSE
  )
  write.table(out, out_tsv, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  quit(save = "no", status = 0)
}

# An empty file has no lengths to judge. The workflow routes it to 'empty'
# (FAIL); here it is simply uncertain.
if (any(n_reads == 0)) {
  why <- "At least one file contains no reads, so the platform cannot be judged."
  write_result("uncertain", "uncertain", "uncertain")
}


# 1. Evidence per file -----------------------------------------------------------

core   <- lapply(seq_len(n_files), function(i) core_call(stats[i, ]))
calls  <- vapply(core, `[[`, "", "call")
max_ln <- vapply(core, `[[`, 0, "max_len")
spread <- round(vapply(core, `[[`, 0, "spread"), 3)

probes  <- lapply(reads, probe_fastq, n = PROBE_N)
hdr_ev  <- vapply(probes, function(p) header_evidence(p$headers), "")
q_ev    <- lapply(probes, function(p) quality_evidence(p$qual))
q_evid  <- vapply(q_ev, `[[`, "", "evidence")
q_n     <- vapply(q_ev, `[[`, 0L, "n_distinct")
u_bases <- vapply(probes, function(p) written_as_rna(p$seq), FALSE)
ad_frac <- lapply(probes, function(p) adapter_scan(p$seq))
ad_ev   <- unique(unlist(lapply(seq_len(n_files), function(i) adapter_platforms(ad_frac[[i]], length(probes[[i]]$seq)))))
# Read-start concentration (D65): the share of probed reads (first file)
# whose first 20 bases are one of the 10 most common 20-base starts. Amplicon
# reads begin at primers, so a few starts dominate; reads from fragmented DNA
# start almost anywhere. Used for the sample-type call in the report, and
# only for short reads: nanopore reads start with adapter sequence.
start_share <- function(seq) {
  s <- seq[nchar(seq) >= 20]
  if (!length(s)) return(NA_real_)
  st <- sort(table(substr(s, 1, 20)), decreasing = TRUE)
  round(sum(head(st, 10)) / length(s), 4)
}
evidence <- list(max_len = max_ln, spread = spread, header = hdr_ev, n_q = q_n, u = u_bases,
                 start_share = start_share(probes[[1]]$seq), probe_reads = length(probes[[1]]$seq),
                 adapters = vapply(ad_frac, function(f) paste(sprintf("%s=%.3f", names(f), f), collapse = ","), ""))
# Library structure (D78) -----------------------------------------------------
base_frac <- function(seq, b) {
  n <- sum(nchar(seq))
  if (n == 0) return(NA_real_)
  sum(nchar(gsub(sprintf("[^%s]", b), "", seq))) / n
}
lib_short <- max_ln <= INDEX_MAX
if (all(lib_short)) {
  flags <- c(flags, "INDEX_READS")
  why   <- sprintf("Every read is %s bp or shorter (%s): these are index or barcode reads, not reads of the sample's DNA. Give the sample's read files (R1, R2) instead; no QC was run.",
                   INDEX_MAX, paste(sprintf("longest %s bp", max_ln), collapse = "; "))
  write_result("uncertain", if (n_files == 1) "single" else "uncertain", "uncertain", evidence)
}
if (n_files == 2 && lib_short[1] && max_ln[2] >= BIO_MIN) {
  # 10x Genomics-style single-cell: R1 = cell barcode + UMI, R2 = cDNA.
  flags <- c(flags, "BARCODE_READ")
  why   <- c(why, sprintf("R1 is only %s bp and R2 %s bp: R1 is a cell-barcode/UMI read, as 10x Genomics-style single-cell libraries write, so this looks like single-cell data. QC is run on R2, the read of the sample's molecules.",
                          max_ln[1], max_ln[2]))
  evidence$start_share <- start_share(probes[[2]]$seq)
  evidence$probe_reads <- length(probes[[2]]$seq)
} else if (n_files == 2 && lib_short[2] && max_ln[1] >= BIO_MIN) {
  flags <- c(flags, "INDEX_READ_AS_MATE")
  why   <- c(why, sprintf("fastq_2 is only %s bp: it looks like an index (I1) read given as the second mate, not R2. QC is run on fastq_1 alone, as single-end; give R2 as fastq_2 if there is one.",
                          max_ln[2]))
}
bio <- if ("BARCODE_READ" %in% flags) 2 else 1
c_frac <- base_frac(probes[[bio]]$seq, "C"); g_frac <- base_frac(probes[[bio]]$seq, "G")
if (!is.na(c_frac) && min(c_frac, g_frac) < BISULFITE_MAX &&
    min(c_frac, g_frac) < BISULFITE_RATIO * max(c_frac, g_frac)) {
  flags <- c(flags, "BISULFITE")
  why   <- c(why, sprintf("Only %.1f%% of bases are C and %.1f%% G in %s: one strand has lost its C, as bisulfite or EM-seq conversion does, so this looks like a methylation library. Converted reads cannot be classified reliably, so contamination is not judged.",
                          100 * c_frac, 100 * g_frac, if (bio == 2) "R2" else "the reads"))
}
interleaved <- FALSE
if (n_files == 1 && max_ln[1] <= SHORT_MAX) {
  nm <- probes[[1]]$names
  k  <- 2 * (length(nm) %/% 2)
  if (k >= 2 && n_reads[1] %% 2 == 0 &&
      mean(nm[seq(1, k, 2)] == nm[seq(2, k, 2)]) >= INTERLEAVED_MIN) {
    interleaved <- TRUE
    flags <- c(flags, "INTERLEAVED")
    why   <- c(why, sprintf("One file, but each read name appears twice in a row (first %s reads): both mates are in one interleaved file. It is split into R1 and R2 before QC.",
                            format(k, big.mark = ",")))
  }
}

if (any(u_bases)) {
  flags <- c(flags, "U_BASES")
  why   <- c(why, "Reads are written with U instead of T, as direct-RNA basecalls are.")
}

describe <- function(i) sprintf("longest read %s bp, length spread %.2f", format(max_ln[i], big.mark = ","), spread[i])
pct <- function(group) sprintf("%.1f%%", 100 * max(vapply(ad_frac, `[[`, 0, group)))
read_kind <- if (all(max_ln >= LONG_MIN)) "long" else if (all(max_ln <= SHORT_MAX)) "short" else "mixed"

if (length(unique(calls)) > 1) {
  flags <- c(flags, "FILES_DISAGREE")
  why   <- sprintf("The two files look like different data (%s; %s), so neither platform nor layout can be trusted.",
                   describe(1), describe(2))
  write_result("uncertain", "uncertain", "uncertain", evidence)
}


# 2. A platform seq-triage does not handle? (D40) --------------------------------
# Instrument markers (read names, or MGI/PacBio adapters) name the platform
# when they agree with the read lengths. Illumina adapters are not markers
# here: they identify the library kit, not the instrument.

# BAM header software (D55): which platform's own software wrote the file.
sw_hit <- unique(BAM_SOFTWARE[intersect(tolower(bam_programs), names(BAM_SOFTWARE))])
sw_names <- intersect(tolower(bam_programs), names(BAM_SOFTWARE))
if (length(sw_hit) > 1) {
  flags <- c(flags, "SIGNALS_DISAGREE")
  why   <- c(why, sprintf("The BAM header names software of more than one platform (%s), so it is not used as evidence.",
                          paste(sw_names, collapse = ", ")))
  sw_hit <- character(0)
}
sw_text <- if (length(sw_hit)) sprintf("the BAM header (written by %s)", paste(sw_names, collapse = ", ")) else ""

markers <- unique(c(intersect(hdr_ev, names(NOT_HANDLED)), intersect(ad_ev, names(NOT_HANDLED)),
                    intersect(sw_hit, names(NOT_HANDLED))))
platform <- NULL
if (length(markers) > 1) {
  flags <- c(flags, "SIGNALS_DISAGREE")
  why   <- c(sprintf("Markers of more than one platform were found (%s), so the platform cannot be named (%s).",
                     paste(DISPLAY[markers], collapse = ", "), describe(1)), why)
  write_result("uncertain", if (n_files == 1) "single" else "uncertain", "uncertain", evidence)
} else if (length(markers) == 1) {
  how <- c(if (markers %in% hdr_ev) "read names" else NULL,
           if (markers %in% ad_ev) sprintf("%s adapters in %s of reads", DISPLAY[[markers]], pct(markers)) else NULL,
           if (markers %in% sw_hit) sw_text else NULL)
  if (NOT_HANDLED[[markers]] != read_kind) {
    flags <- c(flags, "SIGNALS_DISAGREE")
    why   <- c(sprintf("%s identify %s, but %s reads do not fit it (%s).", paste(how, collapse = " and "),
                       DISPLAY[[markers]], if (read_kind == "long") "long" else "short or mixed-length", describe(1)), why)
    write_result("uncertain", if (n_files == 1) "single" else "uncertain", "uncertain", evidence)
  }
  platform   <- markers
  confidence <- "high"
  # PacBio HiFi (D64) is the one named platform seq-triage QCs. HiFi needs a
  # HiFi-specific marker on top of the PacBio one: read names ending /ccs
  # (also /ccs/fwd, /ccs/rev), or a BAM written by ccs. Other PacBio data
  # (CLR, subreads) keeps PLATFORM_NOT_HANDLED: its quality values are not
  # comparable, and are often placeholders.
  hifi_names <- any(vapply(probes, function(p) any(grepl("/ccs(/|[[:space:]]|$)", p$headers)), FALSE))
  hifi_bam   <- "ccs" %in% tolower(bam_programs)
  if (markers == "pacbio" && (hifi_names || hifi_bam)) {
    hifi_how <- c(if (hifi_names) "read names ending /ccs" else NULL, if (hifi_bam) "a BAM written by ccs" else NULL)
    why <- c(sprintf("%s identify PacBio, %s marks the reads as HiFi, and the read lengths fit (%s).",
                     paste(how, collapse = " and "), paste(hifi_how, collapse = " and "), describe(1)), why)
  } else if (markers %in% c("dnbseq", "element")) {
    # D69: MGI/DNBSEQ and Element are QC'd like Illumina (short, fixed-length
    # reads), each with its own rows in the metric guide.
    why <- c(sprintf("%s identify %s, and the read lengths fit (%s).",
                     paste(how, collapse = " and "), DISPLAY[[markers]], describe(1)), why)
  } else {
    flags <- c(flags, "PLATFORM_NOT_HANDLED")
    why   <- c(sprintf("%s identify %s, and the read lengths fit (%s). seq-triage does not run QC on %s data%s.",
                       paste(how, collapse = " and "), DISPLAY[[markers]], describe(1), DISPLAY[[markers]],
                       if (markers == "pacbio") " other than HiFi (no read names ending /ccs, no ccs in a BAM header)" else ""), why)
  }
}


# 3. Core call: Illumina or ONT ---------------------------------------------------

if (is.null(platform)) {
  call <- calls[1]
  if (call == "long_narrow") {
    flags <- c(flags, "LONG_READS_NARROW_SPREAD")
    why   <- c(sprintf("Reads are long but unusually uniform in length (%s). Nanopore lengths vary widely; this pattern fits PacBio HiFi or size-selected data, which seq-triage does not handle.",
                       describe(1)), why)
    write_result("uncertain", if (n_files == 1) "single" else "uncertain", "uncertain", evidence)
  }
  if (call == "uncertain") {
    why <- c(sprintf("Read lengths fit neither platform clearly (%s). Illumina needs longest read <= %d bp and spread <= %.2f; ONT needs longest read >= %d bp and spread >= %.2f.",
                     describe(1), SHORT_MAX, NARROW_MAX, LONG_MIN, WIDE_MIN), why)
    write_result("uncertain", if (n_files == 1) "single" else "uncertain", "uncertain", evidence)
  }
  platform <- call
  why <- c(sprintf("Read lengths indicate %s (%s).", DISPLAY[[call]], describe(1)), why)

  # Supporting evidence: raises or lowers confidence, never makes the call.
  other <- if (call == "illumina") "ont" else "illumina"
  support <- character(0)
  against <- character(0)
  if (any(hdr_ev == call))  support <- c(support, sprintf("%s-format read names", DISPLAY[[call]]))
  if (any(hdr_ev == other)) against <- c(against, sprintf("%s-format read names", DISPLAY[[other]]))
  if (any(q_evid == call))  support <- c(support, sprintf("binned quality values (%d distinct), as NovaSeq/NextSeq write", max(q_n)))
  if (any(q_evid == other)) against <- c(against, sprintf("binned quality values (%d distinct), typical of Illumina", max(q_n)))
  if (any(sw_hit == call))  support <- c(support, sw_text)
  if (any(sw_hit == other)) against <- c(against, sw_text)
  if (call == "ont" && "ont" %in% ad_ev)       support <- c(support, sprintf("ONT adapters in %s of reads", pct("ont")))
  if (call == "illumina" && "ont" %in% ad_ev)  against <- c(against, sprintf("ONT adapters in %s of reads", pct("ont")))
  if (call == "ont" && "illumina_library" %in% ad_ev) against <- c(against, sprintf("Illumina library adapters in %s of reads", pct("illumina_library")))
  if (call == "illumina" && "illumina_library" %in% ad_ev) {
    why <- c(why, sprintf("Illumina library adapters are in %s of reads. They fit, but do not confirm the instrument: Element and MGI instruments can sequence converted Illumina libraries.",
                          pct("illumina_library")))
  }
  if (any(q_evid == "excluded")) {
    flags <- c(flags, "QUALITY_UNINFORMATIVE")
    why   <- c(why, "Every base has the same quality score, so quality is a placeholder and was ignored.")
  }

  if (length(against)) {
    confidence <- "low"
    flags <- c(flags, "SIGNALS_DISAGREE")
    why   <- c(why, sprintf("But %s point the other way.", paste(against, collapse = " and ")))
  } else if (call == "illumina" && length(support) == 1 && any(q_evid == call)) {
    # Binned quality alone does not name the instrument (D58): validation
    # round 2 found SRR14318705, recorded in ENA as BGISEQ-500 with read
    # names stripped, whose quality uses exactly Illumina's 7-value binning;
    # it was called Illumina at high confidence on that alone.
    confidence <- "medium"
    why   <- c(why, sprintf(paste("Binned quality values (%d distinct) fit Illumina, but they are the only confirmation and",
                                  "do not name the instrument (an archived MGI run was found with the same binning), so",
                                  "without Illumina read names confidence stays medium."), max(q_n)))
  } else if (length(support)) {
    confidence <- "high"
    why   <- c(why, sprintf("Confirmed by %s.", paste(support, collapse = " and ")))
  } else {
    confidence <- "medium"
    why   <- c(why, "No read-name, adapter or quality evidence confirms it, so confidence is medium.")
  }
} else if (any(q_evid == "excluded")) {
  flags <- c(flags, "QUALITY_UNINFORMATIVE")
}

if (any(n_reads < MIN_READS)) {
  confidence <- cap_conf(confidence, "low")
  flags <- c(flags, "FEW_READS")
  why   <- c(why, sprintf("Only %s reads, too few to trust the length pattern; confidence capped at low.",
                          format(min(n_reads), big.mark = ",")))
}


# 4. Layout --------------------------------------------------------------------

if (n_files == 1) {
  layout <- if (interleaved && read_kind == "short") "interleaved" else "single"
} else if (read_kind == "long") {
  layout <- "uncertain"
  flags  <- c(flags, "TWO_FILES_LONG_READS")
  why    <- c(why, "Two files were given, but long reads are not paired. Check the samplesheet; no QC is run on this sample.")
} else if (n_reads[1] != n_reads[2]) {
  layout     <- "uncertain"
  confidence <- cap_conf(confidence, "low")
  flags      <- c(flags, "MATE_COUNT_MISMATCH")
  why        <- c(why, sprintf("The two files have different read counts (%s and %s), so they are not a clean pair.",
                               format(n_reads[1], big.mark = ","), format(n_reads[2], big.mark = ",")))
} else {
  n1 <- probes[[1]]$names
  n2 <- probes[[2]]$names
  k  <- min(length(n1), length(n2))
  if (k == 0 || any(n1[seq_len(k)] != n2[seq_len(k)])) {
    layout     <- "uncertain"
    confidence <- cap_conf(confidence, "low")
    flags      <- c(flags, "MATE_NAMES_DIFFER")
    why        <- c(why, "Read names in the two files do not match, so the files may not be mates of each other.")
  } else {
    layout <- "paired"
    why    <- c(why, sprintf("Paired: equal read counts and matching read names in the first %s reads.",
                             format(k, big.mark = ",")))
  }
}

write_result(platform, layout, confidence, evidence)
