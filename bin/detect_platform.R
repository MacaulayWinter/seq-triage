#!/usr/bin/env Rscript
#*******************************************************************************
# detect_platform.R: infer sequencing platform and layout from the reads
#
# WHAT IT DOES  Calls each sample Illumina, ONT or uncertain, with a confidence
#               level, a layout, warning flags and a plain-language reason.
# HOW           1. CORE CALL from read-length shape (seqkit stats, whole
#                  file). The longest read and the length spread must agree;
#                  one alone is never enough (CLAUDE.md D3).
#               2. SUPPORTING EVIDENCE from the first PROBE_N reads: header
#                  format and the number of distinct quality values. These
#                  raise or lower confidence but never make the call.
#               3. LAYOUT from the number of files, their read counts and
#                  their read names.
# RULE          When the evidence is weak or disagrees, say so. "uncertain"
#               is a result, not an error (D1).
# CALIBRATION   Every limit below was checked against 21 real and test files
#               (CLAUDE.md "Test data survey" and D33). Observed ranges are
#               noted next to each limit.
#
# Usage:  detect_platform.R <sample_id> <seqkit_stats.tsv> <output.tsv> <reads_1> [reads_2]
#*******************************************************************************

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4 || length(args) > 5) {
  stop("Usage: detect_platform.R <sample_id> <seqkit_stats.tsv> <output.tsv> <reads_1> [reads_2]")
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

# WHAT IT DOES  Which platform the header format points to: "illumina",
#               "ont" or "neutral" (unrecognised or mixed).
# NOTE          A trailing /1 or /2 is deliberately ignored: SRA adds it to
#               nanopore reads too (found in the survey), so it is no
#               evidence of Illumina mates.
header_evidence <- function(headers) {
  if (!length(headers)) return("neutral")
  illumina_new <- grepl("(^@| )[A-Za-z0-9_-]+:[0-9]+:[A-Za-z0-9_-]+:[0-9]+:[0-9]+:[0-9]+:[0-9]+", headers)
  illumina_old <- grepl("(^@| )[A-Za-z0-9_-]+:[0-9]+:[0-9]+:[0-9]+:[0-9]+(#[A-Za-z0-9+]*)?(/[12])?( |$)", headers)
  ont_uuid     <- grepl("[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}", headers) |
                  grepl("runid=", headers, fixed = TRUE)
  if (mean(illumina_new | illumina_old) >= HEADER_FRAC) return("illumina")
  if (mean(ont_uuid) >= HEADER_FRAC) return("ont")
  "neutral"
}

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


# 1. Core call per file ----------------------------------------------------------

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
evidence <- list(max_len = max_ln, spread = spread, header = hdr_ev, n_q = q_n, u = u_bases)
if (any(u_bases)) {
  flags <- c(flags, "U_BASES")
  why   <- c(why, "Reads are written with U instead of T, as direct-RNA basecalls are.")
}

describe <- function(i) sprintf("longest read %s bp, length spread %.2f", format(max_ln[i], big.mark = ","), spread[i])

if (length(unique(calls)) > 1) {
  flags <- c(flags, "FILES_DISAGREE")
  why   <- sprintf("The two files look like different data (%s; %s), so neither platform nor layout can be trusted.",
                   describe(1), describe(2))
  write_result("uncertain", "uncertain", "uncertain", evidence)
}

call <- calls[1]
if (call == "long_narrow") {
  flags <- c(flags, "LONG_READS_NARROW_SPREAD")
  why   <- sprintf("Reads are long but unusually uniform in length (%s). Nanopore lengths vary widely; this pattern fits PacBio HiFi or size-selected data, which seq-triage does not handle.",
                   describe(1))
  write_result("uncertain", if (n_files == 1) "single" else "uncertain", "uncertain", evidence)
}
if (call == "uncertain") {
  why <- sprintf("Read lengths fit neither platform clearly (%s). Illumina needs longest read <= %d bp and spread <= %.2f; ONT needs longest read >= %d bp and spread >= %.2f.",
                 describe(1), SHORT_MAX, NARROW_MAX, LONG_MIN, WIDE_MIN)
  write_result("uncertain", if (n_files == 1) "single" else "uncertain", "uncertain", evidence)
}

why <- c(sprintf("Read lengths indicate %s (%s).", if (call == "illumina") "Illumina" else "ONT", describe(1)), why)


# 2. Supporting evidence -------------------------------------------------------

other <- if (call == "illumina") "ont" else "illumina"
platform_name <- function(p) if (p == "illumina") "Illumina" else "ONT"

support <- character(0)
against <- character(0)
if (any(hdr_ev == call))  support <- c(support, sprintf("%s-format read headers", platform_name(call)))
if (any(hdr_ev == other)) against <- c(against, sprintf("%s-format read headers", platform_name(other)))
if (any(q_evid == call))  support <- c(support, sprintf("binned quality values (%d distinct), as NovaSeq/NextSeq write", max(q_n)))
if (any(q_evid == other)) against <- c(against, sprintf("binned quality values (%d distinct), typical of Illumina", max(q_n)))
if (any(q_evid == "excluded")) {
  flags <- c(flags, "QUALITY_UNINFORMATIVE")
  why   <- c(why, "Every base has the same quality score, so quality is a placeholder and was ignored.")
}

if (length(against)) {
  confidence <- "low"
  flags <- c(flags, "SIGNALS_DISAGREE")
  why   <- c(why, sprintf("But %s point the other way.", paste(against, collapse = " and ")))
} else if (length(support)) {
  confidence <- "high"
  why   <- c(why, sprintf("Confirmed by %s.", paste(support, collapse = " and ")))
} else {
  confidence <- "medium"
  why   <- c(why, "No header or quality evidence confirms it, so confidence is medium.")
}

if (any(n_reads < MIN_READS)) {
  confidence <- cap_conf(confidence, "low")
  flags <- c(flags, "FEW_READS")
  why   <- c(why, sprintf("Only %s reads, too few to trust the length pattern; confidence capped at low.",
                          format(min(n_reads), big.mark = ",")))
}


# 3. Layout --------------------------------------------------------------------

if (n_files == 1) {
  layout <- "single"
} else if (call == "ont") {
  layout <- "uncertain"
  flags  <- c(flags, "TWO_FILES_LONG_READS")
  why    <- c(why, "Two files were given, but nanopore reads are not paired. Check the samplesheet; no QC is run on this sample.")
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

write_result(call, layout, confidence, evidence)
