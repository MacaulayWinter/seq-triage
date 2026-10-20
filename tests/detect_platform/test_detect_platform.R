#!/usr/bin/env Rscript
# Test for bin/detect_platform.R
#
# Builds small SYNTHETIC samples in a temporary folder (random bases, chosen
# read lengths, headers and quality strings), writes the seqkit-style stats
# for each, runs the detection script, and compares platform, layout,
# confidence and flags with the expected answer. One case per rule. No real
# data is used or stored.
#
# Run from the repository root:
#   Rscript tests/detect_platform/test_detect_platform.R
# Exits non-zero if any case differs, so CI can use it.

set.seed(1)
tmp <- tempfile("detect_test_")
dir.create(tmp)

# Header styles seen in the test-data survey
uuid <- function(n) vapply(seq_len(n), function(i) {
  h <- sample(c(0:9, letters[1:6]), 32, replace = TRUE)
  paste(c(h[1:8], "-", h[9:12], "-", h[13:16], "-", h[17:20], "-", h[21:32]), collapse = "")
}, "")
HEADERS <- list(
  casava    = function(n, mate) sprintf("@A00611:237:H5335DRX5:2:2101:%d:1000 %d:N:0:ACGT", seq_len(n), mate),
  old_hiseq = function(n, mate) sprintf("@HS11_6401:7:2206:%d:167822#14/%d", seq_len(n), mate),
  ont       = function(n, mate) sprintf("@%s runid=811ac9822692f2545ea99b1a0b2f9336ecf68626", uuid(n)),
  ont_sra   = function(n, mate) sprintf("@SRR1.%d %s/1", seq_len(n), uuid(n)),   # SRA adds /1 to nanopore reads
  plain     = function(n, mate) sprintf("@read%d/%d", seq_len(n), mate),
  renamed   = function(n, mate) sprintf("@other%d/%d", seq_len(n), mate)
)
# Quality alphabets: binned (NovaSeq, 4 values), full range, placeholder (1 value)
QUAL <- list(binned = c("F", ":", ",", "#"), full = intToUtf8(35:74, multiple = TRUE), placeholder = "I")
LENGTHS <- list(
  illumina    = function(n) sample(c(rep(151, 9), 150), n, replace = TRUE),
  ont         = function(n) pmax(100, round(rlnorm(n, log(4000), 1))),
  pacbio      = function(n) round(rnorm(n, 18000, 1500)),
  middle      = function(n) sample(300:800, n, replace = TRUE)          # fits neither platform
)

# WHAT IT DOES  Writes one FASTQ and returns its seqkit-style stats row.
make_fastq <- function(path, n, lengths, header, qual, mate = 1, alphabet = "dna") {
  len <- LENGTHS[[lengths]](n)
  hdr <- HEADERS[[header]](n, mate)
  bases <- if (alphabet == "rna") c("A", "C", "G", "U") else c("A", "C", "G", "T")
  seq <- vapply(len, function(l) paste(sample(bases, l, TRUE), collapse = ""), "")
  q   <- vapply(len, function(l) paste(sample(QUAL[[qual]], l, TRUE), collapse = ""), "")
  con <- gzfile(path, "w")
  if (n) writeLines(c(rbind(hdr, seq, "+", q)), con)
  close(con)
  qs  <- if (n) quantile(len, c(.25, .5, .75), names = FALSE) else c(0, 0, 0)
  data.frame(file = basename(path), num_seqs = n, max_len = if (n) max(len) else 0,
             Q1 = qs[1], Q2 = qs[2], Q3 = qs[3])
}

# Each case: files as "n:lengths:header:qual[:rna]", and the expected result
cases <- read.table(header = TRUE, sep = "|", strip.white = TRUE, stringsAsFactors = FALSE, text = "
id                   | file1                         | file2                          | platform  | layout    | confidence | flags
ILLUMINA_PAIRED      | 500:illumina:casava:binned    | 500:illumina:casava:binned     | illumina  | paired    | high       |
ILLUMINA_OLD_HEADER  | 500:illumina:old_hiseq:full   |                                | illumina  | single    | high       |
ILLUMINA_NO_EVIDENCE | 500:illumina:plain:full       |                                | illumina  | single    | medium     |
ONT_HEADER           | 500:ont:ont:full              |                                | ont       | single    | high       |
ONT_SRA_SLASH1       | 500:ont:ont_sra:full          |                                | ont       | single    | high       |
ONT_PLACEHOLDER_Q    | 500:ont:plain:placeholder     |                                | ont       | single    | medium     | QUALITY_UNINFORMATIVE
PACBIO_LIKE          | 200:pacbio:plain:full         |                                | uncertain | single    | uncertain  | LONG_READS_NARROW_SPREAD
MIDDLE_LENGTHS       | 500:middle:plain:full         |                                | uncertain | single    | uncertain  |
SHORT_WITH_ONT_IDS   | 500:illumina:ont:full         |                                | illumina  | single    | low        | SIGNALS_DISAGREE
LONG_WITH_BINNED_Q   | 500:ont:ont:binned            |                                | ont       | single    | low        | SIGNALS_DISAGREE
FEW_READS            | 50:illumina:casava:binned     |                                | illumina  | single    | low        | FEW_READS
MATE_COUNTS_DIFFER   | 500:illumina:casava:binned    | 450:illumina:casava:binned     | illumina  | uncertain | low        | MATE_COUNT_MISMATCH
MATE_NAMES_DIFFER    | 500:illumina:plain:binned     | 500:illumina:renamed:binned    | illumina  | uncertain | low        | MATE_NAMES_DIFFER
ONT_TWO_FILES        | 500:ont:ont:full              | 500:ont:ont:full               | ont       | uncertain | high       | TWO_FILES_LONG_READS
FILES_DISAGREE       | 500:illumina:casava:binned    | 500:ont:ont:full               | uncertain | uncertain | uncertain  | FILES_DISAGREE
EMPTY                | 0:illumina:casava:binned      |                                | uncertain | uncertain | uncertain  |
ONT_DIRECT_RNA       | 500:ont:ont:full:rna          |                                | ont       | single    | high       | U_BASES
")

run_case <- function(i) {
  cs    <- cases[i, ]
  specs <- Filter(nzchar, c(cs$file1, cs$file2))
  rows  <- list(); paths <- character(0)
  for (m in seq_along(specs)) {
    p <- strsplit(specs[m], ":")[[1]]
    path  <- file.path(tmp, sprintf("%s_R%d.fastq.gz", cs$id, m))
    rows[[m]] <- make_fastq(path, as.integer(p[1]), p[2], p[3], p[4], mate = m,
                            alphabet = if (length(p) > 4) p[5] else "dna")
    paths <- c(paths, path)
  }
  stats <- file.path(tmp, paste0(cs$id, ".stats.tsv"))
  write.table(do.call(rbind, rows), stats, sep = "\t", quote = FALSE, row.names = FALSE)
  out <- file.path(tmp, paste0(cs$id, ".detection.tsv"))
  status <- system2("Rscript", c("--vanilla", "bin/detect_platform.R", cs$id, stats, out, paths))
  if (status != 0) return(c(got = sprintf("script exited %d", status), expected = ""))
  d <- read.delim(out, colClasses = "character", na.strings = "")
  sorted <- function(f) if (is.na(f) || f == "") "" else paste(sort(strsplit(f, ",")[[1]]), collapse = ",")
  c(got      = paste(d$detected_platform, d$detected_layout, d$confidence, sorted(d$detection_flags)),
    expected = paste(cs$platform, cs$layout, cs$confidence, sorted(cs$flags)))
}

res <- data.frame(case = cases$id, t(vapply(seq_len(nrow(cases)), run_case, c(got = "", expected = ""))),
                  stringsAsFactors = FALSE)
res$result <- ifelse(res$got == res$expected, "PASS", "FAIL")
print(res, row.names = FALSE, right = FALSE)

n_fail <- sum(res$result == "FAIL")
cat(sprintf("\n%d / %d cases as expected\n", nrow(res) - n_fail, nrow(res)))
unlink(tmp, recursive = TRUE)
if (n_fail > 0) quit(status = 1)
