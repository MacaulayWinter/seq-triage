#!/usr/bin/env Rscript
#*******************************************************************************
# make_faults.R: build deliberately damaged samples with known answers
#                (CLAUDE.md D37, layer C)
#
# WHAT IT DOES  Starts from two known-good samples (E. coli ECE107: NovaSeq
#               120k pairs and R10.4.1 12k reads, both PASS today), damages
#               them in one controlled way each, and writes the files, the
#               samplesheets and expected.tsv: the verdict and flags each
#               sample MUST get. Expected answers are written HERE, before the
#               pipeline runs, and are never edited to match a result.
# BATCHES       "main" runs together. "isolated_*" cases may stop a whole run
#               by design (D18), so each runs on its own.
# SEED          Fixed, so the same files are built every time.
#
# Usage:  Rscript validation/make_faults.R <test_data_dir> <other_illumina_R2.fastq.gz> <out_dir>
#         test_data_dir = a checkout of MacaulayWinter/seq-triage-test-data
#*******************************************************************************

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) stop("Usage: make_faults.R <test_data_dir> <other_illumina_R2.fastq.gz> <out_dir>")
script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(script_dir, "lib_fastq.R"))
td    <- args[1]
other <- args[2]
out   <- normalizePath(args[3], mustWork = FALSE)
dir.create(file.path(out, "fastq"), recursive = TRUE, showWarnings = FALSE)
set.seed(20261020)

message("Reading base samples ...")
ilm1 <- read_fq(file.path(td, "illumina/SRR36800568_120k_R1.fastq.gz"))
ilm2 <- read_fq(file.path(td, "illumina/SRR36800568_120k_R2.fastq.gz"))
ont  <- read_fq(file.path(td, "nanopore/SRR34410547_12k.fastq.gz"))
n_ilm <- length(ilm1$h)
n_ont <- length(ont$h)

rows <- list()
path <- function(id, mate = NULL) file.path(out, "fastq", paste0(id, if (!is.null(mate)) paste0("_R", mate), ".fastq.gz"))

# WHAT IT DOES  Writes one case's files and records its expected answer.
add <- function(id, r1, r2 = NULL, platform = "", sample_type = "isolate", molecule = "dna",
                verdict, flags = "", batch = "main", note, files = NULL) {
  if (is.null(files)) {
    files <- c(write_fq(r1, path(id, if (is.null(r2)) NULL else 1)),
               if (!is.null(r2)) write_fq(r2, path(id, 2)))
  }
  rows[[id]] <<- data.frame(sample = id, batch = batch, fastq_1 = files[1],
                            fastq_2 = if (length(files) > 1) files[2] else "",
                            platform = platform, sample_type = sample_type, molecule = molecule,
                            expected_verdict = verdict, expected_flags = flags, note = note)
  message(sprintf("  %-22s expect %-9s %s", id, verdict, flags))
}

# Base-level edits ---------------------------------------------------------------

# Replace a random fraction p of characters equal to `from` with `to`.
swap_chars <- function(x, from, to, p) {
  vapply(x, function(s) {
    v <- strsplit(s, "")[[1]]
    hit <- which(v == from)
    v[hit[runif(length(hit)) < p]] <- to
    paste(v, collapse = "")
  }, "", USE.NAMES = FALSE)
}
# Set a random fraction p of bases to N (quality '#', as instruments write).
add_ns <- function(fq, p) {
  for (k in seq_along(fq$s)) {
    s <- strsplit(fq$s[k], "")[[1]]; q <- strsplit(fq$q[k], "")[[1]]
    hit <- runif(length(s)) < p
    s[hit] <- "N"; q[hit] <- "#"
    fq$s[k] <- paste(s, collapse = ""); fq$q[k] <- paste(q, collapse = "")
  }
  fq
}
# Keep the first fraction f of every read (lengths scale, their spread stays).
shorten <- function(fq, f) {
  keep <- pmax(1L, round(nchar(fq$s) * f))
  list(h = fq$h, s = substr(fq$s, 1, keep), q = substr(fq$q, 1, keep))
}
# Lower every base quality by d (minimum Q1).
lower_q <- function(fq, d) {
  fq$q <- vapply(fq$q, function(q) intToUtf8(pmax(1L, phred(q) - d) + 33L), "", USE.NAMES = FALSE)
  fq
}
mean_q <- function(fq) mean(vapply(fq$q, read_mean_q, 0, USE.NAMES = FALSE))
revcomp <- function(s) intToUtf8(rev(utf8ToInt(chartr("ACGTN", "TGCAN", s))))

# Controls: the unmodified base samples -----------------------------------------
add("CTRL_ILM", ilm1, ilm2, platform = "illumina", verdict = "PASS",
    note = "unmodified ECE107 NovaSeq: the reference every Illumina case is compared with")
add("CTRL_ONT", ont, platform = "ont", verdict = "PASS",
    note = "unmodified ECE107 R10.4.1: the reference every ONT case is compared with")

# Read counts --------------------------------------------------------------------
add("ILM_READS_40K", subset_fq(ilm1, 1:20000), subset_fq(ilm2, 1:20000), platform = "illumina",
    verdict = "WARN", note = "20,000 pairs = 40,000 reads, under the 100,000 warn limit")
add("ILM_READS_400", subset_fq(ilm1, 1:200), subset_fq(ilm2, 1:200), platform = "illumina",
    verdict = "FAIL", note = "200 pairs = 400 reads, under the 1,000 fail limit")
add("ONT_READS_5K", subset_fq(ont, 1:5000), platform = "ont",
    verdict = "WARN", note = "5,000 reads, under the 10,000 warn limit")
add("ONT_READS_500", subset_fq(ont, 1:500), platform = "ont",
    verdict = "FAIL", note = "500 reads, under the 1,000 fail limit")

# Read length --------------------------------------------------------------------
# Scaling every read keeps the length SPREAD, so detection stays confident and
# only N50 changes. Base N50 15,463: x0.30 -> ~4,600 (warn), x0.045 -> ~700 (fail).
s30 <- shorten(ont, 0.30); s045 <- shorten(ont, 0.045)
add("ONT_N50_WARN", s30, platform = "ont", verdict = "WARN",
    note = sprintf("every read cut to 30%%: read N50 %d, under 5,000", n50(nchar(s30$s))))
add("ONT_N50_FAIL", s045, platform = "ont", verdict = "FAIL",
    note = sprintf("every read cut to 4.5%%: read N50 %d, under 1,000", n50(nchar(s045$s))))
add("ONT_N50_AMPLICON", s045, platform = "ont", sample_type = "amplicon", verdict = "PASS",
    note = "same reads as ONT_N50_FAIL, declared amplicon: N50 not judged")

# Quality ------------------------------------------------------------------------
# NovaSeq bins: F = Q37, ':' = Q25. Turning a fraction p of F into ':' lowers
# Q30 without making any base low enough (< Q15) for fastp to discard reads.
q30_base <- mean(unlist(lapply(c(ilm1$q, ilm2$q), function(q) phred(q) >= 30)))
for (target in c(0.72, 0.60)) {
  p  <- 1 - target / q30_base
  r1 <- ilm1; r2 <- ilm2
  r1$q <- swap_chars(r1$q, "F", ":", p); r2$q <- swap_chars(r2$q, "F", ":", p)
  add(sprintf("ILM_Q30_%d", round(100 * target)), r1, r2, platform = "illumina",
      verdict = if (target < 0.70) "FAIL" else "WARN",
      note = sprintf("Q30 lowered from %.1f%% to about %.0f%% (warn below 75%%, fail below 70%%)",
                     100 * q30_base, 100 * target))
}
# ONT: lower every quality by d, choosing d to land the NanoPlot-style mean
# read quality (base ~16.6) near the target.
for (target in c(8.5, 5.5)) {
  d  <- which.min(abs(vapply(0:40, function(dd) mean_q(lower_q(subset_fq(ont, 1:2000), dd)), 0) - target)) - 1L
  lq <- lower_q(ont, d)
  achieved <- mean_q(lq)
  # The case must land on the side of the limit its expected verdict says.
  if ((target < 7) != (achieved < 7) || achieved >= 10) stop(sprintf("ONT quality case missed its target: %.1f", achieved))
  add(sprintf("ONT_QUAL_%s", sub("\\.", "_", target)), lq, platform = "ont",
      verdict = if (target < 7) "FAIL" else "WARN",
      note = sprintf("every base quality lowered by %d: mean read quality %.1f (warn below 10, fail below 7)", d, mean_q(lq)))
}
ph <- ont; ph$q <- vapply(nchar(ont$q), function(n) strrep("5", n), "")
add("ONT_PLACEHOLDER_Q", ph, platform = "ont", verdict = "WARN", flags = "QUALITY_UNINFORMATIVE",
    note = paste("every base Q20: placeholder quality. Headers still confirm ONT (high confidence).",
                 "DESIRED answer WARN: quality cannot be assessed, so no PASS (D1). Current rules may give PASS",
                 "(quality metrics 'not judged', nothing else wrong); if so, that is a design finding for the user."))

# Duplication --------------------------------------------------------------------
half <- n_ilm %/% 2
dup <- function(fq) { j <- (half + 1):(2 * half); fq$s[j] <- fq$s[1:half]; fq$q[j] <- fq$q[1:half]; fq }
d1 <- dup(ilm1); d2 <- dup(ilm2)
add("ILM_DUP_DNA", d1, d2, platform = "illumina", verdict = "WARN",
    note = "second half of the pairs copies the first half: duplication about 50%, over the 30% warn limit")
add("ILM_DUP_RNA", d1, d2, platform = "illumina", molecule = "rna", verdict = "PASS",
    note = "same reads, declared rna: duplication not judged")

# Adapters -----------------------------------------------------------------------
# Real adapter read-through: for 40% of pairs, make the DNA fragment 80 bp, so
# each read runs off the fragment into its adapter (then NovaSeq's poly-G).
A1 <- "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG"
A2 <- "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT"
a1 <- ilm1; a2 <- ilm2
# Every read keeps its own length (ENA serves some already shortened), so
# only reads longer than the fragment can run into adapter.
hit <- which(runif(n_ilm) < 0.4 & nchar(ilm1$s) > 90 & nchar(ilm2$s) > 90)
for (k in hit) {
  frag <- substr(ilm1$s[k], 1, 80)
  a1$s[k] <- substr(paste0(frag, A1, strrep("G", 151)), 1, nchar(ilm1$s[k]))
  a2$s[k] <- substr(paste0(revcomp(frag), A2, strrep("G", 151)), 1, nchar(ilm2$s[k]))
}
if (length(hit) / n_ilm < 0.35) stop(sprintf("Adapter case missed its target: %.0f%% of pairs", 100 * length(hit) / n_ilm))
add("ILM_ADAPTER", a1, a2, platform = "illumina", verdict = "WARN",
    note = sprintf("%d of %d pairs (%.0f%%) have an 80 bp fragment, so reads run into adapter: adapter-trimmed reads over the 30%% warn limit",
                   length(hit), n_ilm, 100 * length(hit) / n_ilm))

# N bases ------------------------------------------------------------------------
add("ILM_N_WARN", add_ns(ilm1, 0.015), add_ns(ilm2, 0.015), platform = "illumina", verdict = "WARN",
    note = "1.5% of bases set to N: N rate over the 1% warn limit; fastp discards only the few reads with > 5 N, so passed-filter stays above 90%")
add("ILM_N_FAIL", add_ns(ilm1, 0.06), add_ns(ilm2, 0.06), platform = "illumina", verdict = "FAIL",
    note = "6% of bases set to N: N rate over the 5% fail limit (fastp also discards most reads)")

# Files and metadata -------------------------------------------------------------
add("ILM_EMPTY", list(h = character(0), s = character(0), q = character(0)), platform = "illumina",
    verdict = "FAIL", note = "empty file: no reads")
ctrl_ilm <- c(path("CTRL_ILM", 1), path("CTRL_ILM", 2))
add("ILM_DECLARED_ONT", files = ctrl_ilm, platform = "ont", verdict = "WARN",
    flags = "LAYOUT_CONFLICT,FASTQ_SHARED",
    note = "control Illumina reads declared ont: metadata contradicts reads (screen also flags ont + 2 files and the shared file)")
add("ILM_UNDECLARED", files = ctrl_ilm, platform = "", verdict = "WARN", flags = "FASTQ_SHARED",
    note = "control Illumina reads, no platform declared: platform inferred, so at most WARN")
o2 <- read_fq(other, n_ilm)
add("ILM_MIXED_MATES", ilm1, o2, platform = "illumina", verdict = "WARN", flags = "MATE_NAMES_DIFFER",
    note = "R1 from ECE107, R2 from a different isolate (bacass A1403KPN): read names do not match")
add("ILM_SWAPPED_FILES", files = c(path("ILM_SWAP", 2), path("ILM_SWAP", 1)), platform = "illumina",
    verdict = "PASS", flags = "MATES_POSSIBLY_SWAPPED",
    note = "R2 file given as fastq_1: flagged by the screen; metadata flags do not lower the verdict today (open question)")
write_fq(ilm1, path("ILM_SWAP", 1)); write_fq(ilm2, path("ILM_SWAP", 2))
add("ONT_TWO_FILES", subset_fq(ont, 1:(n_ont %/% 2)), subset_fq(ont, (n_ont %/% 2 + 1):n_ont), platform = "ont",
    verdict = "UNCERTAIN", flags = "TWO_FILES_LONG_READS,LAYOUT_CONFLICT",
    note = "nanopore reads split into two files: kept out of NanoPlot, no QC")
u <- ont; u$s <- chartr("T", "U", u$s)
add("ONT_U_DNA", u, platform = "ont", molecule = "dna", verdict = "WARN", flags = "U_BASES",
    note = "T written as U, as direct RNA is, but molecule declared dna: contradiction")
add("ONT_U_RNA", files = path("ONT_U_DNA"), platform = "ont", molecule = "rna", verdict = "PASS",
    flags = "U_BASES,FASTQ_SHARED", note = "same reads, declared rna: consistent")

# Isolated cases: each may stop a run, so each runs alone ------------------------
r2s <- subset_fq(ilm2, 1:100000)
add("ILM_R2_SHORT", ilm1, r2s, platform = "illumina", verdict = "WARN", flags = "MATE_COUNT_MISMATCH",
    batch = "isolated_r2_short",
    note = "R2 has 100,000 reads, R1 120,000. Desired: WARN with MATE_COUNT_MISMATCH. If fastp stops the run instead, that is a finding")
trunc_src <- path("CTRL_ILM", 1)
trunc_dst <- path("ILM_TRUNCATED_GZ")
bytes <- readBin(trunc_src, "raw", file.size(trunc_src))
writeBin(bytes[seq_len(length(bytes) %/% 2)], trunc_dst)
add("ILM_TRUNCATED_GZ", files = trunc_dst, platform = "illumina", verdict = "FAIL", flags = "FILE_UNREADABLE",
    note = paste("gzip cut in half. Since D45 the integrity check catches it: a FAIL row, and the rest of the batch",
                 "carries on (rounds 1 and 1-re-run ran it alone, when it still stopped the run: 'RUN_STOPS')."))

# Write --------------------------------------------------------------------------
exp <- do.call(rbind, rows)
write.table(exp, file.path(out, "expected.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
for (b in unique(exp$batch)) {
  e <- exp[exp$batch == b, ]
  write.csv(data.frame(sample = e$sample, fastq_1 = e$fastq_1, fastq_2 = e$fastq_2, platform = e$platform,
                       sample_type = e$sample_type, molecule = e$molecule),
            file.path(out, sprintf("samplesheet_%s.csv", b)), row.names = FALSE, quote = FALSE, na = "")
}
message(sprintf("\n%d cases written to %s (batches: %s)", nrow(exp), out, paste(unique(exp$batch), collapse = ", ")))
