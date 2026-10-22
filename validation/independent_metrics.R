#!/usr/bin/env Rscript
#*******************************************************************************
# independent_metrics.R: second opinion on every QC number (CLAUDE.md D37, B)
#
# WHAT IT DOES  For each sample, recomputes counts and statistics straight
#               from its FASTQ in base R (validation/lib_fastq.R: no seqkit,
#               fastp or NanoPlot), then compares them with
#                 (1) what each tool reported (results/qc_values.tsv), and
#                 (2) what the report shows (results/triage_metrics.tsv).
#               (1) checks the tools were read correctly (right field, right
#               units, both mates counted); (2) checks the report's own
#               arithmetic and formatting.
# LIMITS        fastp's AFTER-filtering numbers cannot be recomputed without
#               re-implementing fastp's filters; they are checked for
#               consistency only (after <= before; passed = after).
# PASS MARKS    Fixed before the run: counts exact; NanoPlot mean read quality
#               within 0.1; fastp Q30 within 0.0001 (fastp's own rounding);
#               report values within the precision they are displayed at.
#
# Usage:  Rscript validation/independent_metrics.R <samplesheet.csv> <results_dir> <out.tsv>
#*******************************************************************************

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) stop("Usage: independent_metrics.R <samplesheet.csv> <results_dir> <out.tsv>")
script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(script_dir, "lib_fastq.R"))
sheet <- read.csv(args[1], colClasses = "character", na.strings = "")
res   <- args[2]
rt    <- function(f) read.delim(file.path(res, f), colClasses = "character", na.strings = "", quote = "")
qc    <- rt("qc_values.tsv")
summ  <- rt("triage_summary.tsv")
mets  <- rt("triage_metrics.tsv")

num <- function(x) as.numeric(gsub("[,%Q]", "", x))
tool_value <- function(s, tool, field, sum_files = FALSE) {
  v <- num(qc$value[qc$sample == s & qc$tool == tool & qc$field == field])
  if (!length(v)) NA_real_ else if (sum_files) sum(v) else v[1]
}
report_value <- function(s, metric) {
  v <- mets$value[mets$sample == s & mets$metric == metric]
  if (!length(v) || is.na(v)) NA_real_ else num(v)
}

out <- list()
check <- function(s, quantity, independent, pipeline, source, tol) {
  diff <- abs(independent - pipeline)
  out[[length(out) + 1]] <<- data.frame(
    sample = s, quantity = quantity, independent = signif(independent, 8), pipeline = signif(pipeline, 8),
    source = source, tolerance = tol, difference = signif(diff, 4),
    result = if (is.na(independent) || is.na(pipeline)) "MISSING" else if (diff <= tol) "PASS" else "FAIL")
}

for (i in seq_len(nrow(sheet))) {
  s     <- sheet$sample[i]
  files <- na.omit(c(sheet$fastq_1[i], sheet$fastq_2[i]))
  route <- summ$route[summ$sample == s]
  if (!length(route)) next
  # A file the integrity check rejected (D45) has no QC numbers to recalculate;
  # that it is rejected, and how, is layer C's test.
  if (route == "unreadable") next
  fqs   <- lapply(files, read_fq)
  len   <- unlist(lapply(fqs, function(f) nchar(f$s)))
  n_rd  <- length(len)
  n_bp  <- sum(as.numeric(len))
  n_N   <- sum(vapply(unlist(lapply(fqs, `[[`, "s")), function(x) lengths(regmatches(x, gregexpr("N", x, fixed = TRUE))), 0))
  message(sprintf("%-28s %s, %s reads", s, route, format(n_rd, big.mark = ",")))

  # seqkit runs on every sample
  check(s, "reads (all files)", n_rd, tool_value(s, "seqkit", "num_seqs", TRUE), "seqkit num_seqs", 0)
  check(s, "bases (all files)", n_bp, tool_value(s, "seqkit", "sum_len", TRUE), "seqkit sum_len",  0)
  check(s, "N bases",           n_N,  tool_value(s, "seqkit", "sum_n",   TRUE), "seqkit sum_n",    0)
  if (n_bp > 0 && route %in% c("illumina", "ont")) {
    check(s, "N rate (%)", 100 * n_N / n_bp, report_value(s, "n_base_rate"), "report n_base_rate",
          if (100 * n_N / n_bp < 1) 0.005 else 0.05)
  }

  if (route == "ont") {
    check(s, "reads",              n_rd,              tool_value(s, "NanoPlot", "number_of_reads"),    "NanoPlot number_of_reads", 0)
    check(s, "bases",              n_bp,              tool_value(s, "NanoPlot", "number_of_bases"),    "NanoPlot number_of_bases", 0)
    check(s, "read N50",           n50(len),          tool_value(s, "NanoPlot", "n50"),                "NanoPlot n50", 0)
    check(s, "median read length", median(len),       tool_value(s, "NanoPlot", "median_read_length"), "NanoPlot median_read_length", 0)
    # NanoPlot's definition (lib_fastq.R nanoplot_mean_q). Round 1 first used
    # the plain mean of per-read Q values, which disagreed by up to 1.9 Q on
    # R10 data; the corrected definition was confirmed against NanoPlot's
    # mean and median before being adopted (docs/validation.md, finding B1).
    check(s, "mean read quality",  nanoplot_mean_q(fqs[[1]]$q), tool_value(s, "NanoPlot", "mean_qual"), "NanoPlot mean_qual", 0.1)
    check(s, "reads (report)",     n_rd,              report_value(s, "number_of_reads"),              "report number_of_reads", 0)
    check(s, "read N50 (report)",  n50(len),          report_value(s, "read_n50"),                     "report read_n50", 0.5)
  }
  if (route == "illumina") {
    q30 <- sum(vapply(unlist(lapply(fqs, `[[`, "q")), function(q) sum(phred(q) >= 30), 0)) / n_bp
    check(s, "reads before filtering", n_rd, tool_value(s, "fastp", "summary.before_filtering.total_reads"), "fastp before_filtering.total_reads", 0)
    check(s, "bases before filtering", n_bp, tool_value(s, "fastp", "summary.before_filtering.total_bases"), "fastp before_filtering.total_bases", 0)
    check(s, "Q30 before filtering",   q30,  tool_value(s, "fastp", "summary.before_filtering.q30_rate"),    "fastp before_filtering.q30_rate", 1e-4)
    after  <- tool_value(s, "fastp", "summary.after_filtering.total_reads")
    passed <- tool_value(s, "fastp", "filtering_result.passed_filter_reads")
    check(s, "after <= before (0 = consistent)", 0, as.numeric(after > n_rd), "fastp after_filtering.total_reads", 0)
    check(s, "passed = after",                   after, passed, "fastp filtering_result.passed_filter_reads", 0)
    check(s, "reads after filtering (report)",   after, report_value(s, "total_reads"), "report total_reads", 0)
  }
}

o <- do.call(rbind, out)
write.table(o, args[3], sep = "\t", quote = FALSE, row.names = FALSE, na = "")
tab <- table(o$result)
message(sprintf("\nLayer B: %d checks: %s", nrow(o), paste(names(tab), tab, sep = " ", collapse = ", ")))
