#!/usr/bin/env Rscript
#*******************************************************************************
# triage_report.R: per-sample verdicts and a plain-language summary
#
# WHAT IT DOES  Joins the metadata screen (every samplesheet row) with the
#               per-sample detection and routing state, applies the verdict
#               ceiling (CLAUDE.md D2), and writes a TSV and a Markdown report.
# GATE          Stops if any samplesheet row is missing from the per-sample
#               results, rather than writing a shorter table. Same idea as
#               the author's 116_check_rerun_complete.sh: a sample must never
#               vanish silently between stages.
# STATUS        WALKING SKELETON. The verdict ceiling is real. QC metrics
#               from fastp/NanoPlot and the metric guide are not applied yet;
#               the report says so.
#
# Usage:  triage_report.R <metadata_screen.tsv> <samples.tsv> <metric_guide.tsv> <out_prefix>
#*******************************************************************************

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: triage_report.R <metadata_screen.tsv> <samples.tsv> <metric_guide.tsv> <out_prefix>")
}
read_tsv <- function(f) read.delim(f, colClasses = "character", na.strings = "", check.names = FALSE)
screen  <- read_tsv(args[1])
samples <- read_tsv(args[2])
guide   <- read_tsv(args[3])
prefix  <- args[4]


# Completeness gate ------------------------------------------------------------
missing <- setdiff(screen$sample, samples$sample)
extra   <- setdiff(samples$sample, screen$sample)
if (length(missing) || length(extra)) {
  stop("Sample sets do not match between the samplesheet and the per-sample results.",
       if (length(missing)) paste0("\n  Missing results for: ", paste(missing, collapse = ", ")),
       if (length(extra))   paste0("\n  Results for samples not in the sheet: ", paste(extra, collapse = ", ")))
}
d <- merge(screen, samples, by = "sample", all = TRUE, sort = FALSE,
           suffixes = c("", ".pipeline"))
d <- d[match(screen$sample, d$sample), ]   # keep samplesheet order


# Verdict ceiling (D2) -------------------------------------------------------------
# WHAT IT DOES  The best verdict a sample can reach, given how well its reads
#               confirm its metadata. QC metrics may only lower it later.
# RULE          No PASS unless the reads independently confirm the metadata.
ceiling_for <- function(route, routed_by, confidence, declared, detected) {
  if (route == "empty")      return(c("FAIL",      "No reads in the file."))
  if (routed_by == "none")   return(c("UNCERTAIN", "Platform could not be determined from the reads, and none was declared. No QC was run."))
  if (routed_by == "declared") {
    if (confidence == "uncertain")
      return(c("WARN", sprintf("Reads were inconclusive; QC used the declared platform (%s), which is unconfirmed.", declared)))
    return(c("WARN", sprintf("Declared platform (%s) is not confirmed by the reads (detection confidence: %s).", declared, confidence)))
  }
  # routed_by == "detected": only high-confidence detection routes this way.
  if (is.na(declared))     return(c("WARN", sprintf("Platform inferred from the reads (%s); none was declared.", detected)))
  if (declared == detected) return(c("PASS", sprintf("Reads confirm the declared platform (%s).", declared)))
  c("WARN", sprintf("Metadata says %s but the reads indicate %s. Review this sample.", declared, detected))
}

v <- t(mapply(ceiling_for, d$route, d$routed_by, d$confidence, d$declared_platform, d$detected_platform))
d$verdict        <- v[, 1]
d$verdict_reason <- v[, 2]


# Write ------------------------------------------------------------------------
cols <- c("sample", "verdict", "verdict_reason", "declared_platform", "declared_platform_raw",
          "detected_platform", "detected_layout", "confidence", "route", "routed_by",
          "n_files", "n_reads", "declared_type", "n_flags", "flags", "detection_flags", "reasons")
write.table(d[, cols], paste0(prefix, ".tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")

show <- function(x) ifelse(is.na(x) | x == "", "-", x)
md <- c(
  "# seq-triage summary", "",
  sprintf("%d sample(s): %s.", nrow(d),
          paste(sprintf("%d %s", table(factor(d$verdict, c("PASS", "WARN", "UNCERTAIN", "FAIL"))),
                        c("PASS", "WARN", "UNCERTAIN", "FAIL")), collapse = ", ")), "",
  "> **Walking-skeleton build.** Platform detection is a stub that reports every sample as",
  "> *uncertain*, and QC metrics are not yet applied to the verdicts. Verdicts below reflect",
  "> metadata and routing only.", "",
  "| Sample | Verdict | Declared | Detected | Confidence | QC route | Metadata flags |",
  "|---|---|---|---|---|---|---|",
  sprintf("| %s | **%s** | %s | %s | %s | %s | %s |", d$sample, d$verdict, show(d$declared_platform),
          show(d$detected_platform), show(d$confidence), show(d$route), show(d$n_flags)), "",
  "## Per-sample detail", "",
  unlist(lapply(seq_len(nrow(d)), function(i) c(
    sprintf("### %s: %s", d$sample[i], d$verdict[i]), "",
    paste0("- ", d$verdict_reason[i]),
    if (!is.na(d$flags[i])) paste0("- Metadata: ", strsplit(d$flags[i], " | ", fixed = TRUE)[[1]]),
    "")))
)
writeLines(md, paste0(prefix, ".md"))
message(sprintf("Report written for %d sample(s).", nrow(d)))
