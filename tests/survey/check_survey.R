#!/usr/bin/env Rscript
# Check bin/detect_platform.R against the real test-data survey
#
# The survey is 12 public samples (Illumina from four instrument
# generations, old and new nanopore, PacBio HiFi) chosen to cover the cases
# detection must get right. The DATA stays outside the repository; only
# the expected answers (expected_calls.tsv) are kept here. See CLAUDE.md
# "Test data survey" for where each file came from.
#
# Not run in CI (the data is ~600 MB). Run locally from the repository root:
#   Rscript tests/survey/check_survey.R ~/seqtriage-scratch/survey
# The survey folder must hold fastq/*.fastq.gz and seqkit_stats.tsv
# (seqkit stats -a -T over all of them).
# Exits non-zero if any sample differs.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: check_survey.R <survey_dir>")
survey <- path.expand(args[1])
stats  <- file.path(survey, "seqkit_stats.tsv")
if (!file.exists(stats)) stop("No seqkit_stats.tsv in ", survey)

exp <- read.delim("tests/survey/expected_calls.tsv", colClasses = "character", na.strings = "")
tmp <- tempfile("survey_")
dir.create(tmp)

run_one <- function(i) {
  e     <- exp[i, ]
  files <- file.path(survey, "fastq", strsplit(e$files, ",")[[1]])
  out   <- file.path(tmp, paste0(e$sample, ".tsv"))
  status <- system2("Rscript", c("--vanilla", "bin/detect_platform.R", e$sample, stats, out, files))
  if (status != 0) return(c(got = sprintf("script exited %d", status), expected = "", reason = ""))
  d <- read.delim(out, colClasses = "character", na.strings = "")
  sorted <- function(f) if (is.na(f) || f == "") "" else paste(sort(strsplit(f, ",")[[1]]), collapse = ",")
  c(got      = paste(d$detected_platform, d$detected_layout, d$confidence, sorted(d$detection_flags)),
    expected = paste(e$platform, e$layout, e$confidence, sorted(e$flags)),
    reason   = d$reasons)
}

r <- t(vapply(seq_len(nrow(exp)), run_one, c(got = "", expected = "", reason = "")))
res <- data.frame(sample = exp$sample, got = r[, "got"], expected = r[, "expected"],
                  result = ifelse(r[, "got"] == r[, "expected"], "PASS", "FAIL"), stringsAsFactors = FALSE)
print(res, row.names = FALSE, right = FALSE)
cat("\nReasons given:\n")
cat(sprintf("  %-17s %s\n", exp$sample, r[, "reason"]), sep = "")

n_fail <- sum(res$result == "FAIL")
cat(sprintf("\n%d / %d survey samples as expected\n", nrow(res) - n_fail, nrow(res)))
unlink(tmp, recursive = TRUE)
if (n_fail > 0) quit(status = 1)
