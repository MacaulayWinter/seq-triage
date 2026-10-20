#!/usr/bin/env Rscript
# End-to-end check for `-profile test`
#
# Run the pipeline first, then point this at its results folder:
#   nextflow run main.nf -profile test,docker --outdir results
#   Rscript tests/e2e/check_test_profile.R results
#
# Checks, each printed as PASS/FAIL; exits 1 if any fails, so CI can gate on it:
#   1. every samplesheet row has a verdict (no sample silently lost)
#   2. each verdict matches tests/e2e/expected_verdicts.tsv
#   3. no FASTQ was published to results/ (reads never belong there)
#   4. every tool reported its version (reproducibility layer 4)
# Run from the repository root.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: check_test_profile.R <results_dir>")
res <- args[1]
read_tsv <- function(f) read.delim(f, colClasses = "character", na.strings = "", check.names = FALSE)

results <- character(0)
check <- function(name, ok, detail = "") {
  results[name] <<- if (ok) "PASS" else "FAIL"
  cat(sprintf("[%s] %s%s\n", results[name], name, if (!ok && nzchar(detail)) paste0("\n       ", detail) else ""))
}

sheet    <- read.csv("assets/samplesheet_test.csv", colClasses = "character")
expected <- read_tsv("tests/e2e/expected_verdicts.tsv")
summary  <- read_tsv(file.path(res, "triage_summary.tsv"))

# 1. Completeness
lost <- setdiff(sheet$sample, summary$sample)
check("every samplesheet row has a verdict", length(lost) == 0 && nrow(summary) == nrow(sheet),
      paste("missing:", paste(lost, collapse = ", ")))

# 2. Verdicts
stopifnot(setequal(expected$sample, sheet$sample))
got  <- summary$verdict[match(expected$sample, summary$sample)]
diff <- which(is.na(got) | got != expected$expected_verdict)
check("verdicts match expected_verdicts.tsv", length(diff) == 0,
      paste(sprintf("%s: got %s, expected %s", expected$sample[diff], got[diff],
                    expected$expected_verdict[diff]), collapse = "; "))

# 3. No reads published
fq <- list.files(res, pattern = "\\.f(ast)?q(\\.gz)?$", recursive = TRUE)
check("no FASTQ files published", length(fq) == 0, paste(fq, collapse = ", "))

# 4. Versions recorded for every tool
v <- read_tsv(file.path(res, "pipeline_info", "software_versions.tsv"))
need <- c("seqkit", "fastp", "fastqc", "NanoPlot", "multiqc", "R")
absent <- setdiff(need, v$tool[!is.na(v$version) & nzchar(v$version)])
check("every tool reported a version", length(absent) == 0, paste("missing:", paste(absent, collapse = ", ")))

n_fail <- sum(results == "FAIL")
cat(sprintf("\n%d / %d checks passed\n", length(results) - n_fail, length(results)))
if (n_fail > 0) quit(status = 1)
