#!/usr/bin/env Rscript
# =============================================================================
# rerun_report.R  (D96 validation)
# Re-runs bin/triage_report.R on a saved validation run (its triage_summary.tsv
# and qc_values.tsv), with a pathogen list, so the flags can be scored on runs
# whose composition is known without re-running Kraken2. The summary carries
# every column the report reads from the metadata screen and from detection,
# so the two are split back out of it.
#
# Usage:  rerun_report.R <outcomes_dir> <pathogen_list.tsv> <out_dir>
# Writes: <out_dir>/triage_summary.tsv, _metrics.tsv, _pathogens.tsv, _summary.md
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) stop("Usage: rerun_report.R <outcomes_dir> <pathogen_list.tsv> <out_dir>")
run <- args[1]; plist <- args[2]; out <- args[3]
dir.create(out, showWarnings = FALSE, recursive = TRUE)
repo <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "..", ".."))

s <- read.delim(file.path(run, "triage_summary.tsv"), colClasses = "character", na.strings = "", check.names = FALSE, quote = "")
SCREEN_COLS <- c("sample", "declared_platform", "declared_platform_raw", "declared_type", "declared_molecule",
                 "declared_library", "declared_target", "n_flags", "flags")
DETECT_COLS <- c("sample", "detected_platform", "detected_layout", "confidence", "route", "routed_by",
                 "n_files", "n_reads", "detection_flags", "reasons")
w <- function(cols, f) write.table(s[, intersect(cols, names(s)), drop = FALSE], f, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
w(SCREEN_COLS, file.path(out, "screen.tsv"))
w(DETECT_COLS, file.path(out, "samples.tsv"))

status <- system2("Rscript", c("--vanilla", file.path(repo, "bin/triage_report.R"),
                               file.path(out, "screen.tsv"), file.path(out, "samples.tsv"),
                               file.path(repo, "assets/metric_guide.tsv"), file.path(run, "qc_values.tsv"),
                               file.path(repo, "assets/species_groups.tsv"), file.path(out, "triage"),
                               file.path(repo, "assets/species_genome_size.tsv"), plist,
                               shQuote("validation re-run")))
if (status != 0) stop("triage_report.R failed on ", run)
cat("re-ran the report on", run, "\n")
