#!/usr/bin/env Rscript
#*******************************************************************************
# score_layer_E_round7.R: score layer E round 7 against its pre-registered
#                         pass marks (CLAUDE.md D94)
#
# PASS MARKS    B. No library-structure flag (BISULFITE, BARCODE_READ,
#                  INDEX_READS, INTERLEAVED) on any run: every run is a
#                  G-dominated (G >= 50% in R1) genomic run, not bisulfite.
#               D. At least 8 such runs tested; fewer is INCONCLUSIVE.
# REPORTED      How many the rule before D90 would have flagged (computed from
#               the screened base shares), so the test is seen to have teeth.
# Usage:  Rscript score_layer_E_round7.R <panel.tsv> <unpacked_results_dir> <out_prefix>
#*******************************************************************************
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) stop("Usage: score_layer_E_round7.R <panel.tsv> <unpacked_results_dir> <out_prefix>")
rt <- function(f) read.delim(f, colClasses = "character", na.strings = "", quote = "")
panel <- rt(args[1]); summ <- rt(file.path(args[2], "results", "triage_summary.tsv"))
summ$run <- regmatches(summ$sample, regexpr("(SRR|ERR|DRR)[0-9]+", summ$sample))
m <- merge(summ, panel[, c("run_accession", "a_share", "c_share", "g_share", "t_share")], by.x = "run", by.y = "run_accession")
fl <- function(i) if (is.na(m$detection_flags[i])) character(0) else strsplit(m$detection_flags[i], ",")[[1]]
lib <- vapply(seq_len(nrow(m)), function(i) any(fl(i) %in% c("BISULFITE", "BARCODE_READ", "INDEX_READS", "INTERLEAVED")), TRUE)
cs <- as.numeric(m$c_share); gs <- as.numeric(m$g_share)
old_rule <- pmin(cs, gs) < 0.10 & pmin(cs, gs) < 0.5 * pmax(cs, gs)
write.table(m[, c("run", "g_share", "c_share", "detection_flags", "verdict")], paste0(args[3], "_scores.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "")
cat(sprintf("B. Library-structure flags on G-dominated genomic runs: %d (mark 0) -> %s\n", sum(lib), if (sum(lib) == 0) "PASS" else "FAIL"))
if (any(lib)) print(m[lib, c("run", "g_share", "detection_flags")], row.names = FALSE)
cat(sprintf("D. G-dominated runs tested: %d (mark 8) -> %s\n", nrow(m), if (nrow(m) >= 8) "PASS" else "INCONCLUSIVE"))
cat(sprintf("\nReported: the rule before D90 would have flagged %d of %d of these runs as bisulfite\n", sum(old_rule, na.rm = TRUE), nrow(m)))
if (any(lib)) quit(status = 1)
if (nrow(m) < 8) quit(status = 2)
