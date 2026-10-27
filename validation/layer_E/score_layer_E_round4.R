#!/usr/bin/env Rscript
#*******************************************************************************
# score_layer_E_round4.R: score layer E round 4 (D87) or round 5 (D90)
#                         against its pre-registered pass marks
#
# PASS MARKS    A. Sample type: zero wrong calls at high confidence (every
#                  stratum whose type_truth is not "none").
#               B. Library-structure specificity (D78): no BISULFITE or
#                  BARCODE_READ flag outside the bisulfite / single-cell
#                  strata, and no INDEX_READS or INTERLEAVED flag anywhere.
#               C. Bisulfite sensitivity (D78): BISULFITE on at least 7 of 8
#                  bisulfite runs (87.5% if fewer were fetched).
# REPORTED      BARCODE_READ on the single-cell runs (only 10x-style runs
#               carry a barcode read); the rRNA share of RNA strata (D81);
#               correct metagenome calls held at medium by D86 (its cost).
# Usage:  Rscript score_layer_E_round4.R <panel.tsv> <unpacked_results_dir> <out_prefix>
#*******************************************************************************
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) stop("Usage: score_layer_E_round4.R <panel.tsv> <unpacked_results_dir> <out_prefix>")
rt <- function(f) read.delim(f, colClasses = "character", na.strings = "", quote = "")
panel <- rt(args[1]); summ <- rt(file.path(args[2], "results", "triage_summary.tsv"))
mets <- rt(file.path(args[2], "results", "triage_metrics.tsv"))
summ$run <- regmatches(summ$sample, regexpr("(SRR|ERR|DRR)[0-9]+", summ$sample))
m <- merge(summ, panel[, c("run_accession", "stratum", "type_truth", "scientific_name")], by.x = "run", by.y = "run_accession")
flags <- function(i) if (is.na(m$detection_flags[i])) character(0) else strsplit(m$detection_flags[i], ",")[[1]]
has_flag <- function(f) vapply(seq_len(nrow(m)), function(i) f %in% flags(i), TRUE)

# A. sample type
typed <- m$type_truth != "none"
m$outcome <- ifelse(!typed, "not scored", ifelse(m$detected_type == "uncertain", "uncertain",
             ifelse(m$detected_type == m$type_truth, paste("correct,", m$type_confidence), paste("WRONG,", m$type_confidence))))
wrong_high <- typed & m$outcome == "WRONG, high"

# B. specificity of the library-structure flags
bs_stratum <- grepl("BISULFITE", m$stratum); sc_stratum <- grepl("SINGLE_CELL", m$stratum)
false_flag <- (has_flag("BISULFITE") & !bs_stratum) | (has_flag("BARCODE_READ") & !sc_stratum) |
              has_flag("INDEX_READS") | has_flag("INTERLEAVED")

# C. bisulfite sensitivity
n_bs <- sum(bs_stratum); hit_bs <- sum(bs_stratum & has_flag("BISULFITE"))
# 87.5% of the bisulfite runs fetched: 7 of 8 (round 4), 13 of 14 (round 5, D90).
need_bs <- ceiling(0.875 * n_bs)

# Reported
rrna <- mets[mets$metric == "rrna_fraction", c("sample", "value", "result")]
d86 <- grepl("k-mers recur|narrow GC range|too little sequence", m$type_reason)
m$d86_held_medium <- d86

write.table(m[order(m$stratum, m$run), c("stratum", "run", "scientific_name", "type_truth", "detected_type", "type_confidence", "outcome",
                                         "d86_held_medium", "detection_flags", "type_reason")],
            paste0(args[3], "_scores.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
write.table(rrna, paste0(args[3], "_rrna.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")

print(table(m$stratum[typed], m$outcome[typed]))
cat(sprintf("\nA. Wrong sample-type calls at high confidence: %d (mark 0) -> %s\n", sum(wrong_high), if (sum(wrong_high) == 0) "PASS" else "FAIL"))
if (any(wrong_high)) print(m[wrong_high, c("run", "stratum", "scientific_name", "detected_type")], row.names = FALSE)
cat(sprintf("B. Library-structure flags where they do not belong: %d (mark 0) -> %s\n", sum(false_flag), if (sum(false_flag) == 0) "PASS" else "FAIL"))
if (any(false_flag)) print(m[false_flag, c("run", "stratum", "detection_flags")], row.names = FALSE)
cat(sprintf("C. BISULFITE on bisulfite runs: %d of %d (mark %d) -> %s\n", hit_bs, n_bs, need_bs, if (n_bs > 0 && hit_bs >= need_bs) "PASS" else "FAIL"))
cat(sprintf("\nReported: BARCODE_READ on %d of %d single-cell runs (only 10x-style runs have a barcode read)\n",
            sum(sc_stratum & has_flag("BARCODE_READ")), sum(sc_stratum)))
cat(sprintf("Reported: D86 held %d correct metagenome call(s) at medium, and %d wrong one(s)\n",
            sum(d86 & typed & grepl("^correct", m$outcome)), sum(d86 & typed & grepl("^WRONG", m$outcome))))
cat(sprintf("Reported: rRNA share measured on %d RNA sample(s); see %s_rrna.tsv\n", sum(rrna$result != "not_available"), args[3]))
# D. Round 5 only (D90): the fix must meet the case it fixes. Runs with GC
# content of 75% or more stand in for no-signal (poly-G) reads; fewer than
# 2 of them among the tiny two-colour runs means the fix was not tested on
# fresh data, so the round is INCONCLUSIVE (not a pass) and more are drawn.
tiny <- grepl("TINY_2COLOUR", m$stratum)
if (any(tiny)) {
  gc <- suppressWarnings(as.numeric(sub("%", "", mets$value[match(paste(m$sample, "gc_content"), paste(mets$sample, mets$metric))])))
  n_polyg <- sum(tiny & !is.na(gc) & gc >= 75)
  cat(sprintf("D. G-dominated (GC >= 75%%) tiny two-colour runs: %d (mark 2) -> %s\n", n_polyg, if (n_polyg >= 2) "PASS" else "INCONCLUSIVE"))
  if (n_polyg < 2) quit(status = 2)
}
if (any(wrong_high) || any(false_flag) || n_bs == 0 || hit_bs < need_bs) quit(status = 1)
