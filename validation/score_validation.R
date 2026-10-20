#!/usr/bin/env Rscript
#*******************************************************************************
# score_validation.R: score layers A, B and C against their pass marks
#                     (CLAUDE.md D37). Pass marks were fixed before any run.
#
# A  Platform calls vs ENA's instrument platform. Pass mark: ZERO wrong calls
#    at high confidence. "uncertain" is never wrong; for platforms seq-triage
#    does not handle (PacBio, DNBSEQ, Element, Ion Torrent) it is the right answer.
# B  Independent recalculation (independent_metrics.R output). Pass mark:
#    every check PASS.
# C  Damaged samples vs expected.tsv. Pass mark: every expected verdict met and
#    every expected flag present (extra flags are allowed and listed).
#
# Usage:  Rscript validation/score_validation.R <panel.tsv> <panel_results_dir>
#           <faults_dir> <faults_results_dir> <layer_B.tsv> <out_dir>
#*******************************************************************************

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) stop("Usage: score_validation.R <panel.tsv> <panel_results> <faults_dir> <faults_results> <layer_B.tsv> <out_dir>")
rt <- function(f) read.delim(f, colClasses = "character", na.strings = "", quote = "")
panel  <- rt(args[1])
pres   <- rt(file.path(args[2], "triage_summary.tsv"))
exp    <- rt(file.path(args[3], "expected.tsv"))
fres   <- rt(file.path(args[4], "triage_summary.tsv"))
layerB <- rt(args[5])
out    <- args[6]
dir.create(out, showWarnings = FALSE, recursive = TRUE)

# Layer A ----------------------------------------------------------------------
pres$run <- sub("^.*_", "", pres$sample)
pres$stratum <- sub("_[^_]+$", "", pres$sample)
p <- panel[panel$role %in% c("panel", "backup"), ]
a <- merge(pres[, c("sample", "stratum", "run", "detected_platform", "confidence", "detected_layout", "detection_flags")],
           p[, c("run_accession", "truth", "library_layout", "instrument_model", "scientific_name")],
           by.x = "run", by.y = "run_accession")
handled <- c("illumina", "ont")
a$outcome <- ifelse(a$detected_platform == "uncertain",
                    ifelse(a$truth %in% handled, "uncertain (missed)", "uncertain (correct: not handled)"),
             ifelse(a$detected_platform == a$truth, paste("correct,", a$confidence),
                    paste("WRONG,", a$confidence)))
a$wrong_high <- grepl("^WRONG", a$outcome) & a$confidence == "high"
a$layout_ok  <- ifelse(a$detected_platform %in% handled & a$truth %in% handled,
                       (a$library_layout == "PAIRED") == (a$detected_layout == "paired"), NA)
a <- a[order(a$stratum, a$run), c("stratum", "run", "instrument_model", "scientific_name", "truth", "detected_platform",
                                  "confidence", "outcome", "wrong_high", "library_layout", "detected_layout",
                                  "layout_ok", "detection_flags")]
write.table(a, file.path(out, "layer_A.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")

# Layer C ----------------------------------------------------------------------
codes <- function(screen_flags, det_flags) {
  sc <- if (is.na(screen_flags)) character(0) else sub(":.*", "", strsplit(screen_flags, " | ", fixed = TRUE)[[1]])
  dt <- if (is.na(det_flags)) character(0) else strsplit(det_flags, ",")[[1]]
  unique(c(sc, dt))
}
c_rows <- lapply(seq_len(nrow(exp)), function(i) {
  e <- exp[i, ]
  if (e$batch != "main") return(NULL)   # isolated cases are scored from their own run logs
  r <- fres[fres$sample == e$sample, ]
  got_flags <- if (nrow(r)) codes(r$flags, r$detection_flags) else character(0)
  want <- if (is.na(e$expected_flags)) character(0) else strsplit(e$expected_flags, ",")[[1]]
  data.frame(sample = e$sample, expected_verdict = e$expected_verdict,
             got_verdict = if (nrow(r)) r$verdict else "MISSING",
             missing_flags = paste(setdiff(want, got_flags), collapse = ","),
             extra_flags = paste(setdiff(got_flags, want), collapse = ","),
             result = if (nrow(r) && r$verdict == e$expected_verdict && all(want %in% got_flags)) "PASS" else "FAIL",
             qc_findings = if (nrow(r)) r$qc_findings else "", note = e$note)
})
cc <- do.call(rbind, c_rows)
write.table(cc, file.path(out, "layer_C.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")

# Summary ----------------------------------------------------------------------
line <- function(...) cat(sprintf(...), "\n", sep = "")
line("LAYER A  %d runs. Wrong at high confidence: %d (pass mark 0) -> %s", nrow(a), sum(a$wrong_high),
     if (sum(a$wrong_high) == 0) "PASS" else "FAIL")
print(table(a$stratum, a$outcome))
line("         Layout correct for %d of %d called runs", sum(a$layout_ok, na.rm = TRUE), sum(!is.na(a$layout_ok)))
line("LAYER B  %d checks: %d PASS, %d FAIL, %d MISSING -> %s", nrow(layerB), sum(layerB$result == "PASS"),
     sum(layerB$result == "FAIL"), sum(layerB$result == "MISSING"),
     if (all(layerB$result == "PASS")) "PASS" else "FAIL")
line("LAYER C  %d main-batch cases: %d PASS, %d FAIL -> %s", nrow(cc), sum(cc$result == "PASS"),
     sum(cc$result == "FAIL"), if (all(cc$result == "PASS")) "PASS" else "FAIL")
if (any(cc$result == "FAIL")) print(cc[cc$result == "FAIL", c("sample", "expected_verdict", "got_verdict", "missing_flags")], row.names = FALSE)
