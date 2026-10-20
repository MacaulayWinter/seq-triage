#!/usr/bin/env Rscript
#*******************************************************************************
# score_layer_D.R: score the contamination screen against known answers
#                  (CLAUDE.md D43). Pass marks were fixed before any run:
#
# ESTIMATE  The judged metric (host_fraction or off_target_fraction) must be
#           within max(20% of the true value, 0.2 percentage points) of the
#           true fraction, and any true contamination must be DETECTED
#           (estimate > 0). Controls (true 0%): off-target below 5%.
#           Corrected before any real run: the first version used a 2 pp
#           floor, and the local rehearsal showed it let a complete miss at
#           0.5% count as a success (CLAUDE.md D43).
# VERDICT   Must equal expected_verdict where score_verdict = "yes".
#           "no": estimate scored, verdict not (sits on a limit).
#           "limit": the method is EXPECTED to miss it (another strain of
#           the same species); a PASS documents the limit.
#
# Usage:  Rscript validation/layer_D/score_layer_D.R <unpacked_results_dir> <out.tsv>
#         (the folder layer_D_results.tgz unpacks to: results/, mixtures/)
#*******************************************************************************

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("Usage: score_layer_D.R <unpacked_results_dir> <out.tsv>")
d  <- args[1]
rt <- function(f) read.delim(file.path(d, f), colClasses = "character", na.strings = "", quote = "")
truth <- rt("mixtures/truth.tsv")
summ  <- rt("results/triage_summary.tsv")
mets  <- rt("results/triage_metrics.tsv")
pct   <- function(x) as.numeric(sub("%", "", x)) / 100

rows <- lapply(seq_len(nrow(truth)), function(i) {
  t <- truth[i, ]
  m <- mets[mets$sample == t$sample, ]
  est_raw <- m$value[m$metric == t$metric]
  est <- if (length(est_raw) && !is.na(est_raw) && nzchar(est_raw)) pct(est_raw) else NA_real_
  tf  <- as.numeric(t$true_fraction)
  tol <- max(0.002, 0.2 * tf)
  est_ok <- if (is.na(est)) FALSE else if (tf == 0) est < 0.05 else est > 0 && abs(est - tf) <= tol
  got <- summ$verdict[summ$sample == t$sample]
  got <- if (length(got)) got else "MISSING"
  verdict_ok <- switch(t$score_verdict, yes = got == t$expected_verdict, no = NA, limit = got == "PASS")
  result <- if (t$score_verdict == "limit") {
    if (isTRUE(verdict_ok)) "LIMIT_CONFIRMED" else "LIMIT_NOT_SEEN"
  } else if (est_ok && !isFALSE(verdict_ok)) "PASS" else "FAIL"
  if (t$score_verdict == "limit") est_ok <- NA   # the estimate is expected to be wrong here
  data.frame(sample = t$sample, mixture = t$mixture, platform = t$platform, metric = t$metric,
             true = sprintf("%.1f%%", 100 * tf), estimate = if (is.na(est)) "n/a" else sprintf("%.1f%%", 100 * est),
             tolerance = if (tf == 0) "< 5%" else sprintf("+/- %.2f pp", 100 * tol),
             estimate_ok = est_ok, expected = t$expected_verdict, got = got, verdict_ok = verdict_ok,
             dominant = { x <- m$value[m$metric == "dominant_species"]; if (length(x)) x else "" },
             result = result)
})
s <- do.call(rbind, rows)
write.table(s, args[2], sep = "\t", quote = FALSE, row.names = FALSE, na = "")
print(s[, c("sample", "true", "estimate", "tolerance", "expected", "got", "result")], row.names = FALSE, right = FALSE)
scored <- s$result %in% c("PASS", "FAIL")
cat(sprintf("\nLayer D: %d / %d scored samples PASS; limit cases: %s\n", sum(s$result == "PASS"), sum(scored),
            paste(sprintf("%s=%s", s$sample[!scored], s$result[!scored]), collapse = ", ")))
if (any(s$result == "FAIL")) quit(status = 1)
