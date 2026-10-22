#!/usr/bin/env Rscript
#*******************************************************************************
# score_layer_E.R: sample-type calls against ENA's library records (D65)
#
# PASS MARK  (fixed before any run) ZERO wrong sample-type calls at high
#            confidence. Truth: ENA library_strategy/source of the stratum
#            (WGS GENOMIC = isolate, WGS METAGENOMIC = metagenome,
#            AMPLICON = amplicon), corrected only by a documented hand
#            review (validation/adjudications_type.tsv, if present).
#            "uncertain" is never wrong.
# REPORTED   Correct high / medium calls per stratum; medium calls that are
#            wrong (allowed, but counted).
#
# Usage:  Rscript validation/layer_E/score_layer_E.R <panel.tsv> <unpacked_results_dir> <out.tsv>
#*******************************************************************************
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) stop("Usage: score_layer_E.R <panel.tsv> <unpacked_results_dir> <out.tsv>")
rt <- function(f) read.delim(f, colClasses = "character", na.strings = "", quote = "")
panel <- rt(args[1]); summ <- rt(file.path(args[2], "results", "triage_summary.tsv"))
summ$run <- sub("^.*_", "", summ$sample)
adj_file <- file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "..", "adjudications_type.tsv")
adj <- if (file.exists(adj_file)) rt(adj_file) else data.frame(run = character(0), adjudicated_truth = character(0))
m <- merge(summ[, c("run", "sample", "route", "detected_type", "type_confidence", "type_reason")],
           panel[, c("run_accession", "stratum", "type_truth", "instrument_model", "scientific_name")],
           by.x = "run", by.y = "run_accession")
m$truth <- m$type_truth; h <- match(m$run, adj$run); m$truth[!is.na(h)] <- adj$adjudicated_truth[h[!is.na(h)]]
m$outcome <- ifelse(m$detected_type == "uncertain", "uncertain",
             ifelse(m$detected_type == m$truth, paste("correct,", m$type_confidence), paste("WRONG,", m$type_confidence)))
m$wrong_high <- m$outcome == "WRONG, high"
m <- m[order(m$stratum, m$run), c("stratum", "run", "instrument_model", "scientific_name", "type_truth", "truth",
                                  "route", "detected_type", "type_confidence", "outcome", "wrong_high", "type_reason")]
write.table(m, args[3], sep = "\t", quote = FALSE, row.names = FALSE, na = "")
print(table(m$stratum, m$outcome))
cat(sprintf("\nLayer E: %d runs scored. Wrong at high confidence: %d (pass mark 0) -> %s\n", nrow(m), sum(m$wrong_high),
            if (sum(m$wrong_high) == 0) "PASS" else "FAIL"))
if (any(m$wrong_high)) print(m[m$wrong_high, c("run", "truth", "detected_type", "type_reason")], row.names = FALSE)
if (any(m$wrong_high)) quit(status = 1)
