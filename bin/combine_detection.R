#!/usr/bin/env Rscript
#*******************************************************************************
# combine_detection.R: one platform call for a sample in several files (D75)
#
# WHAT IT DOES  Joins the per-lane detection rows of one sample (a lane is a
#               mate pair, or one file of single-end reads), given in lane
#               order, into one row in detect_platform.R's format.
# RULES         Every lane must name the same platform; the sample's
#               confidence is the LOWEST lane's (one weak lane weakens the
#               call). Lanes that disagree, including one lane that cannot be
#               identified, make the sample uncertain with FILES_DISAGREE:
#               a mixed-up file is more likely than a mixed platform.
#               Layout must agree too. Read and file counts add up; the
#               read-start share (D65) is averaged, weighted by probed reads.
#
# Usage:  combine_detection.R <sample> <out.tsv> <lane1.tsv> <lane2.tsv> ...
#*******************************************************************************
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) stop("Usage: combine_detection.R <sample> <out.tsv> <lane1.tsv> <lane2.tsv> ...")
sample_id <- args[1]; out_tsv <- args[2]
lanes <- do.call(rbind, lapply(args[-(1:2)], function(f)
  read.delim(f, colClasses = "character", na.strings = "", quote = "", check.names = FALSE)))
k <- nrow(lanes)

CONF_ORDER <- c("high", "medium", "low", "uncertain")
# "lane" for mate pairs or mixed parts; "file" when every part is one file
# (a nanopore barcode folder, single-end lanes), so the reasons read right.
unit <- if (all(lanes$n_files == "1")) "file" else "lane"
Unit <- paste0(toupper(substr(unit, 1, 1)), substring(unit, 2))
num <- function(x) suppressWarnings(as.numeric(x))
split_flags <- function(x) if (is.na(x) || !nzchar(x)) character(0) else strsplit(x, ",")[[1]]

platforms <- unique(lanes$detected_platform)
layouts   <- unique(lanes$detected_layout)
flags     <- unique(unlist(lapply(lanes$detection_flags, split_flags)))
calls     <- sprintf("%s %d: %s (%s)", unit, seq_len(k), lanes$detected_platform, lanes$confidence)

if (length(platforms) == 1) {
  platform   <- platforms
  confidence <- CONF_ORDER[max(match(lanes$confidence, CONF_ORDER))]
  same_why   <- length(unique(lanes$reasons)) == 1
  reasons    <- paste0(if (unit == "lane") sprintf("Detected per lane (%d lanes, %d files); every lane says %s. ", k, sum(num(lanes$n_files)), platform)
                       else sprintf("Detected per file (%d files); every file says %s. ", k, platform),
                       if (same_why) lanes$reasons[1]
                       else paste(sprintf("%s %d: %s", Unit, seq_len(k), lanes$reasons), collapse = " "))
  if (confidence != CONF_ORDER[min(match(lanes$confidence, CONF_ORDER))])
    reasons <- paste(reasons, sprintf("Confidence is the weakest %s's: %s.", unit, paste(calls, collapse = "; ")))
} else {
  platform   <- "uncertain"
  confidence <- "uncertain"
  flags      <- unique(c(flags, "FILES_DISAGREE"))
  reasons    <- sprintf("The %ss of this sample look like different data (%s), so the platform cannot be trusted. Check that every file of the sample belongs to it.", unit,
                        paste(calls, collapse = "; "))
}
layout <- if (length(layouts) == 1) layouts else "uncertain"
if (length(layouts) > 1)
  reasons <- paste(reasons, sprintf("The %ss do not agree on the layout (%s).", unit,
                                    paste(sprintf("%s %d: %s", unit, seq_len(k), lanes$detected_layout), collapse = "; ")))

probe <- num(lanes$probe_reads); start <- num(lanes$start_share)
ok    <- !is.na(probe) & !is.na(start) & probe > 0
out <- data.frame(
  sample            = sample_id,
  n_files           = sum(num(lanes$n_files)),
  n_reads           = format(sum(num(lanes$n_reads)), scientific = FALSE),
  detected_platform = platform,
  detected_layout   = layout,
  confidence        = confidence,
  detection_flags   = paste(flags, collapse = ","),
  reasons           = reasons,
  stringsAsFactors  = FALSE
)
# Per-file evidence columns: the lanes' values in lane order, " / " between lanes.
for (col in c("max_len", "length_spread", "header_evidence", "quality_values", "u_bases", "adapters"))
  out[[col]] <- paste(ifelse(is.na(lanes[[col]]), "", lanes[[col]]), collapse = " / ")
out$start_share <- if (any(ok)) format(sum(start[ok] * probe[ok]) / sum(probe[ok]), digits = 4) else ""
out$probe_reads <- if (any(ok)) format(sum(probe[ok]), scientific = FALSE) else ""
write.table(out, out_tsv, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
