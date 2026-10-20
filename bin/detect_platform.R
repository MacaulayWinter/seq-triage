#!/usr/bin/env Rscript
#*******************************************************************************
# detect_platform.R: infer sequencing platform and layout from seqkit stats
#
# STATUS        WALKING-SKELETON STUB. Every sample is reported as
#               "uncertain", the conservative default (CLAUDE.md D1), so the
#               pipeline's routing and reporting can be tested end to end
#               before the detection rules exist. The output columns are
#               final: the real rules will fill the same columns.
#               Only n_reads is computed for real, because the "empty"
#               branch depends on it.
#
# Usage:  detect_platform.R <sample_id> <seqkit_stats.tsv> <output.tsv>
#*******************************************************************************

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: detect_platform.R <sample_id> <seqkit_stats.tsv> <output.tsv>")
}
sample_id <- args[1]
stats_tsv <- args[2]
out_tsv   <- args[3]

# seqkit stats --tabular --all: one row per FASTQ file. check.names = FALSE
# keeps column names such as "Q30(%)" exactly as seqkit writes them.
stats <- read.delim(stats_tsv, check.names = FALSE, stringsAsFactors = FALSE)
if (nrow(stats) < 1) stop("No rows in ", stats_tsv)

out <- data.frame(
  sample            = sample_id,
  n_files           = nrow(stats),
  n_reads           = stats$num_seqs[1],     # reads in fastq_1 (pairs, if paired)
  detected_platform = "uncertain",
  detected_layout   = "uncertain",
  confidence        = "uncertain",
  detection_flags   = "",
  reasons           = "Detection rules not implemented yet (walking-skeleton stub); reported as uncertain, the conservative default.",
  stringsAsFactors  = FALSE
)
write.table(out, out_tsv, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
