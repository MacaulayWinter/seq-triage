#*******************************************************************************
# lib_fastq.R: small FASTQ helpers shared by the validation scripts
#
# Base R only, and deliberately independent of every tool the pipeline uses
# (seqkit, fastp, NanoPlot), so the validation's own numbers are a genuine
# second opinion (CLAUDE.md D37, layer B).
#*******************************************************************************

# WHAT IT DOES  Reads up to n records of a FASTQ (gzipped or not).
# RETURNS       list(h = headers, s = sequences, q = quality strings)
read_fq <- function(path, n = -1L) {
  con <- gzfile(path, "r")
  on.exit(close(con))
  x <- readLines(con, n = if (n < 0) -1L else 4L * n, warn = FALSE)
  if (length(x) %% 4 != 0) stop("Truncated FASTQ record in ", path)
  i <- seq_len(length(x) %/% 4)
  list(h = x[4 * i - 3], s = x[4 * i - 2], q = x[4 * i])
}

# WHAT IT DOES  Writes records to a gzipped FASTQ.
write_fq <- function(fq, path) {
  con <- gzfile(path, "w")
  on.exit(close(con))
  if (length(fq$h)) writeLines(c(rbind(fq$h, fq$s, "+", fq$q)), con)
  invisible(path)
}

subset_fq <- function(fq, i) list(h = fq$h[i], s = fq$s[i], q = fq$q[i])

# Phred scores of one quality string (Sanger / Illumina 1.8+ encoding).
phred <- function(q) utf8ToInt(q) - 33L

# WHAT IT DOES  One read's average quality, the way NanoPlot defines it:
#               average the per-base ERROR probabilities, then convert that
#               mean back to a Phred score. Always lower than the plain mean
#               of the Q values.
read_mean_q <- function(q) -10 * log10(mean(10^(-phred(q) / 10)))

# WHAT IT DOES  A file's mean read quality as NanoPlot reports it (mean_qual):
#               each read's quality (read_mean_q), then THOSE are averaged as
#               error probabilities and converted back. Found in validation
#               round 1 (2026-10-20): it matches NanoPlot to the decimal it
#               prints, while the plain mean of per-read Q values does not
#               (up to 1.9 Q higher on R10 data, where read quality varies a lot).
nanoplot_mean_q <- function(quals) {
  rq <- vapply(quals, read_mean_q, 0, USE.NAMES = FALSE)
  -10 * log10(mean(10^(-rq / 10)))
}

# WHAT IT DOES  Length N50: the length L such that reads of length >= L hold
#               at least half of all bases.
n50 <- function(len) {
  len <- sort(len, decreasing = TRUE)
  len[which(cumsum(as.numeric(len)) >= sum(as.numeric(len)) / 2)[1]]
}
