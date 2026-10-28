#!/usr/bin/env Rscript
# Test for bin/screen_metadata.R
#
# Runs the screen on a synthetic samplesheet in which every row is built to
# trigger specific flags (or none), then compares the flag codes produced
# against expected_flags.tsv. The FASTQ paths point at no real files: the
# screen only reads the samplesheet, never the reads.
#
# Run from the repository root:
#   Rscript tests/screen_metadata/test_screen_metadata.R
# Exits non-zero if any sample differs, so CI can use it.

here     <- "tests/screen_metadata"
`%|%` <- function(a, b) if (length(a) && !is.na(a[1])) a[1] else b
out_tsv  <- tempfile(fileext = ".tsv")
status   <- system2("Rscript", c("--vanilla", "bin/screen_metadata.R",
                                 file.path(here, "samplesheet.csv"), out_tsv))
if (status != 0) stop("screen_metadata.R exited with status ", status)

got <- read.delim(out_tsv, colClasses = "character", na.strings = "")
exp <- read.delim(file.path(here, "expected_flags.tsv"), colClasses = "character", na.strings = "")

# Reduce "CODE: message | CODE: message" to sorted codes, e.g. "A,B".
codes <- function(f) {
  if (is.na(f)) return("")
  paste(sort(sub(":.*", "", strsplit(f, " | ", fixed = TRUE)[[1]])), collapse = ",")
}

stopifnot(identical(got$sample, exp$sample))
res <- data.frame(
  sample   = got$sample,
  got      = unname(vapply(got$flags, codes, "")),
  expected = ifelse(is.na(exp$expected), "", exp$expected),
  stringsAsFactors = FALSE
)
# D93: every name for a methylation library is read as "bisulfite".
for (k in c("LIB_BS", "LIB_WGBS"))
  res <- rbind(res, data.frame(sample = paste0(k, "_declared_library"), got = got$declared_library[got$sample == k] %|% "",
                               expected = "bisulfite", stringsAsFactors = FALSE))
res$result <- ifelse(res$got == res$expected, "PASS", "FAIL")
print(res, row.names = FALSE, right = FALSE)

n_fail <- sum(res$result == "FAIL")
cat(sprintf("\n%d / %d samples as expected\n", nrow(res) - n_fail, nrow(res)))
if (n_fail > 0) quit(status = 1)
