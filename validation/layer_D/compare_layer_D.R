#!/usr/bin/env Rscript
#*******************************************************************************
# compare_layer_D.R: the same layer D samples run against two Kraken2
#                    databases, side by side (CLAUDE.md D49)
#
# WHAT IT DOES  Scores each run with score_layer_D.R (the same pre-registered
#               pass marks), then lines the two up sample by sample: estimate,
#               verdict, result, unclassified fraction and dominant species.
#               Answers: which results change when only the database changes?
# CHECKS        Both runs must hold the same samples with the same true
#               fractions (same mixtures, only the database differs), or it
#               stops: comparing different mixtures would be meaningless.
# WRITES        <out.tsv>, and prints the differences and a one-line count
#               per database. Judges nothing new: no pass mark here.
#
# Usage:  Rscript validation/layer_D/compare_layer_D.R <label_A>=<dir_A> <label_B>=<dir_B> <out.tsv>
#         (each dir is an unpacked layer_D_results.tgz: results/, mixtures/)
#*******************************************************************************

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3 || !all(grepl("=", args[1:2])))
  stop("Usage: compare_layer_D.R <label_A>=<dir_A> <label_B>=<dir_B> <out.tsv>")
script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
runs <- lapply(args[1:2], function(a) list(label = sub("=.*", "", a), dir = sub("^[^=]*=", "", a)))
if (runs[[1]]$label == runs[[2]]$label) stop("The two labels must differ")

options(width = 220)
rt <- function(f) read.delim(f, colClasses = "character", na.strings = "", quote = "")

score_run <- function(r) {
  out <- tempfile(fileext = ".tsv")
  # The scorer exits 1 when any sample fails; that is a result, not an error here.
  system2("Rscript", c(file.path(script_dir, "score_layer_D.R"), shQuote(r$dir), shQuote(out)), stdout = FALSE)
  if (!file.exists(out)) stop(sprintf("score_layer_D.R could not score %s (%s)", r$label, r$dir))
  s <- rt(out)
  m <- rt(file.path(r$dir, "results", "triage_metrics.tsv"))
  s$unclassified <- vapply(s$sample, function(id) {
    v <- m$value[m$sample == id & m$metric == "unclassified_fraction"]
    if (length(v)) v else "n/a"
  }, "")
  db <- file.path(r$dir, "results", "database_info.tsv")
  r$db <- if (file.exists(db)) { d <- rt(db); setNames(d[[2]], d[[1]]) } else NULL
  r$scores <- s
  r
}
runs <- lapply(runs, score_run)
a <- runs[[1]]; b <- runs[[2]]

key <- function(s) paste(s$sample, s$true)
if (!setequal(key(a$scores), key(b$scores)))
  stop("The two runs do not hold the same samples and true fractions: not the same mixtures")

keep <- c("sample", "estimate", "got", "result", "unclassified", "dominant")
x <- merge(a$scores[, c("sample", "metric", "true", "tolerance", "expected")],
           setNames(a$scores[, keep], c("sample", paste0(keep[-1], "_", a$label))), by = "sample")
x <- merge(x, setNames(b$scores[, keep], c("sample", paste0(keep[-1], "_", b$label))), by = "sample")
x <- x[match(a$scores$sample, x$sample), ]          # design order, as in truth.tsv
x$changed <- x[[paste0("result_", a$label)]] != x[[paste0("result_", b$label)]] |
             x[[paste0("got_", a$label)]]    != x[[paste0("got_", b$label)]]
write.table(x, args[3], sep = "\t", quote = FALSE, row.names = FALSE, na = "")

for (r in runs) if (!is.null(r$db))
  cat(sprintf("%-8s database: hash.k2d %s bytes, files dated %s\n", r$label,
              format(as.numeric(r$db[["hash_k2d_bytes"]]), big.mark = ",", scientific = FALSE), r$db[["files_dated"]]))
cat("\n")
show <- c("sample", "true", paste0("estimate_", c(a$label, b$label)), paste0("unclassified_", c(a$label, b$label)),
          paste0("got_", c(a$label, b$label)), paste0("result_", c(a$label, b$label)))
print(x[, show], row.names = FALSE, right = FALSE)
cat("\n")
for (r in runs) {
  res <- r$scores$result
  cat(sprintf("%-8s %d / %d scored samples PASS; limit cases: %s\n", r$label, sum(res == "PASS"),
              sum(res %in% c("PASS", "FAIL")), paste(sprintf("%s=%s", r$scores$sample[!res %in% c("PASS", "FAIL")],
                                                             res[!res %in% c("PASS", "FAIL")]), collapse = ", ")))
}
cat(sprintf("Samples whose result or verdict differs between the databases: %d%s\n", sum(x$changed),
            if (any(x$changed)) paste0(" (", paste(x$sample[x$changed], collapse = ", "), ")") else ""))
