#!/usr/bin/env Rscript
# Tests for bin/combine_detection.R (D75): per-lane detection rows -> one
# sample row. Run from the repository root:
#   Rscript tests/combine_detection/test_combine_detection.R
lane <- function(platform = "illumina", conf = "high", layout = "paired", flags = "", reads = 50000,
                 files = 2, start = "0.02", probe = "2000", why = "Short, uniform reads with Illumina read names.",
                 gcx = "0.05", krep = "0.50", kb = "10000000") {
  data.frame(sample = "x", n_files = files, n_reads = reads, detected_platform = platform, detected_layout = layout,
             confidence = conf, detection_flags = flags, reasons = why, max_len = "151;151", length_spread = "0;0",
             header_evidence = "illumina;illumina", quality_values = "4;4", u_bases = "0;0", adapters = "",
             start_share = start, probe_reads = probe, gc_excess = gcx, kmer_repeated = krep, kmer_bases = kb,
             kmer_genome_mb = "2.0", stringsAsFactors = FALSE)
}
run <- function(lanes) {
  dir <- tempfile(); dir.create(dir)
  files <- vapply(seq_along(lanes), function(i) {
    f <- file.path(dir, sprintf("S.u%d.detection.tsv", i))
    write.table(lanes[[i]], f, sep = "\t", quote = FALSE, row.names = FALSE, na = ""); f }, "")
  out <- file.path(dir, "S.detection.tsv")
  status <- system2("Rscript", c("--vanilla", "bin/combine_detection.R", "S", out, files))
  if (status != 0) stop("combine_detection.R failed")
  read.delim(out, colClasses = "character", na.strings = "", quote = "")
}
has_flag <- function(r, f) f %in% strsplit(ifelse(is.na(r$detection_flags), "", r$detection_flags), ",")[[1]]

cases <- list(
  list(id = "AGREE", lanes = list(lane(), lane(start = "0.04", probe = "6000")),
       check = function(r) r$detected_platform == "illumina" && r$confidence == "high" && r$detected_layout == "paired" &&
         r$n_files == "4" && r$n_reads == "100000" && abs(as.numeric(r$start_share) - 0.035) < 1e-9 && r$probe_reads == "8000" &&
         grepl("every lane says illumina", r$reasons) && !has_flag(r, "FILES_DISAGREE")),
  list(id = "WEAKEST_LANE", lanes = list(lane(), lane(conf = "medium", why = "Short, uniform reads; no read names.")),
       check = function(r) r$detected_platform == "illumina" && r$confidence == "medium" && grepl("weakest lane", r$reasons) &&
         grepl("Lane 2: Short, uniform reads; no read names.", r$reasons, fixed = TRUE)),
  list(id = "PLATFORMS_DIFFER", lanes = list(lane(), lane(platform = "ont", layout = "single", files = 1)),
       check = function(r) r$detected_platform == "uncertain" && r$confidence == "uncertain" && has_flag(r, "FILES_DISAGREE") &&
         grepl("lane 2: ont", r$reasons) && r$detected_layout == "uncertain"),
  list(id = "ONE_LANE_UNKNOWN", lanes = list(lane(), lane(platform = "uncertain", conf = "uncertain", flags = "FEW_READS")),
       check = function(r) r$detected_platform == "uncertain" && has_flag(r, "FILES_DISAGREE") && has_flag(r, "FEW_READS")),
  list(id = "BROKEN_PAIR", lanes = list(lane(), lane(layout = "uncertain", conf = "low", flags = "MATE_COUNT_MISMATCH")),
       check = function(r) r$detected_platform == "illumina" && r$confidence == "low" && r$detected_layout == "uncertain" &&
         has_flag(r, "MATE_COUNT_MISMATCH") && grepl("do not agree on the layout", r$reasons)),
  # D86: the one-genome measures are averaged over lanes by bases examined.
  list(id = "ONE_GENOME_MEASURES", lanes = list(lane(gcx = "0.04", krep = "0.60", kb = "10000000"), lane(gcx = "0.10", krep = "0.30", kb = "5000000")),
       check = function(r) abs(as.numeric(r$kmer_repeated) - 0.5) < 1e-9 && abs(as.numeric(r$gc_excess) - 0.06) < 1e-9 && r$kmer_bases == "15000000"),
  list(id = "THREE_SINGLE_FILES", lanes = list(lane("ont", layout = "single", files = 1, start = "", probe = ""),
                                               lane("ont", layout = "single", files = 1, start = "", probe = ""),
                                               lane("ont", layout = "single", files = 1, start = "", probe = "")),
       check = function(r) r$detected_platform == "ont" && r$n_files == "3" && r$detected_layout == "single" && is.na(r$start_share) &&
         grepl("Detected per file (3 files)", r$reasons, fixed = TRUE))
)
ok <- vapply(cases, function(cs) isTRUE(cs$check(run(cs$lanes))), TRUE)
for (i in seq_along(cases)) cat(sprintf(" %-20s %s\n", cases[[i]]$id, if (ok[i]) "PASS" else "FAIL"))
cat(sprintf("\n%d / %d cases as expected\n", sum(ok), length(ok)))
if (!all(ok)) quit(status = 1)
