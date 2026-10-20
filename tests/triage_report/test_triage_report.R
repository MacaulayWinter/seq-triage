#!/usr/bin/env Rscript
# Test for bin/triage_report.R
#
# Builds SYNTHETIC report inputs (metadata screen, per-sample state, QC
# values) in a temporary folder, one sample per rule, runs the report with
# the real assets/metric_guide.tsv, and checks each sample's verdict and the
# result of the metric the case is about. Also checks that the report
# refuses to run when a sample sent to QC has no QC numbers.
#
# Run from the repository root:
#   Rscript tests/triage_report/test_triage_report.R
# Exits non-zero if any check differs, so CI can use it.

tmp <- tempfile("report_test_")
dir.create(tmp)
write_tsv <- function(x, f) write.table(x, f, sep = "\t", quote = FALSE, row.names = FALSE, na = "")

# Healthy QC values; each case overrides what it is about.
FASTP_OK <- c("summary.after_filtering.total_reads" = 500000, "summary.after_filtering.q30_rate" = 0.9,
              "filtering_result.passed_filter_reads" = 495000, "summary.before_filtering.total_reads" = 500000,
              "duplication.rate" = 0.05, "adapter_cutting.adapter_trimmed_reads" = 1000,
              "summary.after_filtering.gc_content" = 0.5, "insert_size.peak" = 300,
              "summary.after_filtering.read1_mean_length" = 150)
NANO_OK  <- c(number_of_reads = 20000, number_of_bases = 1e8, n50 = 10000, mean_qual = 15, median_read_length = 5000)
SEQKIT_OK <- c(num_seqs = 250000, sum_len = 75000000, sum_n = 1000)

# Cases: platform route, metadata, detection, QC overrides (NA removes a
# field), and what to expect.
cases <- list(
  list(id = "ILL_GOOD",          route = "illumina", verdict = "PASS"),
  list(id = "ILL_LOW_READS",     route = "illumina", fastp = c("summary.after_filtering.total_reads" = 50000),
       verdict = "WARN", check = c(total_reads = "WARN")),
  list(id = "ILL_AT_WARN_LIMIT", route = "illumina", fastp = c("summary.after_filtering.total_reads" = 100000),
       verdict = "PASS", check = c(total_reads = "PASS")),
  list(id = "ILL_Q30_FAIL",      route = "illumina", fastp = c("summary.after_filtering.q30_rate" = 0.6),
       verdict = "FAIL", check = c(q30_rate = "FAIL")),
  list(id = "ILL_DUP_DNA",       route = "illumina", fastp = c("duplication.rate" = 0.6),
       type = "isolate", molecule = "dna",              # fully declared: no hint
       verdict = "WARN", check = c(duplication_rate = "WARN")),
  list(id = "ILL_DUP_UNLABELLED", route = "illumina", fastp = c("duplication.rate" = 0.6),
       verdict = "WARN", check = c(duplication_rate = "WARN"), hint = TRUE),
  list(id = "ILL_DUP_RNA",       route = "illumina", fastp = c("duplication.rate" = 0.6), molecule = "rna",
       verdict = "PASS", check = c(duplication_rate = "exempt")),
  list(id = "ILL_NO_ADAPTER",    route = "illumina", fastp = c("adapter_cutting.adapter_trimmed_reads" = NA, "insert_size.peak" = NA),
       verdict = "PASS", check = c(adapter_trimmed_pct = "not_available", insert_size_peak = "not_available")),
  list(id = "ILL_N_RATE",        route = "illumina", seqkit = c(sum_n = 1500000),
       verdict = "WARN", check = c(n_base_rate = "WARN")),
  list(id = "ILL_FASTQC_INFO",   route = "illumina", fastqc = "FAIL",
       verdict = "PASS", check = c(fastqc_per_base_sequence_content = "info")),
  list(id = "ONT_GOOD",          route = "ont", verdict = "PASS"),
  list(id = "ONT_AMPLICON",      route = "ont", nano = c(n50 = 600), type = "amplicon",
       verdict = "PASS", check = c(read_n50 = "exempt")),
  list(id = "ONT_UNLABELLED",    route = "ont", nano = c(n50 = 600),
       verdict = "FAIL", check = c(read_n50 = "FAIL"), hint = TRUE),
  list(id = "ONT_PLACEHOLDER_Q", route = "ont", nano = c(mean_qual = 40), flags = "QUALITY_UNINFORMATIVE",
       confidence = "medium", routed_by = "declared",
       verdict = "WARN", check = c(mean_read_quality = "not_meaningful")),
  list(id = "CEILING_WARN_QC_FAIL", route = "ont", nano = c(number_of_reads = 500), confidence = "medium",
       routed_by = "declared", verdict = "FAIL", check = c(number_of_reads = "FAIL")),
  list(id = "ONT_U_DECLARED_DNA", route = "ont", flags = "U_BASES", molecule = "dna", verdict = "WARN"),
  list(id = "UNROUTED",          route = "unrouted", routed_by = "none", confidence = "uncertain",
       declared = NA, verdict = "UNCERTAIN")
)

build <- function(cases, dir) {
  screen <- samples <- qc <- list()
  for (cs in cases) {
    route    <- cs$route
    declared <- if ("declared" %in% names(cs)) cs$declared else if (route %in% c("illumina", "ont")) route else NA
    screen[[cs$id]] <- data.frame(sample = cs$id, declared_platform_raw = declared, declared_platform = declared,
      declared_type = if (is.null(cs$type)) NA else cs$type,
      declared_molecule = if (is.null(cs$molecule)) NA else cs$molecule, n_flags = "0", flags = NA)
    samples[[cs$id]] <- data.frame(sample = cs$id, single_end = "false", n_files = "1", n_reads = "1000",
      detected_platform = if (route %in% c("illumina", "ont")) route else "uncertain",
      detected_layout = "single", confidence = if (is.null(cs$confidence)) "high" else cs$confidence,
      detection_flags = if (is.null(cs$flags)) NA else cs$flags, reasons = "synthetic",
      route = route, routed_by = if (is.null(cs$routed_by)) "detected" else cs$routed_by)
    add <- function(tool, base, over) {
      v <- base; if (length(over)) v[names(over)] <- over
      v <- v[!is.na(v)]
      data.frame(sample = cs$id, tool = tool, file = "f", field = names(v), value = format(v, scientific = FALSE, trim = TRUE))
    }
    qc[[paste0(cs$id, "s")]] <- add("seqkit", SEQKIT_OK, cs$seqkit)
    if (route == "illumina") {
      qc[[paste0(cs$id, "f")]] <- add("fastp", FASTP_OK, cs$fastp)
      st <- if (is.null(cs$fastqc)) "PASS" else cs$fastqc
      qc[[paste0(cs$id, "q")]] <- data.frame(sample = cs$id, tool = "FastQC", file = "f",
        field = c("Per base sequence content", "Overrepresented sequences", "Per sequence GC content"), value = st)
    }
    if (route == "ont") qc[[paste0(cs$id, "n")]] <- add("NanoPlot", NANO_OK, cs$nano)
  }
  write_tsv(do.call(rbind, screen),  file.path(dir, "screen.tsv"))
  write_tsv(do.call(rbind, samples), file.path(dir, "samples.tsv"))
  write_tsv(do.call(rbind, qc),      file.path(dir, "qc.tsv"))
}

run_report <- function(dir) {
  system2("Rscript", c("--vanilla", "bin/triage_report.R", file.path(dir, "screen.tsv"), file.path(dir, "samples.tsv"),
                       "assets/metric_guide.tsv", file.path(dir, "qc.tsv"), file.path(dir, "out")),
          stdout = FALSE, stderr = FALSE)
}

# Main run
build(cases, tmp)
status <- run_report(tmp)
if (status != 0) stop("triage_report.R exited with status ", status)
summ <- read.delim(file.path(tmp, "out_summary.tsv"), colClasses = "character", na.strings = "")
mets <- read.delim(file.path(tmp, "out_metrics.tsv"), colClasses = "character", na.strings = "")

res <- do.call(rbind, lapply(cases, function(cs) {
  got_v <- summ$verdict[summ$sample == cs$id]
  got_c <- if (length(cs$check)) paste(names(cs$check), mets$result[match(paste(cs$id, names(cs$check)), paste(mets$sample, mets$metric))], sep = "=", collapse = " ") else ""
  exp_c <- if (length(cs$check)) paste(names(cs$check), cs$check, sep = "=", collapse = " ") else ""
  got_h <- !is.na(summ$hint[summ$sample == cs$id])
  ok <- identical(got_v, cs$verdict) && got_c == exp_c && got_h == isTRUE(cs$hint)
  data.frame(case = cs$id, got = paste(got_v, got_c, if (got_h) "+hint" else ""),
             expected = paste(cs$verdict, exp_c, if (isTRUE(cs$hint)) "+hint" else ""),
             result = if (ok) "PASS" else "FAIL")
}))

# Gate: a sample routed to QC with no QC numbers must stop the report.
gate <- file.path(tmp, "gate"); dir.create(gate)
build(list(list(id = "NO_QC", route = "illumina", verdict = "")), gate)
qcg <- read.delim(file.path(gate, "qc.tsv"), colClasses = "character")
write_tsv(qcg[qcg$tool == "seqkit", ], file.path(gate, "qc.tsv"))
gate_status <- run_report(gate)
res <- rbind(res, data.frame(case = "GATE_MISSING_QC", got = sprintf("exit %d", gate_status),
                             expected = "non-zero exit", result = if (gate_status != 0) "PASS" else "FAIL"))

print(res, row.names = FALSE, right = FALSE)
n_fail <- sum(res$result == "FAIL")
cat(sprintf("\n%d / %d cases as expected\n", nrow(res) - n_fail, nrow(res)))
unlink(tmp, recursive = TRUE)
if (n_fail > 0) quit(status = 1)
