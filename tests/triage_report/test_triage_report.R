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

# Kraken2 report lines as the workflow flattens them: "rank|taxid|name" = clade reads.
# A clean E. coli isolate: 1,000,000 reads, 2% unclassified, a little Shigella.
KRAKEN_OK <- c("U|0|unclassified" = 20000, "R|1|root" = 980000, "S|562|Escherichia coli" = 900000,
               "S|623|Shigella flexneri" = 30000, "S|573|Klebsiella pneumoniae" = 5000, "S|9606|Homo sapiens" = 2000)

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
       declared = NA, verdict = "UNCERTAIN"),
  # D45: a damaged file is a FAIL row, not a stopped run.
  list(id = "UNREADABLE",        route = "unreadable", routed_by = "none", confidence = "uncertain",
       flags = "FILE_UNREADABLE", verdict = "FAIL"),
  # Decision C1 (validation round 1): headers confirm ONT, but quality is a
  # placeholder, so it cannot be assessed: WARN, never PASS.
  list(id = "ONT_PLACEHOLDER_HIGH", route = "ont", nano = c(mean_qual = 40), flags = "QUALITY_UNINFORMATIVE",
       verdict = "WARN", check = c(mean_read_quality = "not_meaningful")),
  # Contamination screen (D41). Shigella counts with E. coli (species_groups.tsv).
  list(id = "SCREEN_CLEAN",       route = "illumina", screen = TRUE, verdict = "PASS",
       check = c(host_fraction = "PASS", off_target_fraction = "PASS", dominant_species = "info")),
  list(id = "SCREEN_HOST_WARN",   route = "ont", screen = TRUE, kraken = c("S|9606|Homo sapiens" = 50000),
       type = "isolate", verdict = "WARN", check = c(host_fraction = "WARN")),          # declared: no hint
  list(id = "SCREEN_HOST_FAIL",   route = "ont", screen = TRUE, kraken = c("S|9606|Homo sapiens" = 200000),
       verdict = "FAIL", check = c(host_fraction = "FAIL"), hint = TRUE),              # blank type: metagenome hint
  list(id = "SCREEN_HOST_META",   route = "ont", screen = TRUE, kraken = c("S|9606|Homo sapiens" = 200000), type = "metagenome",
       verdict = "PASS", check = c(host_fraction = "exempt")),
  list(id = "SCREEN_MIXED",       route = "ont", screen = TRUE, kraken = c("S|573|Klebsiella pneumoniae" = 300000),
       verdict = "FAIL", check = c(off_target_fraction = "FAIL"), hint = TRUE),
  list(id = "SCREEN_SAME_GENUS",  route = "ont", screen = TRUE, kraken = c("S|564|Escherichia fergusonii" = 90000),
       verdict = "WARN", check = c(off_target_fraction = "WARN"), hint = TRUE),
  list(id = "SCREEN_BRACKEN",     route = "illumina", screen = TRUE, bracken = c("S|562|Escherichia coli" = 700000, "S|573|Klebsiella pneumoniae" = 300000),
       verdict = "FAIL", check = c(off_target_fraction = "FAIL"), hint = TRUE),
  list(id = "SCREEN_NO_SPECIES",  route = "ont", screen = TRUE, kraken = c("U|0|unclassified" = 999995, "R|1|root" = 5, "S|562|Escherichia coli" = 5, "S|623|Shigella flexneri" = 0, "S|573|Klebsiella pneumoniae" = 0, "S|9606|Homo sapiens" = 0),
       verdict = "PASS", check = c(off_target_fraction = "not_available", unclassified_fraction = "info")),
  # D47: trace species (< 1% each) do not count as off-target, however many;
  # a real contaminant above 1% still does. 40 trace species of 1,000 reads
  # each (~0.1% each, ~4% in total) on a clean culture:
  list(id = "SCREEN_TRACE_TAIL",   route = "ont", screen = TRUE, type = "isolate",
       kraken = setNames(rep(1000, 40), sprintf("S|%d|Trace species %d", 900001:900040, 1:40)),
       verdict = "PASS", check = c(off_target_fraction = "PASS")),
  list(id = "SCREEN_TAIL_PLUS_REAL", route = "ont", screen = TRUE, type = "isolate",
       kraken = c(setNames(rep(1000, 40), sprintf("S|%d|Trace species %d", 900001:900040, 1:40)),
                  "S|573|Klebsiella pneumoniae" = 60000),
       verdict = "WARN", check = c(off_target_fraction = "WARN")),
  # D40: a platform seq-triage names but has no QC for.
  list(id = "UNSUPPORTED_PACBIO", route = "unrouted", routed_by = "unsupported", detected = "pacbio",
       declared = NA, flags = "PLATFORM_NOT_HANDLED", verdict = "UNCERTAIN"),
  # D55: BAM header @RG PL. Blank samplesheet platform, header names it, reads
  # confirm: the header counts as the declaration, so PASS is possible.
  list(id = "BAM_HEADER_DECLARED", route = "illumina", declared = NA, hdr_declared = "illumina",
       declared_from = "bam_header", verdict = "PASS"),
  # Same reads, no platform anywhere: inferred, at most WARN (control).
  list(id = "BAM_NO_PLATFORM",     route = "illumina", declared = NA, verdict = "WARN"),
  # D64: PacBio HiFi QC. Same NanoPlot fields as nanopore, own quality limit
  # (Q20 / Q10): Q15 passes for nanopore but WARNs for HiFi, which shows the
  # report uses the platform's own guide row, not the first row of that name.
  list(id = "HIFI_GOOD",       route = "pacbio", declared = "pacbio", nano = c(mean_qual = 25),
       verdict = "PASS", check = c(mean_read_quality = "PASS")),
  list(id = "HIFI_Q15",        route = "pacbio", declared = "pacbio", nano = c(mean_qual = 15),
       verdict = "WARN", check = c(mean_read_quality = "WARN")),
  list(id = "ONT_Q15",         route = "ont", nano = c(mean_qual = 15),
       verdict = "PASS", check = c(mean_read_quality = "PASS")),
  list(id = "HIFI_INFERRED",   route = "pacbio", declared = NA, nano = c(mean_qual = 25), verdict = "WARN"),
  list(id = "PACBIO_NOT_HIFI", route = "unrouted", routed_by = "pacbio_not_hifi", declared = "pacbio",
       detected = "uncertain", confidence = "uncertain", verdict = "UNCERTAIN"),
  # D60: a medium call contradicting the declaration is not routed.
  list(id = "CONTRADICTED_MEDIUM", route = "unrouted", routed_by = "contradicted", declared = "ont",
       detected = "illumina", confidence = "medium", verdict = "UNCERTAIN"),
  # Samplesheet and header disagree: two metadata sources conflict, no PASS.
  list(id = "BAM_HEADER_DISAGREES", route = "illumina", flags = "BAM_HEADER_DISAGREES", verdict = "WARN")
)

build <- function(cases, dir) {
  screen <- samples <- qc <- list()
  for (cs in cases) {
    route    <- cs$route
    declared <- if ("declared" %in% names(cs)) cs$declared else if (route %in% c("illumina", "ont", "pacbio")) route else NA
    screen[[cs$id]] <- data.frame(sample = cs$id, declared_platform_raw = declared, declared_platform = declared,
      declared_type = if (is.null(cs$type)) NA else cs$type,
      declared_molecule = if (is.null(cs$molecule)) NA else cs$molecule, n_flags = "0", flags = NA)
    samples[[cs$id]] <- data.frame(sample = cs$id, single_end = "false", n_files = "1", n_reads = "1000",
      detected_platform = if (!is.null(cs$detected)) cs$detected else if (route %in% c("illumina", "ont", "pacbio")) route else "uncertain",
      detected_layout = "single", confidence = if (is.null(cs$confidence)) "high" else cs$confidence,
      detection_flags = if (is.null(cs$flags)) NA else cs$flags, reasons = "synthetic",
      route = route, routed_by = if (is.null(cs$routed_by)) "detected" else cs$routed_by,
      declared_platform = if (!is.null(cs$hdr_declared)) cs$hdr_declared else declared,
      declared_from = if (!is.null(cs$declared_from)) cs$declared_from else if (!is.na(declared)) "samplesheet" else NA)
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
    if (route %in% c("ont", "pacbio")) qc[[paste0(cs$id, "n")]] <- add("NanoPlot", NANO_OK, cs$nano)
    if (isTRUE(cs$screen)) {
      kr <- add("Kraken2", KRAKEN_OK, cs$kraken)
      if (!is.null(cs$kraken) && any(names(cs$kraken) != "U|0|unclassified")) {
        # keep root consistent: root = classified total
        v <- as.numeric(kr$value); n <- kr$field
        v[n == "R|1|root"] <- sum(v[grepl("^S\\|", n)])
        kr$value <- format(v, scientific = FALSE, trim = TRUE)
        if ("R|1|root" %in% names(cs$kraken)) kr$value[n == "R|1|root"] <- format(cs$kraken[["R|1|root"]], scientific = FALSE)
      }
      qc[[paste0(cs$id, "k")]] <- kr
      if (!is.null(cs$bracken)) qc[[paste0(cs$id, "b")]] <- add("Bracken", cs$bracken, NULL)
    }
  }
  write_tsv(do.call(rbind, screen),  file.path(dir, "screen.tsv"))
  write_tsv(do.call(rbind, samples), file.path(dir, "samples.tsv"))
  write_tsv(do.call(rbind, qc),      file.path(dir, "qc.tsv"))
}

run_report <- function(dir) {
  system2("Rscript", c("--vanilla", "bin/triage_report.R", file.path(dir, "screen.tsv"), file.path(dir, "samples.tsv"),
                       "assets/metric_guide.tsv", file.path(dir, "qc.tsv"), "assets/species_groups.tsv", file.path(dir, "out")),
          stdout = FALSE, stderr = FALSE)
}

# Two report runs: without the contamination screen, and with it (when the
# screen is on, EVERY sample sent to QC must have Kraken2 results, so the two
# kinds cannot share one report).
screened <- vapply(cases, function(cs) isTRUE(cs$screen), TRUE)
summ <- mets <- NULL
for (part in list(cases[!screened], cases[screened])) {
  d <- tempfile("run_", tmpdir = tmp); dir.create(d)
  build(part, d)
  status <- run_report(d)
  if (status != 0) stop("triage_report.R exited with status ", status)
  summ <- rbind(summ, read.delim(file.path(d, "out_summary.tsv"), colClasses = "character", na.strings = ""))
  mets <- rbind(mets, read.delim(file.path(d, "out_metrics.tsv"), colClasses = "character", na.strings = ""))
}

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

# Gate: screen on, but one sample sent to QC has no Kraken2 results.
gate2 <- file.path(tmp, "gate2"); dir.create(gate2)
build(list(list(id = "HAS_SCREEN", route = "ont", screen = TRUE, verdict = ""),
           list(id = "NO_SCREEN",  route = "ont", verdict = "")), gate2)
gate2_status <- run_report(gate2)
res <- rbind(res, data.frame(case = "GATE_MISSING_KRAKEN2", got = sprintf("exit %d", gate2_status),
                             expected = "non-zero exit", result = if (gate2_status != 0) "PASS" else "FAIL"))

print(res, row.names = FALSE, right = FALSE)
n_fail <- sum(res$result == "FAIL")
cat(sprintf("\n%d / %d cases as expected\n", nrow(res) - n_fail, nrow(res)))
unlink(tmp, recursive = TRUE)
if (n_fail > 0) quit(status = 1)
