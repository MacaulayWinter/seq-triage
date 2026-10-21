#!/usr/bin/env Rscript
#*******************************************************************************
# triage_report.R: per-sample verdicts and a plain-language summary
#
# WHAT IT DOES  1. Joins the metadata screen (every samplesheet row) with the
#                  per-sample detection and routing state.
#               2. Sets each sample's verdict CEILING from how well the reads
#                  confirm its metadata (CLAUDE.md D2).
#               3. Judges each QC number against assets/metric_guide.tsv,
#                  explains it, and lets the worst result LOWER the verdict.
#                  QC can never raise a verdict (D34).
#               4. Writes <prefix>_summary.tsv, <prefix>_metrics.tsv (one row
#                  per sample and metric) and <prefix>_summary.md.
# GATES         Stops, rather than writing a shorter report, if a samplesheet
#               row has no results, or if a sample was sent to QC but has no
#               QC numbers. Same idea as the author's
#               116_check_rerun_complete.sh: nothing vanishes silently.
#
# Usage:  triage_report.R <metadata_screen.tsv> <samples.tsv> <metric_guide.tsv> <qc_values.tsv> <species_groups.tsv> <out_prefix>
#*******************************************************************************

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) {
  stop("Usage: triage_report.R <metadata_screen.tsv> <samples.tsv> <metric_guide.tsv> <qc_values.tsv> <species_groups.tsv> <out_prefix>")
}
read_tsv <- function(f) read.delim(f, colClasses = "character", na.strings = "", check.names = FALSE, quote = "")
screen  <- read_tsv(args[1])
samples <- read_tsv(args[2])
guide   <- read_tsv(args[3])
qc      <- read_tsv(args[4])
groups  <- read_tsv(args[5])
prefix  <- args[6]
options(scipen = 999)


# 1. Join, with a completeness gate --------------------------------------------

missing <- setdiff(screen$sample, samples$sample)
extra   <- setdiff(samples$sample, screen$sample)
if (length(missing) || length(extra)) {
  stop("Sample sets do not match between the samplesheet and the per-sample results.",
       if (length(missing)) paste0("\n  Missing results for: ", paste(missing, collapse = ", ")),
       if (length(extra))   paste0("\n  Results for samples not in the sheet: ", paste(extra, collapse = ", ")))
}
d <- merge(screen, samples, by = "sample", all = TRUE, sort = FALSE,
           suffixes = c("", ".pipeline"))
d <- d[match(screen$sample, d$sample), ]   # keep samplesheet order
if (!"declared_molecule" %in% names(d)) d$declared_molecule <- NA_character_


# 2. Verdict ceiling (D2) --------------------------------------------------------
# WHAT IT DOES  The best verdict a sample can reach, given how well its reads
#               confirm its metadata. QC metrics may only lower it.
# RULE          No PASS unless the reads independently confirm the metadata.
ceiling_for <- function(route, routed_by, confidence, declared, detected, reasons) {
  if (route == "unreadable") return(c("FAIL", sprintf("File could not be read, so nothing else was run: %s. Re-download or re-copy the file.", reasons)))
  if (route == "empty")      return(c("FAIL",      "No reads in the file."))
  if (routed_by == "blocked")
    return(c("UNCERTAIN", "Long (nanopore-type) reads were given as two files. Nanopore reads are never paired, so the samplesheet is probably wrong. No QC was run."))
  if (routed_by == "unsupported") {
    name <- c(pacbio = "PacBio", dnbseq = "MGI/DNBSEQ", element = "Element AVITI")[detected]
    return(c("UNCERTAIN", sprintf("Reads come from %s%s. seq-triage v0.1 has no QC for this platform, so no QC was run.",
                                  if (is.na(name)) detected else name,
                                  if (!is.na(declared)) sprintf(" (declared: %s)", declared) else "")))
  }
  if (routed_by == "none")   return(c("UNCERTAIN", "Platform could not be determined from the reads, and none was declared. No QC was run."))
  if (routed_by == "declared") {
    if (confidence == "uncertain")
      return(c("WARN", sprintf("Reads were inconclusive; QC used the declared platform (%s), which is unconfirmed.", declared)))
    return(c("WARN", sprintf("Declared platform (%s) is not confirmed by the reads (detection confidence: %s).", declared, confidence)))
  }
  # routed_by == "detected": only high-confidence detection routes this way.
  if (is.na(declared))      return(c("WARN", sprintf("Platform inferred from the reads (%s); none was declared.", detected)))
  if (declared == detected) return(c("PASS", sprintf("Reads confirm the declared platform (%s).", declared)))
  c("WARN", sprintf("Metadata says %s but the reads indicate %s. Review this sample.", declared, detected))
}
v <- t(mapply(ceiling_for, d$route, d$routed_by, d$confidence, d$declared_platform, d$detected_platform, d$reasons))
d$ceiling        <- v[, 1]
d$verdict_reason <- v[, 2]

# Molecule check: reads written with U are RNA (direct-RNA basecalls). If the
# sheet says dna, the metadata contradicts the reads, so cap at WARN. The
# reverse (rna declared, reads written with T) proves nothing: cDNA is
# written with T.
u_dna <- grepl("U_BASES", ifelse(is.na(d$detection_flags), "", d$detection_flags)) &
         d$declared_molecule %in% "dna"
d$ceiling[u_dna & d$ceiling == "PASS"] <- "WARN"
d$verdict_reason[u_dna] <- paste(d$verdict_reason[u_dna],
  "Molecule is declared dna, but the reads are written with U, as direct RNA is. Review this sample.")


# 3. QC metrics (D34) --------------------------------------------------------------

# WHAT IT DOES  The numeric value of one source field for one sample.
#               seqkit reports one row per file, so its fields are summed
#               over the sample's files (R1 + R2). Every other tool writes one
#               file per sample. "a / b" in the guide means a ratio.
field_value <- function(rows, field, tool) {
  if (grepl(" / ", field, fixed = TRUE)) {
    parts <- strsplit(field, " / ", fixed = TRUE)[[1]]
    num <- field_value(rows, parts[1], tool)
    den <- field_value(rows, parts[2], tool)
    return(if (is.na(num) || is.na(den) || den == 0) NA_real_ else num / den)
  }
  x <- suppressWarnings(as.numeric(rows$value[rows$field == field]))
  if (!length(x) || anyNA(x)) return(NA_real_)
  if (tool == "seqkit") return(sum(x))
  if (length(x) > 1) stop("More than one value of ", field, " from ", tool, " for one sample")
  x
}

# FastQC writes one status per file; report the worst of R1 and R2.
status_value <- function(rows, field) {
  s <- rows$value[rows$field == field]
  if (!length(s)) return(NA_character_)
  c("PASS", "WARN", "FAIL")[max(match(s, c("PASS", "WARN", "FAIL")), na.rm = TRUE)]
}

QUALITY_METRICS <- c("q30_rate", "mean_read_quality")   # meaningless on placeholder quality (D26)

# Contamination screen (D41) ----------------------------------------------------
# On only if Kraken2 ran (a database was given); then every sample sent to QC
# must have Kraken2 results.
SCREEN_ON   <- any(qc$tool == "Kraken2")
SCREEN_TOOL <- "Kraken2"
HOST_TAXID  <- "9606"           # Homo sapiens
MIN_SPECIES_READS <- 10         # below this, no species-level estimate (as Bracken)
# A species counts towards off-target only at >= 1% of species-level reads
# (D47). Validation layer D showed pure cultures carry a long tail of trace
# misclassifications (8-53 species, each < 1%, 2.4-3.1% in total) that no
# real contaminant resembles. The tail is reported, not hidden.
MIN_SPECIES_SHARE <- 0.01

# WHAT IT DOES  Contamination numbers for one sample from its Kraken2 report
#               (host, unclassified) and species-level counts (Bracken for
#               Illumina when it ran, Kraken2's own species counts otherwise).
#               Look-alike species in species_groups.tsv count as one.
# RETURNS       named list: host_fraction, unclassified_fraction,
#               off_target_fraction, dominant_species, top_species, and
#               species_note when no species-level estimate was possible.
contamination <- function(rows, route) {
  k <- rows[rows$tool == "Kraken2", ]
  part <- do.call(rbind, strsplit(k$field, "|", fixed = TRUE))
  k$rank <- part[, 1]; k$taxid <- part[, 2]; k$name <- part[, 3]
  k$n <- as.numeric(k$value)
  total <- sum(k$n[k$taxid %in% c("0", "1")])          # unclassified + root
  out <- list(host_fraction = if (total > 0) sum(k$n[k$taxid == HOST_TAXID]) / total else NA_real_,
              unclassified_fraction = if (total > 0) sum(k$n[k$taxid == "0"]) / total else NA_real_,
              off_target_fraction = NA_real_, dominant_species = NA_character_, top_species = NA_character_)

  b <- rows[rows$tool == "Bracken", ]
  if (route == "illumina" && nrow(b)) {
    bp <- do.call(rbind, strsplit(b$field, "|", fixed = TRUE))
    sp <- data.frame(rank = bp[, 1], taxid = bp[, 2], name = bp[, 3], n = as.numeric(b$value))
    source <- "Bracken"
  } else {
    sp <- k[k$rank == "S", c("rank", "taxid", "name", "n")]
    source <- "Kraken2"
  }
  sp <- sp[sp$taxid != HOST_TAXID & sp$rank == "S", ]
  if (sum(sp$n) < MIN_SPECIES_READS) {
    out$species_note <- sprintf("species-level estimate not possible: fewer than %d reads were classified to a species", MIN_SPECIES_READS)
    return(out)
  }
  sp$unit <- ifelse(sp$taxid %in% groups$taxid, groups$group[match(sp$taxid, groups$taxid)], sp$name)
  by_unit <- sort(tapply(sp$n, sp$unit, sum), decreasing = TRUE)
  share <- by_unit / sum(by_unit)
  others <- share[-1]
  counted <- others[others >= MIN_SPECIES_SHARE]
  trace   <- others[others <  MIN_SPECIES_SHARE]
  out$off_target_fraction <- sum(counted)
  out$dominant_species <- sprintf("%s (%.1f%% of species-level reads, %s)", names(share)[1], 100 * share[[1]], source)
  top <- head(share, 3)
  out$top_species <- paste(c(sprintf("%s %.1f%%", names(top), 100 * top),
                             if (length(trace)) sprintf("%d trace species under 1%% each, %.1f%% in total (not counted as off-target)",
                                                        length(trace), 100 * sum(trace))),
                           collapse = "; ")
  out
}

# WHAT IT DOES  Judges every guide metric that applies to one sample.
# RETURNS       data.frame(sample, metric, value, unit, warn, fail, result, note)
#               result: PASS | WARN | FAIL | info | exempt | not_meaningful | not_available
judge_sample <- function(s) {
  route <- s$route
  if (!route %in% c("illumina", "ont")) return(NULL)
  tag  <- paste0(route, "_reads")
  rows <- qc[qc$sample == s$sample, ]
  main_tool <- if (route == "illumina") "fastp" else "NanoPlot"
  if (!any(rows$tool == main_tool)) {
    stop(sprintf("Sample %s was sent to %s QC but no %s results reached the report.", s$sample, route, main_tool))
  }
  declared  <- na.omit(c(s$declared_type, s$declared_molecule))
  flags     <- if (is.na(s$detection_flags)) character(0) else strsplit(s$detection_flags, ",")[[1]]
  contam    <- NULL
  if (SCREEN_ON) {
    if (!any(rows$tool == "Kraken2")) {
      stop(sprintf("Sample %s was sent to QC and the contamination screen is on, but no Kraken2 results reached the report.", s$sample))
    }
    contam <- contamination(rows, route)
  }

  applies <- which(vapply(strsplit(guide$applies_to, ","), function(a) tag %in% a, TRUE) &
                   (SCREEN_ON | guide$source_tool != SCREEN_TOOL))
  out <- lapply(applies, function(j) {
    g      <- guide[j, ]
    trows  <- rows[rows$tool == g$source_tool, ]
    exempt <- intersect(strsplit(ifelse(is.na(g$exempt_types), "", g$exempt_types), ",")[[1]], declared)
    warn   <- suppressWarnings(as.numeric(g$warn))
    fail   <- suppressWarnings(as.numeric(g$fail))
    note   <- ""

    if (g$source_tool == SCREEN_TOOL) {
      value <- contam[[g$metric]]
      if (g$direction == "label") {
        return(data.frame(metric = g$metric, value = if (is.na(value)) "" else value,
                          result = if (is.na(value)) "not_available" else "info",
                          note = if (is.na(value) && !is.null(contam$species_note)) contam$species_note else ""))
      }
    }

    if (g$direction == "status") {
      value  <- status_value(trows, g$source_field)
      result <- if (is.na(value)) "not_available" else "info"
      if (length(exempt) && !is.na(value)) note <- sprintf("known false alarm for %s libraries", exempt[1])
      return(data.frame(metric = g$metric, value = value, result = result, note = note))
    }

    if (g$source_tool != SCREEN_TOOL) value <- field_value(trows, g$source_field, g$source_tool)
    result <-
      if (is.na(value))                                 "not_available"
      else if (length(exempt))                          "exempt"
      else if (g$metric %in% QUALITY_METRICS &&
               "QUALITY_UNINFORMATIVE" %in% flags)      "not_meaningful"
      else if (g$direction == "higher_is_better") {
        if (!is.na(fail) && value < fail) "FAIL" else if (!is.na(warn) && value < warn) "WARN" else "PASS"
      } else if (g$direction == "lower_is_better") {
        if (!is.na(fail) && value > fail) "FAIL" else if (!is.na(warn) && value > warn) "WARN" else "PASS"
      } else "info"
    note <- switch(result,
      not_available  = if (g$source_tool == SCREEN_TOOL && !is.null(contam$species_note)) contam$species_note
                       else sprintf("not reported by %s", g$source_tool),
      exempt         = sprintf("not judged: %s library", exempt[1]),
      not_meaningful = "not judged: quality scores are a placeholder (every base the same)",
      "")
    data.frame(metric = g$metric, value = format_value(value, g$unit), result = result, note = note)
  })
  m <- do.call(rbind, out)
  m <- cbind(sample = s$sample, m)
  gi <- match(m$metric, guide$metric)
  m$unit <- guide$unit[gi]
  m$warn <- vapply(gi, function(j) format_value(guide$warn[j], guide$unit[j]), "")
  m$fail <- vapply(gi, function(j) format_value(guide$fail[j], guide$unit[j]), "")
  m[, c("sample", "metric", "value", "unit", "warn", "fail", "result", "note")]
}

# WHAT IT DOES  Formats a number for people: 87.5%, 49,986, Q16.6.
format_value <- function(x, unit) {
  x <- suppressWarnings(as.numeric(x))
  if (is.na(x)) return("")
  switch(unit,
    fraction = sprintf("%s%%", formatC(100 * x, format = "f", digits = if (100 * x < 1) 2 else 1)),
    Phred    = sprintf("Q%.1f", x),
    reads    = ,
    bp       = formatC(x, format = "f", digits = 0, big.mark = ","),
    format(x))
}

metrics <- do.call(rbind, lapply(seq_len(nrow(d)), function(i) judge_sample(d[i, ])))
if (is.null(metrics)) {
  metrics <- data.frame(sample = character(0), metric = character(0), value = character(0), unit = character(0),
                        warn = character(0), fail = character(0), result = character(0), note = character(0))
}


# 4. Final verdict: the worse of the ceiling and the worst metric ----------------

RANK <- c(PASS = 1, WARN = 2, FAIL = 3)
label <- setNames(guide$metric, guide$metric)
label[] <- gsub("_", " ", label)
qc_line <- function(m) {
  limit <- if (m$result == "FAIL") m$fail else m$warn
  dir   <- guide$direction[guide$metric == m$metric]
  sprintf("%s is %s (%s limit: %s %s).", label[[m$metric]], m$value, tolower(m$result),
          if (dir == "higher_is_better") "below" else "above", limit)
}

d$qc_result   <- ""
d$qc_findings <- ""
d$verdict     <- d$ceiling
for (i in seq_len(nrow(d))) {
  m <- metrics[metrics$sample == d$sample[i] & metrics$result %in% c("WARN", "FAIL"), ]
  if (!nrow(m)) {
    if (d$route[i] %in% c("illumina", "ont")) d$qc_result[i] <- "PASS"
    next
  }
  worst <- names(RANK)[max(RANK[m$result])]
  d$qc_result[i]   <- worst
  d$qc_findings[i] <- paste(vapply(seq_len(nrow(m)), function(k) qc_line(m[k, ]), ""), collapse = " ")
  if (d$ceiling[i] %in% names(RANK) && RANK[[worst]] > RANK[[d$ceiling[i]]]) d$verdict[i] <- worst
}

# Quality that could not be judged (placeholder scores, D26) means QC is
# incomplete, so the sample cannot PASS (D1; user decision after validation
# round 1, finding C1).
for (i in seq_len(nrow(d))) {
  if (d$verdict[i] != "PASS") next
  if (any(metrics$sample == d$sample[i] & metrics$result == "not_meaningful")) {
    d$verdict[i]     <- "WARN"
    d$qc_result[i]   <- "WARN"
    d$qc_findings[i] <- trimws(paste(d$qc_findings[i],
      "Quality could not be assessed: every base has the same score. The sample cannot PASS until real quality values are available."))
  }
}

# Hint: a metric failed that would NOT be judged for some library type the
# sample has not declared (e.g. short reads in an amplicon library). The data
# may be fine; the metadata is what is missing. Only types whose column is
# blank are suggested: sample_type for amplicon/metagenome, molecule for rna.
TYPE_COLUMN <- c(amplicon = "sample_type", metagenome = "sample_type", isolate = "sample_type",
                 unknown = "sample_type", rna = "molecule", dna = "molecule")
d$hint <- ""
for (i in seq_len(nrow(d))) {
  m <- metrics[metrics$sample == d$sample[i] & metrics$result %in% c("WARN", "FAIL"), ]
  if (!nrow(m)) next
  ex <- guide$exempt_types[match(m$metric, guide$metric)]
  types <- unique(unlist(strsplit(ex[!is.na(ex)], ",")))
  blank <- c(sample_type = is.na(d$declared_type[i]), molecule = is.na(d$declared_molecule[i]))
  types <- types[blank[TYPE_COLUMN[types]]]
  if (!length(types)) next
  which_metrics <- label[m$metric[!is.na(ex) & vapply(strsplit(ifelse(is.na(ex), "", ex), ","), function(t) any(t %in% types), TRUE)]]
  d$hint[i] <- sprintf("If this is %s data, declare it in the samplesheet (%s). The %s limit%s %s not applied to such libraries.",
                       paste(types, collapse = " or "),
                       paste(unique(sprintf("%s = %s", TYPE_COLUMN[types], types)), collapse = ", "),
                       paste(which_metrics, collapse = " and "),
                       if (length(which_metrics) > 1) "s" else "", if (length(which_metrics) > 1) "are" else "is")
}


# 5. Write -----------------------------------------------------------------------

cols <- c("sample", "verdict", "verdict_reason", "qc_result", "qc_findings", "hint",
          "declared_platform", "declared_platform_raw", "detected_platform", "detected_layout",
          "confidence", "route", "routed_by", "n_files", "n_reads", "declared_type",
          "declared_molecule", "n_flags", "flags", "detection_flags", "reasons")
write.table(d[, cols], paste0(prefix, "_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
write.table(metrics, paste0(prefix, "_metrics.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")

show <- function(x) ifelse(is.na(x) | x == "", "-", x)
SYMBOL <- c(PASS = "pass", WARN = "**warn**", FAIL = "**FAIL**", info = "info", exempt = "not judged",
            not_meaningful = "not judged", not_available = "n/a")

metric_table <- function(s) {
  m <- metrics[metrics$sample == s, ]
  if (!nrow(m)) return(character(0))
  side  <- ifelse(guide$direction[match(m$metric, guide$metric)] == "higher_is_better", "below", "above")
  limit <- ifelse(m$warn == "" & m$fail == "", "-",
           ifelse(m$fail == "", sprintf("warn %s %s", side, m$warn),
                  sprintf("warn %s %s; fail %s %s", side, m$warn, side, m$fail)))
  c("", "| Metric | Value | Limits | Result | Note |", "|---|---|---|---|---|",
    sprintf("| %s | %s | %s | %s | %s |", label[m$metric], show(m$value), limit,
            SYMBOL[m$result], show(m$note)))
}

# "How to read these numbers": only the metrics that appear in this report.
present <- unique(metrics$metric[metrics$result != "not_available"])
present <- guide$metric[guide$metric %in% present]   # guide order
guide_section <- unlist(lapply(present, function(mt) {
  g <- guide[guide$metric == mt, ]
  limits <- if (g$direction %in% c("higher_is_better", "lower_is_better"))
    sprintf("- **Limits:** %s%s. %s",
            if (!is.na(g$warn)) paste0("warn ", if (g$direction == "higher_is_better") "below " else "above ", format_value(g$warn, g$unit)) else "",
            if (!is.na(g$fail)) paste0(", fail ", if (g$direction == "higher_is_better") "below " else "above ", format_value(g$fail, g$unit)) else "",
            if (!is.na(g$exempt_types)) sprintf("Not applied to: %s.", gsub(",", ", ", g$exempt_types)) else "")
  else if (g$direction == "status") "- **Information only:** never changes the verdict."
  else if (g$source_tool == SCREEN_TOOL) "- **Information only:** never changes the verdict."
  else "- **Information only:** no limit, because judging it needs to know the organism."
  c(sprintf("### %s", label[[mt]]), "",
    sprintf("- **What it measures:** %s", g$what_it_measures),
    sprintf("- **How to read it:** %s", g$how_to_read_it),
    limits,
    if (!is.na(g$caveat)) sprintf("- **Caveat:** %s", g$caveat),
    sprintf("- **Source of the limits:** %s", g$provenance), "")
}))

counts <- table(factor(d$verdict, c("PASS", "WARN", "UNCERTAIN", "FAIL")))
md <- c(
  "# seq-triage summary", "",
  sprintf("%d sample(s): %s.", nrow(d), paste(sprintf("%d %s", counts, names(counts)), collapse = ", ")), "",
  "A sample can only PASS when its reads confirm its metadata and every QC number is within its limits.",
  if (SCREEN_ON) "The contamination screen ran (Kraken2): host and contaminant fractions are judged below." else
    "The contamination screen did not run (no Kraken2 database was given): host DNA and contaminants were not checked.",
  "QC numbers can lower a verdict but never raise it. The limits are provisional defaults; see",
  "*How to read these numbers* at the end.", "",
  "| Sample | Verdict | Declared | Detected | Confidence | QC | Metadata flags |",
  "|---|---|---|---|---|---|---|",
  sprintf("| %s | **%s** | %s | %s | %s | %s | %s |", d$sample, d$verdict, show(d$declared_platform),
          show(d$detected_platform), show(d$confidence), show(d$qc_result), show(d$n_flags)), "",
  "## Per-sample detail", "",
  unlist(lapply(seq_len(nrow(d)), function(i) c(
    sprintf("### %s: %s", d$sample[i], d$verdict[i]), "",
    paste0("- ", d$verdict_reason[i]),
    if (nzchar(d$qc_findings[i])) paste0("- QC: ", d$qc_findings[i]),
    if (nzchar(d$hint[i])) paste0("- Hint: ", d$hint[i]),
    paste0("- Detection (", show(d$confidence[i]), " confidence): ", show(d$reasons[i])),
    if (!SCREEN_ON && d$route[i] %in% c("illumina", "ont")) "- Contamination: not screened (no Kraken2 database was given).",
    if (!is.na(d$detection_flags[i])) paste0("- Detection flags: ", gsub(",", ", ", d$detection_flags[i], fixed = TRUE)),
    if (!is.na(d$flags[i])) paste0("- Metadata: ", strsplit(d$flags[i], " | ", fixed = TRUE)[[1]]),
    metric_table(d$sample[i]),
    ""))),
  if (length(guide_section)) c("## How to read these numbers", "", guide_section)
)
writeLines(md, paste0(prefix, "_summary.md"))
message(sprintf("Report written for %d sample(s), %d QC value(s) judged.", nrow(d), nrow(metrics)))
