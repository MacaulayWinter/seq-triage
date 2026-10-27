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
if (!length(args) %in% 6:7) {
  stop("Usage: triage_report.R <metadata_screen.tsv> <samples.tsv> <metric_guide.tsv> <qc_values.tsv> <species_groups.tsv> <out_prefix> [species_genome_size.tsv]")
}
read_tsv <- function(f) read.delim(f, colClasses = "character", na.strings = "", check.names = FALSE, quote = "")
screen  <- read_tsv(args[1])
samples <- read_tsv(args[2])
guide   <- read_tsv(args[3])
qc      <- read_tsv(args[4])
groups  <- read_tsv(args[5])
prefix  <- args[6]
# D80: NCBI expected genome sizes per species taxid (optional 7th argument).
gsizes  <- if (length(args) >= 7 && file.exists(args[7])) read_tsv(args[7]) else data.frame(species_taxid = character(0), expected_ungapped_length = character(0))
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
if (!"declared_target" %in% names(d))   d$declared_target   <- NA_character_
if (!"declared_genome_size" %in% names(d)) d$declared_genome_size <- NA_character_
# BAM input (D55): with no platform in the samplesheet, the workflow may take
# one from the BAM header (@RG PL); it is then the declared platform here too.
if (all(c("declared_from", "declared_platform.pipeline") %in% names(d))) {
  from_hdr <- d$declared_from %in% "bam_header"
  d$declared_platform[from_hdr] <- d$declared_platform.pipeline[from_hdr]
}


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
    return(c("UNCERTAIN", sprintf("Reads come from %s%s. seq-triage has no QC for this platform, so no QC was run.",
                                  if (is.na(name)) detected else name,
                                  if (!is.na(declared)) sprintf(" (declared: %s)", declared) else "")))
  }
  if (routed_by == "none")   return(c("UNCERTAIN", "Platform could not be determined from the reads, and none was declared. No QC was run."))
  # D78: only index or barcode reads were given.
  if (routed_by == "index_reads")
    return(c("UNCERTAIN", sprintf("Only index or barcode reads were given (%s). Give the sample's read files (R1, R2) instead; no QC was run.", reasons)))
  # D64: only PacBio HiFi is QC'd; a declaration alone cannot route PacBio.
  if (routed_by == "pacbio_not_hifi")
    return(c("UNCERTAIN", "Declared PacBio, but the reads do not show PacBio HiFi markers (read names ending /ccs, or a BAM written by ccs). seq-triage QCs PacBio HiFi only, so no QC was run."))
  # D60: the reads lean one way (not confidently), the metadata says another.
  if (routed_by == "contradicted")
    return(c("UNCERTAIN", sprintf(paste("Metadata says %s but the reads suggest %s (detection confidence: %s). Neither is",
                                        "certain enough to choose the QC, so no QC was run. Review this sample."),
                                  declared, detected, confidence)))
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
# The samplesheet and the BAM header name different platforms (D55): two
# sources of metadata disagree, so the sample cannot PASS (D1).
hdr_conflict <- grepl("BAM_HEADER_DISAGREES", ifelse(is.na(d$detection_flags), "", d$detection_flags))
lift <- hdr_conflict & d$ceiling == "PASS"
d$ceiling[lift] <- "WARN"
d$verdict_reason[lift] <- paste(d$verdict_reason[lift],
  "But the BAM header names a different platform from the samplesheet. Review this sample.")

# D78: an index read given as fastq_2 means the samplesheet is probably
# wrong (R2 missing, or I1 given in its place): at most WARN, review.
idx_mate <- grepl("INDEX_READ_AS_MATE", ifelse(is.na(d$detection_flags), "", d$detection_flags))
d$ceiling[idx_mate & d$ceiling == "PASS"] <- "WARN"
d$verdict_reason[idx_mate] <- paste(d$verdict_reason[idx_mate],
  "fastq_2 looks like an index read, not R2, so QC used fastq_1 only. Review the samplesheet.")

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
# Routes that get QC. PacBio means HiFi only (D64): other PacBio reads are
# recognised but not routed.
QC_ROUTES <- c("illumina", "dnbseq", "element", "ont", "pacbio")
# Short-read routes (D69): fastp + FastQC, Bracken, read-start amplicon signal.
SHORT_ROUTES <- c("illumina", "dnbseq", "element")
HOST_TAXID  <- "9606"           # Homo sapiens
HOST_GENUS_TAXID <- "9605"      # Homo
BACTERIA_TAXID   <- "2"
# D76 (provisional): with no target declared, a sample whose classified reads
# are at least half human is either a human sample or one swamped by human
# DNA; the reads cannot say which, so it is UNCERTAIN until target_taxid says.
HOST_DOMINANT_SHARE <- 0.50
VIRUS_TAXID <- "10239"          # Viruses (domain)
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
#               D76: with a declared target (NCBI taxid), also
#               target_fraction (reads under the target / classified reads,
#               Kraken2 clade count), target_seen, and off-target measured
#               against the TARGET instead of the most common species; a human
#               target keeps human reads in the composition. Also
#               human_share and bacterial_fraction (of classified reads).
contamination <- function(rows, route, target = NA_character_) {
  k <- rows[rows$tool == "Kraken2", ]
  part <- do.call(rbind, strsplit(k$field, "|", fixed = TRUE))
  k$rank <- part[, 1]; k$taxid <- part[, 2]; k$name <- part[, 3]
  k$n <- as.numeric(k$value)
  total <- sum(k$n[k$taxid %in% c("0", "1")])          # unclassified + root
  classified <- sum(k$n[k$taxid == "1"])
  out <- list(host_fraction = if (total > 0) sum(k$n[k$taxid == HOST_TAXID]) / total else NA_real_,
              unclassified_fraction = if (total > 0) sum(k$n[k$taxid == "0"]) / total else NA_real_,
              # D67: share of CLASSIFIED reads under the Viruses domain (taxid
              # 10239), from Kraken2 for every platform. Recognises viruses by
              # taxonomy, not by species name (names change between versions).
              # By taxid ONLY, whatever the rank code: NCBI moved Viruses to an
              # "acellular root" in 2024, so the Oct 2025 Kraken2 databases
              # report it as R1, not D. Matching rank D gave 0% viral on real
              # SARS-CoV-2 runs (layer E re-run, 2026-10-22).
              viral_fraction = if (classified > 0) sum(k$n[k$taxid == VIRUS_TAXID]) / classified else NA_real_,
              off_target_fraction = NA_real_, dominant_species = NA_character_, top_species = NA_character_)
  clade <- function(id) sum(k$n[k$taxid == id])
  has_target <- !is.na(target) && nzchar(target)
  human_target <- has_target && target %in% c(HOST_TAXID, HOST_GENUS_TAXID)
  out$human_share        <- if (classified > 0) clade(HOST_TAXID) / classified else NA_real_
  out$bacterial_fraction <- if (classified > 0) clade(BACTERIA_TAXID) / classified else NA_real_
  out$target_fraction    <- if (has_target && classified > 0) clade(target) / classified else NA_real_
  out$target_seen        <- has_target && clade(target) > 0
  # D73: the most common genus and its share of genus-level reads (Kraken2
  # clade counts at rank G, human excluded), for every platform. A genome the
  # database knows only through relatives (a new phage, a poorly-known
  # bacterium) is spread over sister species but stays in one genus.
  g <- k[k$rank == "G" & k$taxid != HOST_GENUS_TAXID & k$n > 0, ]
  if (nrow(g) && sum(g$n) > 0) {
    i <- which.max(g$n)
    out$top_genus <- g$name[i]; out$top_genus_share <- g$n[i] / sum(g$n)
  }

  b <- rows[rows$tool == "Bracken", ]
  if (route %in% SHORT_ROUTES && nrow(b)) {
    bp <- do.call(rbind, strsplit(b$field, "|", fixed = TRUE))
    sp <- data.frame(rank = bp[, 1], taxid = bp[, 2], name = bp[, 3], n = as.numeric(b$value))
    source <- "Bracken"
  } else {
    sp <- k[k$rank == "S", c("rank", "taxid", "name", "n")]
    source <- "Kraken2"
  }
  sp <- sp[(human_target | sp$taxid != HOST_TAXID) & sp$rank == "S", ]
  if (sum(sp$n) < MIN_SPECIES_READS) {
    out$species_note <- sprintf("species-level estimate not possible: fewer than %d reads were classified to a species", MIN_SPECIES_READS)
    return(out)
  }
  sp$unit <- ifelse(sp$taxid %in% groups$taxid, groups$group[match(sp$taxid, groups$taxid)], sp$name)
  by_unit <- sort(tapply(sp$n, sp$unit, sum), decreasing = TRUE)
  share <- by_unit / sum(by_unit)
  # D76: off-target is everything that is not the declared target, when the
  # target is a species the screen saw; otherwise the most common species.
  target_unit <- if (out$target_seen) unique(sp$unit[sp$taxid == target]) else character(0)
  others <- if (length(target_unit) == 1) share[names(share) != target_unit] else share[-1]
  out$target_not_dominant <- length(target_unit) == 1 && names(share)[1] != target_unit
  counted <- others[others >= MIN_SPECIES_SHARE]
  trace   <- others[others <  MIN_SPECIES_SHARE]
  out$off_target_fraction <- sum(counted)
  if (has_target && !out$target_seen) {
    out$off_target_fraction <- NA_real_
    out$species_note <- sprintf("no reads were classified to the declared target (taxid %s): the database may not contain it (Standard-8 holds no eukaryote but human), or the reads are not from it, so contamination cannot be judged against the target", target)
  }
  # For the sample-type call (D65): the top unit's share, and how many units
  # (the top one included) reach the 1% floor.
  out$dominant_share <- share[[1]]
  # D80: the top unit's taxid (the most-read species in a look-alike group).
  top_rows <- sp[sp$unit == names(share)[1], ]
  out$dominant_taxid <- top_rows$taxid[which.max(top_rows$n)]
  out$n_species_1pct <- 1 + length(counted)
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
  if (!route %in% QC_ROUTES) return(NULL)
  tag  <- paste0(route, "_reads")
  rows <- qc[qc$sample == s$sample, ]
  main_tool <- if (route %in% SHORT_ROUTES) "fastp" else "NanoPlot"
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
    contam <- contamination(rows, route, s$declared_target)
  }
  human_target <- s$declared_target %in% c(HOST_TAXID, HOST_GENUS_TAXID)

  applies <- which(vapply(strsplit(guide$applies_to, ","), function(a) tag %in% a, TRUE) &
                   (SCREEN_ON | guide$source_tool != SCREEN_TOOL) & guide$source_tool != "computed" &
                   # D81: the rRNA share only for RNA samples (or where SortMeRNA ran).
                   (guide$source_tool != "SortMeRNA" | any(rows$tool == "SortMeRNA") | s$declared_molecule %in% "rna"))
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
                          note = if (is.na(value) && !is.null(contam$species_note)) contam$species_note else "", gi = j))
      }
    }

    if (g$direction == "status") {
      value  <- status_value(trows, g$source_field)
      result <- if (is.na(value)) "not_available" else "info"
      if (length(exempt) && !is.na(value)) note <- sprintf("known false alarm for %s libraries", exempt[1])
      return(data.frame(metric = g$metric, value = value, result = result, note = note, gi = j))
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
    if (g$metric == "host_fraction" && human_target && result %in% c("PASS", "WARN", "FAIL")) result <- "exempt"
    # D78: bisulfite-converted reads (C read as T) classify poorly, so the
    # screen's contamination numbers are not judged.
    bisulfite <- "BISULFITE" %in% flags && g$metric %in% c("host_fraction", "off_target_fraction")
    if (bisulfite && result %in% c("PASS", "WARN", "FAIL")) result <- "exempt"
    note <- switch(result,
      not_available  = if (g$source_tool == "SortMeRNA") "not screened: give --rrna_db (an rRNA reference) to measure it"
                       else if (g$metric == "target_fraction") (if (is.na(s$declared_target)) "no target_taxid declared"
                                                           else "not computed: no reads were classified at all")
                       else if (g$source_tool == SCREEN_TOOL && !is.null(contam$species_note)) contam$species_note
                       else sprintf("not reported by %s", g$source_tool),
      exempt         = if (bisulfite) "not judged: bisulfite-converted reads cannot be classified reliably"
                       else if (g$metric == "host_fraction" && human_target) "not judged: human is the declared target"
                       else sprintf("not judged: %s library", exempt[1]),
      not_meaningful = "not judged: quality scores are a placeholder (every base the same)",
      "")
    data.frame(metric = g$metric, value = format_value(value, g$unit), result = result, note = note, gi = j)
  })
  m <- do.call(rbind, out)
  m <- cbind(sample = s$sample, m)
  # gi: the guide row each metric was judged by. A metric can have one row
  # per platform (e.g. mean_read_quality for nanopore and for PacBio HiFi,
  # D64), so rows are never looked up by metric name.
  gi <- m$gi
  m$unit <- guide$unit[gi]
  m$warn <- vapply(gi, function(j) format_value(guide$warn[j], guide$unit[j]), "")
  m$fail <- vapply(gi, function(j) format_value(guide$fail[j], guide$unit[j]), "")
  m[, c("sample", "metric", "value", "unit", "warn", "fail", "result", "note", "gi")]
}

# WHAT IT DOES  Formats a number for people: 87.5%, 49,986, Q16.6.
format_value <- function(x, unit) {
  x <- suppressWarnings(as.numeric(x))
  if (is.na(x)) return("")
  switch(unit,
    fraction = sprintf("%s%%", formatC(100 * x, format = "f", digits = if (100 * x < 1) 2 else 1)),
    Phred    = sprintf("Q%.1f", x),
    depth    = if (x < 1) sprintf("%.2fx", x) else sprintf("%.1fx", x),
    reads    = ,
    bp       = formatC(x, format = "f", digits = 0, big.mark = ","),
    format(x))
}

metrics <- do.call(rbind, lapply(seq_len(nrow(d)), function(i) judge_sample(d[i, ])))
if (is.null(metrics)) {
  metrics <- data.frame(sample = character(0), metric = character(0), value = character(0), unit = character(0),
                        warn = character(0), fail = character(0), result = character(0), note = character(0),
                        gi = integer(0))
}


# 3b. Sample type (D65) ------------------------------------------------------------
# WHAT IT DOES  Calls isolate, metagenome or amplicon from the reads, with a
#               confidence, or says uncertain. It never relaxes a verdict:
#               only a DECLARED type stops a metric being judged (D27). A
#               high-confidence call that contradicts the declared type caps
#               the sample at WARN (review), and an undeclared amplicon or
#               metagenome strengthens the hint to declare it.
# SIGNALS       Amplicon (Illumina only): share of reads starting with one of
#               the 10 most common 20-base starts (detection). Calibrated on 18
#               runs (validation round 1 + test data): amplicons 0.26-0.35,
#               whole-genome and metagenome <= 0.034, metatranscriptomes
#               0.05-0.13 (rRNA). Nanopore and PacBio reads start with
#               adapter sequence, so amplicons cannot be told from them.
#               Isolate vs metagenome: the contamination screen's species
#               counts (look-alike groups combined, human excluded, 1% floor).
#               Provisional limits from layer D pure cultures (top species
#               97-99%) and mixtures (30% contamination -> ~70%).
for (col in c("start_share", "probe_reads", "gc_excess", "kmer_repeated", "kmer_bases", "kmer_genome_mb")) if (!col %in% names(d)) d[[col]] <- NA_character_
TYPE_MIN_PROBE      <- 1000   # probed reads needed to judge read starts
TYPE_START_AMPLICON <- 0.20
TYPE_START_RANDOM   <- 0.05
TYPE_ISOLATE_TOP    <- 0.90   # and at most 2 other species at >= 1%
TYPE_META_TOP       <- 0.50   # or at least 5 species at >= 1%
# D67 (decided after validation layer E, so a later round must test them):
TYPE_ISOLATE_MAX_UNCLASSIFIED <- 0.25  # isolate HIGH needs most reads classified
                                       # (layer E: Illumina isolates <= 11.3%, metagenomes median 93.7%)
TYPE_VIRAL_DOMINANT <- 0.50            # virus-dominated: isolate at most medium (viral
                                       # genomes are usually amplified; fragmented amplicons
                                       # do not start at primers: layer E SARS-CoV-2)
# D73 (decided after layer E round 2, so a later round must test it):
TYPE_META_MAX_GENUS <- 0.50
# D86 (thresholds chosen on layer E rounds 1-2, tested on round 3; a later
# round must confirm): a metagenome call is high only when the reads
# themselves, with no database, look like many genomes.
TYPE_META_MAX_KMER   <- 0.49    # share of sampled 21-mers seen twice: below this
TYPE_META_MIN_GC     <- 0.054   # excess per-read GC spread: at least this
TYPE_META_MIN_KBASES <- 5e6     # bases examined for the k-mer check: at least this            # metagenome HIGH needs the top genus under 50%
                                       # (round 2: metagenomes <= 37%; phages and poorly-known
                                       # isolates split over sister species 72-100%)
classify_type <- function(s) {
  out <- function(type, conf, why, amplicon_judged) list(type = type, conf = conf, why = why, amplicon_judged = amplicon_judged)
  if (!s$route %in% QC_ROUTES) return(out("uncertain", "uncertain", "not assessed: the sample was not sent to QC.", FALSE))
  # D78: library kinds the isolate / metagenome / amplicon question does not fit.
  lib_flags <- if (is.na(s$detection_flags)) character(0) else strsplit(s$detection_flags, ",")[[1]]
  if ("BARCODE_READ" %in% lib_flags)
    return(out("uncertain", "uncertain", "A single-cell library (cell-barcode R1): isolate, metagenome or amplicon does not describe it.", FALSE))
  if ("BISULFITE" %in% lib_flags)
    return(out("uncertain", "uncertain", "A bisulfite / EM-seq methylation library: converted reads cannot show the composition reliably.", FALSE))
  short <- s$route %in% SHORT_ROUTES
  amp <- "unknown"; amp_note <- "Nanopore and PacBio reads start with adapter sequence, so an amplicon library cannot be recognised from them."
  if (short) {
    st <- suppressWarnings(as.numeric(s$start_share)); n <- suppressWarnings(as.numeric(s$probe_reads))
    if (is.na(st) || is.na(n) || n < TYPE_MIN_PROBE) {
      amp_note <- sprintf("Too few reads (%s) to judge whether reads start at primers.", if (is.na(n)) "?" else format(n, big.mark = ","))
    } else {
      pct <- sprintf("%.0f%% of reads start with one of 10 sequences", 100 * st)
      # D85: RNA libraries start unevenly (a few highly expressed transcripts
      # give many reads with the same start), so for declared RNA the read
      # starts cannot PROVE an amplicon: at most medium.
      if (st >= TYPE_START_AMPLICON && s$declared_molecule %in% "rna")
        return(out("amplicon", "medium", sprintf("%s, as amplicon reads do; but RNA libraries also start unevenly (highly expressed transcripts), so this is not proof of an amplicon.", pct), TRUE))
      if (st >= TYPE_START_AMPLICON)
        return(out("amplicon", "high", sprintf("%s: amplicon reads start at their primers.", pct), TRUE))
      if (st >= TYPE_START_RANDOM)
        return(out("uncertain", "uncertain", sprintf("%s: more than fragmented DNA gives (under %.0f%%), fewer than amplicons (%.0f%% or more). rRNA-rich RNA libraries look like this.",
                                                     pct, 100 * TYPE_START_RANDOM, 100 * TYPE_START_AMPLICON), TRUE))
      amp <- "no"; amp_note <- sprintf("%s, as fragmented DNA gives, so not an amplicon library.", sub("^(\\d+)%", "Only \\1%", pct))
    }
  }
  if (!SCREEN_ON) return(out("uncertain", "uncertain", paste(amp_note, "Isolate or metagenome needs the contamination screen (--kraken2_db)."), amp == "no"))
  co <- contamination(qc[qc$sample == s$sample, ], s$route, s$declared_target)
  # D76: composition cannot speak for a target the database does not see, or
  # for a human-dominated sample nobody has said is human.
  if (!is.na(s$declared_target) && !co$target_seen)
    return(out("uncertain", "uncertain", paste(amp_note, sprintf("No reads were classified to the declared target (taxid %s), so the composition the database sees says nothing about the sample.", s$declared_target)), amp == "no"))
  if (is.na(s$declared_target) && !s$declared_type %in% "metagenome" && !is.na(co$human_share) && co$human_share >= HOST_DOMINANT_SHARE)
    return(out("uncertain", "uncertain", paste(amp_note, sprintf("%.0f%% of classified reads are human and no target_taxid is declared: a human sample, or one swamped by human DNA.", 100 * co$human_share)), amp == "no"))
  if (is.null(co$dominant_share))
    return(out("uncertain", "uncertain", paste(amp_note, "Too few reads were classified to a species to tell isolate from metagenome."), amp == "no"))
  top <- co$dominant_share; k <- co$n_species_1pct
  comp <- sprintf("the most common species holds %.0f%% of species-level reads and %d species %s 1%%", 100 * top, k,
                  if (k == 1) "reaches" else "reach")
  if (top >= TYPE_ISOLATE_TOP && k - 1 <= 2) {
    unc <- co$unclassified_fraction; vir <- co$viral_fraction
    if (!is.na(vir) && vir >= TYPE_VIRAL_DOMINANT)
      return(out("isolate", "medium", sprintf("One organism, a virus (%.0f%% of classified reads are viral): %s. Viral genomes are usually sequenced from amplicons, which do not start at primers when the amplicons were fragmented, so an amplicon cannot be excluded.", 100 * vir, comp), FALSE))
    if (!is.na(unc) && unc > TYPE_ISOLATE_MAX_UNCLASSIFIED)
      return(out("isolate", "medium", sprintf("%s, but %.0f%% of reads are unclassified, so the composition rests on the few reads the database recognises (a poorly-known community also looks like this). %s", sub("^t", "T", comp), 100 * unc, amp_note), amp == "no"))
    if (amp == "no") return(out("isolate", "high", sprintf("One organism: %s. %s", comp, amp_note), TRUE))
    return(out("isolate", "medium", sprintf("One organism: %s. %s It could also be a single-organism amplicon.", comp, amp_note), FALSE))
  }
  if (top < TYPE_META_TOP || k >= 5) {
    # High only when the reads also rule out an amplicon (Illumina): a 16S
    # amplicon of a community looks like a metagenome by composition.
    if (amp == "no") {
      gs <- co$top_genus_share
      if (is.null(gs) || is.na(gs))
        return(out("metagenome", "medium", sprintf("Many species: %s, but no genus-level counts, so one organism split over related species cannot be excluded. %s", comp, amp_note), TRUE))
      if (gs >= TYPE_META_MAX_GENUS)
        return(out("metagenome", "medium", sprintf("Many species: %s, but one genus (%s) holds %.0f%% of genus-level reads: one organism the database knows only through related species (e.g. a new phage) looks like this. %s",
                                                   comp, co$top_genus, 100 * gs, amp_note), TRUE))
      kr <- suppressWarnings(as.numeric(s$kmer_repeated)); kb <- suppressWarnings(as.numeric(s$kmer_bases))
      gx <- suppressWarnings(as.numeric(s$gc_excess))
      if (is.na(kr) || is.na(kb) || is.na(gx) || kb < TYPE_META_MIN_KBASES)
        return(out("metagenome", "medium", sprintf("Many species: %s, but too little sequence (%s bases; %s needed) to check from the reads alone that they come from many genomes. %s",
                                                   comp, if (is.na(kb)) "?" else format(kb, big.mark = ","), format(TYPE_META_MIN_KBASES, big.mark = ",", scientific = FALSE), amp_note), TRUE))
      if (kr >= TYPE_META_MAX_KMER)
        return(out("metagenome", "medium", sprintf("Many species: %s, but %.0f%% of sampled k-mers recur in the reads, as one genome's reads do (a community's rarely do): an organism the database knows only through relatives looks like this. %s",
                                                   comp, 100 * kr, amp_note), TRUE))
      if (gx < TYPE_META_MIN_GC)
        return(out("metagenome", "medium", sprintf("Many species: %s, but the reads share a narrow GC range (excess spread %.3f), as one genome's reads do. %s",
                                                   comp, gx, amp_note), TRUE))
      return(out("metagenome", "high", sprintf("Many organisms: %s; the most common genus holds %.0f%%, and the reads themselves look like many genomes (%.0f%% of k-mers recur, GC spread %.3f). %s",
                                               comp, 100 * gs, 100 * kr, gx, amp_note), TRUE))
    }
    return(out("metagenome", "medium", sprintf("Many organisms: %s. %s It could also be an amplicon (e.g. 16S) of a community.", comp, amp_note), FALSE))
  }
  out("uncertain", "uncertain", sprintf("Between isolate and metagenome: %s (a contaminated isolate also looks like this). %s", comp, amp_note), amp == "no")
}
ty <- lapply(seq_len(nrow(d)), function(i) classify_type(d[i, ]))
d$detected_type   <- vapply(ty, `[[`, "", "type")
d$type_confidence <- vapply(ty, `[[`, "", "conf")
d$type_reason     <- vapply(ty, `[[`, "", "why")
amp_judged        <- vapply(ty, `[[`, TRUE, "amplicon_judged")
# A declared amplicon is never contradicted when the reads cannot show
# whether they are an amplicon (long reads, too few reads): a nanopore 16S
# amplicon of a community would otherwise look like a "metagenome".
declared_t <- ifelse(d$declared_type %in% c("isolate", "metagenome", "amplicon"), d$declared_type, NA)
type_conflict <- !is.na(declared_t) & d$type_confidence == "high" & d$detected_type != declared_t &
                 !(declared_t == "amplicon" & !amp_judged)
d$detection_flags[type_conflict] <- ifelse(is.na(d$detection_flags[type_conflict]), "TYPE_CONTRADICTS",
                                           paste0(d$detection_flags[type_conflict], ",TYPE_CONTRADICTS"))
lift <- type_conflict & d$ceiling == "PASS"
d$ceiling[lift] <- "WARN"
d$verdict_reason[type_conflict] <- paste(d$verdict_reason[type_conflict],
  sprintf("The samplesheet says %s, but the reads look like %s (%s) Review this sample.",
          declared_t[type_conflict], d$detected_type[type_conflict], d$type_reason[type_conflict]))


# 3b'. Genome size and depth of coverage (D80) -------------------------------------
# WHAT IT DOES  For a single organism's DNA (isolate declared or detected, or a
#               target declared): genome size from the samplesheet
#               (genome_size), else NCBI's expected size for the declared
#               target_taxid, else for the most common species; depth =
#               all sequenced bases (every file) / genome size. Information
#               only: thin data is not bad data (it feeds the advice, D82).
#               Not computed for communities, amplicons, RNA or single-cell
#               libraries, where one genome's depth means nothing.
gsize_of <- function(taxid) {
  v <- suppressWarnings(as.numeric(gsizes$expected_ungapped_length[match(taxid, gsizes$species_taxid)]))
  if (length(v) && !is.na(v)) v else NA_real_
}
depth_rows <- list()
gi_size  <- match("genome_size", guide$metric); gi_depth <- match("depth_of_coverage", guide$metric)
if (!is.na(gi_size) && !is.na(gi_depth)) for (i in seq_len(nrow(d))) {
  s <- d[i, ]
  if (!s$route %in% QC_ROUTES) next
  lf <- if (is.na(s$detection_flags)) character(0) else strsplit(s$detection_flags, ",")[[1]]
  one_organism <- s$declared_type %in% "isolate" || !is.na(s$declared_target) || s$detected_type == "isolate"
  why_not <- if (s$declared_molecule %in% "rna") "not computed: RNA (expression, not genome coverage)"
             else if ("BARCODE_READ" %in% lf) "not computed: single-cell library"
             else if (s$declared_type %in% c("metagenome", "amplicon") || s$detected_type %in% c("metagenome", "amplicon")) "not computed: not a single genome"
             else if (!one_organism) "not computed: not known to be one organism"
             else ""
  size <- NA_real_; src <- ""
  if (!nzchar(why_not)) {
    if (!is.na(s$declared_genome_size)) { size <- as.numeric(s$declared_genome_size); src <- "declared genome_size" }
    if (is.na(size) && !is.na(s$declared_target)) { size <- gsize_of(s$declared_target); src <- sprintf("NCBI expected size for the declared target (taxid %s)", s$declared_target) }
    if (is.na(size) && SCREEN_ON) {
      co <- contamination(qc[qc$sample == s$sample, ], s$route, s$declared_target)
      if (!is.null(co$dominant_taxid)) { size <- gsize_of(co$dominant_taxid); src <- sprintf("NCBI expected size for %s", sub(" \\(.*$", "", co$dominant_species)) }
    }
    if (is.na(size)) why_not <- if (SCREEN_ON) "no genome size: declare genome_size or target_taxid (the screen found no species to look one up for)"
                                else "no genome size: declare genome_size or target_taxid (or run the contamination screen)"
  }
  bases <- sum(suppressWarnings(as.numeric(qc$value[qc$sample == s$sample & qc$tool == "seqkit" & qc$field == "sum_len"])), na.rm = TRUE)
  depth_rows[[length(depth_rows) + 1]] <- data.frame(
    sample = s$sample, metric = c("genome_size", "depth_of_coverage"),
    value = if (nzchar(why_not)) c("", "") else c(format_value(size, "bp"), format_value(bases / size, "depth")),
    unit = c("bp", "depth"), warn = "", fail = "",
    result = if (nzchar(why_not)) "not_available" else "info",
    note = if (nzchar(why_not)) why_not else c(src, sprintf("%s bases / %s", format_value(bases, "bp"), format_value(size, "bp"))),
    gi = c(gi_size, gi_depth))
}
if (length(depth_rows)) metrics <- rbind(metrics, do.call(rbind, depth_rows))

# D86: the two one-genome-or-many measures, as information.
gi_k <- match("kmer_repeat_share", guide$metric); gi_g <- match("gc_spread", guide$metric)
if (!is.na(gi_k) && !is.na(gi_g)) for (i in seq_len(nrow(d))) {
  if (!d$route[i] %in% QC_ROUTES) next
  kr <- suppressWarnings(as.numeric(d$kmer_repeated[i])); gx <- suppressWarnings(as.numeric(d$gc_excess[i]))
  kb <- suppressWarnings(as.numeric(d$kmer_bases[i]))
  metrics <- rbind(metrics, data.frame(sample = d$sample[i], metric = c("kmer_repeat_share", "gc_spread"),
    value = c(if (is.na(kr)) "" else format_value(kr, "fraction"), if (is.na(gx)) "" else format(round(gx, 3))),
    unit = c("fraction", "number"), warn = "", fail = "", result = ifelse(is.na(c(kr, gx)), "not_available", "info"),
    note = c(if (is.na(kb)) "" else sprintf("of sampled 21-mers in %s bases", format_value(kb, "bp")), "excess SD of per-read GC beyond chance"),
    gi = c(gi_k, gi_g)))
}


# 3c. Target organism (D76) --------------------------------------------------------
# WHAT IT DOES  With a declared target_taxid: a target the screen did not see
#               (contamination cannot be judged) or one that is not the most
#               common species caps the sample at WARN. With none declared, a
#               human-dominated sample (not declared a metagenome) is not
#               judged on host_fraction and cannot PASS or WARN: UNCERTAIN,
#               declare the target (a FAIL from another metric still stands).
add_flag_d <- function(i, f) d$detection_flags[i] <<- if (is.na(d$detection_flags[i])) f else paste0(d$detection_flags[i], ",", f)
host_dominant <- rep(FALSE, nrow(d))
if (SCREEN_ON) for (i in seq_len(nrow(d))) {
  if (!d$route[i] %in% QC_ROUTES) next
  co <- contamination(qc[qc$sample == d$sample[i], ], d$route[i], d$declared_target[i])
  t  <- d$declared_target[i]
  if (!is.na(t) && !co$target_seen) {
    add_flag_d(i, "TARGET_NOT_SEEN")
    if (d$ceiling[i] == "PASS") d$ceiling[i] <- "WARN"
    d$verdict_reason[i] <- paste(d$verdict_reason[i], sprintf(
      "No reads were classified to the declared target (taxid %s), so contamination cannot be judged: use a Kraken2 database that contains it (e.g. PlusPF for fungi and protists).", t))
  } else if (!is.na(t) && isTRUE(co$target_not_dominant)) {
    add_flag_d(i, "TARGET_NOT_DOMINANT")
    if (d$ceiling[i] == "PASS") d$ceiling[i] <- "WARN"
    d$verdict_reason[i] <- paste(d$verdict_reason[i], sprintf(
      "The declared target (taxid %s) is not the most common species; the most common is %s. Review this sample.", t, co$dominant_species))
  } else if (is.na(t) && !d$declared_type[i] %in% "metagenome" && !is.na(co$human_share) && co$human_share >= HOST_DOMINANT_SHARE) {
    add_flag_d(i, "HOST_DOMINANT")
    host_dominant[i] <- TRUE
    h <- metrics$sample == d$sample[i] & metrics$metric == "host_fraction" & metrics$result %in% c("PASS", "WARN", "FAIL")
    metrics$result[h] <- "exempt"
    metrics$note[h]   <- "not judged: human-dominated sample with no target_taxid declared"
    d$verdict_reason[i] <- paste(d$verdict_reason[i], sprintf(
      "%.0f%% of classified reads are human and no target_taxid is declared: either a human sample or one heavily contaminated with human DNA. Declare target_taxid (9606 if it is human) to say which.", 100 * co$human_share))
  }
}


# 4. Final verdict: the worse of the ceiling and the worst metric ----------------

RANK <- c(PASS = 1, WARN = 2, FAIL = 3)
label <- setNames(guide$metric, guide$metric)
label[] <- gsub("_", " ", label)
qc_line <- function(m) {
  limit <- if (m$result == "FAIL") m$fail else m$warn
  dir   <- guide$direction[m$gi]
  sprintf("%s is %s (%s limit: %s %s).", label[[m$metric]], m$value, tolower(m$result),
          if (dir == "higher_is_better") "below" else "above", limit)
}

d$qc_result   <- ""
d$qc_findings <- ""
d$verdict     <- d$ceiling
for (i in seq_len(nrow(d))) {
  m <- metrics[metrics$sample == d$sample[i] & metrics$result %in% c("WARN", "FAIL"), ]
  if (!nrow(m)) {
    if (d$route[i] %in% QC_ROUTES) d$qc_result[i] <- "PASS"
    next
  }
  worst <- names(RANK)[max(RANK[m$result])]
  d$qc_result[i]   <- worst
  d$qc_findings[i] <- paste(vapply(seq_len(nrow(m)), function(k) qc_line(m[k, ]), ""), collapse = " ")
  if (d$ceiling[i] %in% names(RANK) && RANK[[worst]] > RANK[[d$ceiling[i]]]) d$verdict[i] <- worst
}

# D76: a human-dominated sample with no declared target is UNCERTAIN, unless
# another metric already FAILs it.
d$verdict[host_dominant & d$verdict %in% c("PASS", "WARN")] <- "UNCERTAIN"

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
  ex <- guide$exempt_types[m$gi]
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

# D65: an undeclared amplicon or metagenome, recognised with high confidence,
# gets a hint even when no metric failed: the declaration changes judgements.
for (i in seq_len(nrow(d))) {
  if (!is.na(d$declared_type[i]) || d$type_confidence[i] != "high" || !d$detected_type[i] %in% c("amplicon", "metagenome")) next
  extra <- sprintf("The reads look like a%s %s library (%s) If so, declare sample_type = %s.",
                   if (d$detected_type[i] == "amplicon") "n" else "", d$detected_type[i], d$type_reason[i], d$detected_type[i])
  d$hint[i] <- trimws(paste(d$hint[i], extra))
}


# 5. Write -----------------------------------------------------------------------

cols <- c("sample", "verdict", "verdict_reason", "qc_result", "qc_findings", "hint",
          "declared_platform", "declared_platform_raw", "detected_platform", "detected_layout",
          "confidence", "route", "routed_by", "n_files", "n_reads", "declared_type",
          "declared_molecule", "n_flags", "flags", "detection_flags", "reasons",
          "detected_type", "type_confidence", "type_reason", "declared_target")
write.table(d[, cols], paste0(prefix, "_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
write.table(metrics[, setdiff(names(metrics), "gi")], paste0(prefix, "_metrics.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")

show <- function(x) ifelse(is.na(x) | x == "", "-", x)
SYMBOL <- c(PASS = "pass", WARN = "**warn**", FAIL = "**FAIL**", info = "info", exempt = "not judged",
            not_meaningful = "not judged", not_available = "n/a")

metric_table <- function(s) {
  m <- metrics[metrics$sample == s, ]
  if (!nrow(m)) return(character(0))
  side  <- ifelse(guide$direction[m$gi] == "higher_is_better", "below", "above")
  limit <- ifelse(m$warn == "" & m$fail == "", "-",
           ifelse(m$fail == "", sprintf("warn %s %s", side, m$warn),
                  sprintf("warn %s %s; fail %s %s", side, m$warn, side, m$fail)))
  c("", "| Metric | Value | Limits | Result | Note |", "|---|---|---|---|---|",
    sprintf("| %s | %s | %s | %s | %s |", label[m$metric], show(m$value), limit,
            SYMBOL[m$result], show(m$note)))
}

# "How to read these numbers": only the metrics that appear in this report.
# One entry per guide ROW used, in guide order; a metric with rows for two
# platforms gets one entry per platform, each with its own limits.
present <- sort(unique(metrics$gi[metrics$result != "not_available"]))
PLATFORM_LABEL <- c(illumina_reads = "Illumina", dnbseq_reads = "MGI/DNBSEQ", element_reads = "Element AVITI",
                    ont_reads = "nanopore", pacbio_reads = "PacBio HiFi")
guide_section <- unlist(lapply(present, function(j) {
  g  <- guide[j, ]
  mt <- g$metric
  if (sum(guide$metric == mt) > 1) {
    plats <- PLATFORM_LABEL[strsplit(g$applies_to, ",")[[1]]]
    mt_title <- sprintf("%s (%s)", label[[mt]], paste(plats[!is.na(plats)], collapse = ", "))
  } else mt_title <- label[[mt]]
  limits <- if (g$direction %in% c("higher_is_better", "lower_is_better"))
    sprintf("- **Limits:** %s%s. %s",
            if (!is.na(g$warn)) paste0("warn ", if (g$direction == "higher_is_better") "below " else "above ", format_value(g$warn, g$unit)) else "",
            if (!is.na(g$fail)) paste0(", fail ", if (g$direction == "higher_is_better") "below " else "above ", format_value(g$fail, g$unit)) else "",
            if (!is.na(g$exempt_types)) sprintf("Not applied to: %s.", gsub(",", ", ", g$exempt_types)) else "")
  else if (g$direction == "status") "- **Information only:** never changes the verdict."
  else if (g$source_tool == SCREEN_TOOL) "- **Information only:** never changes the verdict."
  else "- **Information only:** no limit, because judging it needs to know the organism."
  c(sprintf("### %s", mt_title), "",
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
  if (any(qc$tool == "SortMeRNA")) "The rRNA screen ran (SortMeRNA) on RNA samples: their rRNA share is judged below.",
  "QC numbers can lower a verdict but never raise it. The limits are provisional defaults; see",
  "*How to read these numbers* at the end.", "",
  "**What next:** `recommendations/analyses.md` lists, for each sample, the analyses its data supports",
  "(ready, possible or not now, with the reason) and the nf-core pipeline to run.", "",
  "| Sample | Verdict | Declared | Detected | Confidence | Type | QC | Metadata flags |",
  "|---|---|---|---|---|---|---|---|",
  sprintf("| %s | **%s** | %s | %s | %s | %s | %s | %s |", d$sample, d$verdict, show(d$declared_platform),
          show(d$detected_platform), show(d$confidence),
          # The type used for the judgement: declared wins; otherwise the
          # call from the reads, with its confidence.
          ifelse(d$declared_type %in% c("isolate", "metagenome", "amplicon"), sprintf("%s (declared)", d$declared_type),
                 ifelse(d$detected_type == "uncertain", "-", sprintf("%s (%s)", d$detected_type, d$type_confidence))),
          show(d$qc_result), show(d$n_flags)), "",
  "## Per-sample detail", "",
  unlist(lapply(seq_len(nrow(d)), function(i) c(
    sprintf("### %s: %s", d$sample[i], d$verdict[i]), "",
    paste0("- ", d$verdict_reason[i]),
    if (nzchar(d$qc_findings[i])) paste0("- QC: ", d$qc_findings[i]),
    if (nzchar(d$hint[i])) paste0("- Hint: ", d$hint[i]),
    paste0("- Detection (", show(d$confidence[i]), " confidence): ", show(d$reasons[i])),
    if (!SCREEN_ON && d$route[i] %in% QC_ROUTES) "- Contamination: not screened (no Kraken2 database was given).",
    if (!is.na(d$detection_flags[i])) paste0("- Detection flags: ", gsub(",", ", ", d$detection_flags[i], fixed = TRUE)),
    sprintf("- Sample type (from the reads): %s%s. %s", d$detected_type[i],
            if (d$type_confidence[i] %in% c("high", "medium")) sprintf(", %s confidence", d$type_confidence[i]) else "",
            d$type_reason[i]),
    if (!is.na(d$flags[i])) paste0("- Metadata: ", strsplit(d$flags[i], " | ", fixed = TRUE)[[1]]),
    metric_table(d$sample[i]),
    ""))),
  if (length(guide_section)) c("## How to read these numbers", "", guide_section)
)
writeLines(md, paste0(prefix, "_summary.md"))
message(sprintf("Report written for %d sample(s), %d QC value(s) judged.", nrow(d), nrow(metrics)))
