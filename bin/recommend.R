#!/usr/bin/env Rscript
#*******************************************************************************
# recommend.R: which nf-core pipeline each ready sample should go to next,
#              with a ready-to-run samplesheet and params file (CLAUDE.md D66)
#
# WHAT IT DOES  For every sample in the triage summary, either names ONE
#               downstream pipeline (from assets/downstream_pipelines.tsv,
#               pinned versions) or says why there is none. Then writes, per
#               recommended pipeline, a samplesheet in that pipeline's own
#               format, a params.yaml and a README with the command.
# READY MEANS   verdict PASS or WARN; platform confirmed (routed by detection,
#               or by a declaration the reads did not contradict); sample type
#               declared, or detected with high confidence (D65); no metadata
#               conflict flag; FASTQ input (no pipeline here takes BAM).
# FILL_IN       Values only the user knows (primers, databases, reference) are
#               never guessed. Where the pipeline would reject the word
#               FILL_IN (a file that must exist, a number, a fixed list), it
#               is written as `key: FILL_IN`. Free-text fields marked with !
#               in the table (ampliseq's primers) would ACCEPT it and run with
#               nonsense, so they are written commented out instead: the
#               pipeline then stops because the value is missing (checked in
#               ampliseq 2.18.0's code).
# WRITES        <outdir>/recommendations.tsv and <outdir>/<pipeline>/...
#
# Usage:  recommend.R <triage_summary.tsv> <triage_metrics.tsv> <samplesheet.csv> <downstream_pipelines.tsv> <outdir>
#*******************************************************************************

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) stop("Usage: recommend.R <triage_summary.tsv> <triage_metrics.tsv> <samplesheet.csv> <downstream_pipelines.tsv> <outdir>")
rt   <- function(f) read.delim(f, colClasses = "character", na.strings = "", quote = "", check.names = FALSE)
summ <- rt(args[1]); mets <- rt(args[2]); rules <- rt(args[4]); out <- args[5]
sheet <- read.csv(args[3], colClasses = "character", na.strings = "", check.names = FALSE, strip.white = TRUE)
for (col in c("fastq_1", "fastq_2", "bam")) if (!col %in% names(sheet)) sheet[[col]] <- NA_character_
sheet$sample <- trimws(sheet$sample)
dir.create(out, recursive = TRUE, showWarnings = FALSE)

SARSCOV2 <- "Severe acute respiratory syndrome coronavirus 2"
has <- function(x) !is.na(x) & nzchar(x)
in_list <- function(value, list) value %in% strsplit(list, ",")[[1]]

# WHAT IT DOES  Decides for one sample. RETURNS list(rule = row of rules or
#               NULL, reason = text).
decide <- function(s, row) {
  none <- function(why) list(rule = NULL, reason = why)
  if (!s$verdict %in% c("PASS", "WARN"))
    return(none(sprintf("Verdict is %s: resolve the triage findings first.", s$verdict)))
  if (has(row$bam))
    return(none("BAM input: the downstream pipelines take FASTQ. Convert first (samtools collate | samtools fastq)."))
  flags <- if (has(s$detection_flags)) strsplit(s$detection_flags, ",")[[1]] else character(0)
  conflict <- intersect(flags, c("TYPE_CONTRADICTS", "BAM_HEADER_DISAGREES"))
  if (length(conflict))
    return(none(sprintf("Metadata conflict (%s): review the sample first.", paste(conflict, collapse = ", "))))
  platform <- switch(s$routed_by, detected = s$detected_platform, declared = s$declared_platform, NA)
  if (!has(platform)) return(none("Platform not confirmed by the reads."))
  type <- if (has(s$declared_type) && s$declared_type %in% c("isolate", "metagenome", "amplicon")) s$declared_type
          else if (has(s$detected_type) && s$type_confidence %in% "high") s$detected_type else NA
  if (!has(type)) return(none("Sample type not declared and not detected with high confidence: declare sample_type."))
  molecule <- if (has(s$declared_molecule)) s$declared_molecule else "dna"
  dom <- mets$value[mets$sample == s$sample & mets$metric == "dominant_species"]
  target <- if (!length(dom) || !has(dom)) NA else if (startsWith(dom, SARSCOV2)) "sarscov2" else "other"
  fits <- which(vapply(seq_len(nrow(rules)), function(i) {
    r <- rules[i, ]
    in_list(platform, r$platform) && r$sample_type == type && in_list(molecule, r$molecule) &&
      (r$target == "any" || (r$target == "sarscov2" && target %in% "sarscov2") ||
       (r$target == "not_sarscov2" && target %in% "other"))
  }, TRUE))
  if (length(fits)) {
    note <- c(if (!has(s$declared_molecule)) "molecule not declared, so DNA assumed" else NULL,
              if (!has(s$declared_type)) sprintf("sample type %s detected from the reads", type) else NULL)
    return(list(rule = rules[fits[1], ], reason = paste(c(rules$note[fits[1]], if (length(note)) sprintf("(%s)", paste(note, collapse = "; "))), collapse = " ")))
  }
  # Known gaps, said plainly.
  if (type == "amplicon" && platform != "illumina")
    return(none("Long-read amplicons: no pipeline set up here (nf-core/viralrecon expects the sequencer's barcode folders)."))
  if (type == "amplicon" && is.na(target))
    return(none("Amplicon target unknown without the contamination screen (--kraken2_db): marker gene or virus?"))
  if (type == "isolate" && platform == "pacbio")
    return(none("PacBio HiFi isolates: no pipeline set up here yet."))
  if (molecule == "rna" && type != "isolate")
    return(none(sprintf("RNA %s (e.g. metatranscriptome): no pipeline set up here yet.", type)))
  none(sprintf("No pipeline in assets/downstream_pipelines.tsv for %s / %s / %s.", platform, type, molecule))
}

rows <- lapply(seq_len(nrow(summ)), function(i) {
  s   <- summ[i, ]
  row <- sheet[match(s$sample, sheet$sample), ]
  d   <- decide(s, row)
  list(sample = s$sample, row = row, platform = switch(s$routed_by, detected = s$detected_platform, declared = s$declared_platform, NA),
       rule = d$rule, reason = d$reason)
})

rec <- data.frame(
  sample   = vapply(rows, `[[`, "", "sample"),
  pipeline = vapply(rows, function(r) if (is.null(r$rule)) "none" else paste0("nf-core/", r$rule$pipeline), ""),
  version  = vapply(rows, function(r) if (is.null(r$rule)) "" else r$rule$version, ""),
  rule     = vapply(rows, function(r) if (is.null(r$rule)) "" else r$rule$rule, ""),
  reason   = vapply(rows, `[[`, "", "reason"),
  stringsAsFactors = FALSE)
write.table(rec, file.path(out, "recommendations.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")

# Samplesheets, one format per pipeline (checked against each pipeline's own
# assets/schema_input.json at the pinned commit; tests/recommend/).
NAv <- function(x) ifelse(has(x), x, "NA")
blank <- function(x) ifelse(has(x), x, "")
sheet_row <- function(fmt, r) {
  f1 <- r$row$fastq_1; f2 <- r$row$fastq_2
  switch(fmt,
    bacass_short = data.frame(ID = r$sample, R1 = NAv(f1), R2 = NAv(f2), LongFastQ = "NA", Fast5 = "NA", GenomeSize = "NA"),
    bacass_long  = data.frame(ID = r$sample, R1 = "NA", R2 = "NA", LongFastQ = NAv(f1), Fast5 = "NA", GenomeSize = "NA"),
    taxprofiler  = data.frame(sample = r$sample, run_accession = r$sample,
                              instrument_platform = c(illumina = "ILLUMINA", ont = "OXFORD_NANOPORE", pacbio = "PACBIO_SMRT")[[r$platform]],
                              fastq_1 = f1, fastq_2 = blank(f2), fasta = ""),
    ampliseq     = data.frame(sample = r$sample, fastq_1 = f1, fastq_2 = blank(f2)),
    viralrecon   = data.frame(sample = r$sample, fastq_1 = f1, fastq_2 = blank(f2)),
    rnaseq       = data.frame(sample = r$sample, fastq_1 = f1, fastq_2 = blank(f2), strandedness = "auto"),
    stop("unknown samplesheet format: ", fmt))
}

for (p in unique(rec$pipeline[rec$pipeline != "none"])) {
  sel  <- Filter(function(r) !is.null(r$rule) && paste0("nf-core/", r$rule$pipeline) == p, rows)
  rule <- sel[[1]]$rule
  dir  <- file.path(out, rule$pipeline); dir.create(dir, showWarnings = FALSE)
  # bacass short and long reads use one format with different columns filled.
  ss <- do.call(rbind, lapply(sel, function(r) sheet_row(r$rule$format, r)))
  write.csv(ss, file.path(dir, "samplesheet.csv"), row.names = FALSE, quote = FALSE)

  fixed <- if (has(rule$params)) strsplit(rule$params, ";")[[1]] else character(0)
  if (rule$pipeline == "bacass") {
    types <- unique(vapply(sel, function(r) sub("^bacass_", "", r$rule$rule), ""))
    fixed <- c(sprintf("assembly_type=%s", if (length(types) == 1) types else "FILL_IN"))
  }
  if (rule$pipeline == "ampliseq" && all(vapply(sel, function(r) !has(r$row$fastq_2), TRUE))) fixed <- c(fixed, "single_end=true")
  kv <- do.call(rbind, strsplit(fixed, "=", fixed = TRUE))
  fill <- if (has(rule$fill_in)) strsplit(rule$fill_in, ",")[[1]] else character(0)
  free <- endsWith(fill, "!"); fill <- sub("!$", "", fill)
  yaml <- c(sprintf("# nf-core/%s %s, written by seq-triage (recommend.R, D66).", rule$pipeline, rule$version),
            "# Replace every FILL_IN before running: the pipeline refuses to start with them.",
            "input: samplesheet.csv",
            if (length(fixed)) sprintf("%s: %s", kv[, 1], kv[, 2]),
            if (any(!free)) sprintf("%s: FILL_IN", fill[!free]),
            if (any(free)) c("# Free text, so left out until you supply it (the pipeline stops without it):",
                             sprintf("# %s: <your value>", fill[free])))
  writeLines(yaml, file.path(dir, "params.yaml"))

  readme <- c(sprintf("# nf-core/%s %s", rule$pipeline, rule$version), "",
              sprintf("Recommended by seq-triage for %d sample(s): %s.", length(sel), paste(vapply(sel, `[[`, "", "sample"), collapse = ", ")), "",
              rule$note, "",
              if (length(fill)) c(sprintf("**Fill in first** (in `params.yaml`): %s.", paste(sprintf("`%s`", fill), collapse = ", ")), "") else NULL,
              if (any(grepl("FILL_IN", fixed))) c("**Fill in first:** `assembly_type` (short- and long-read samples are mixed here; run them separately).", "") else NULL,
              "Run from this folder:", "", "```bash",
              sprintf("nextflow run nf-core/%s -r %s -profile <docker|singularity|conda> -params-file params.yaml --outdir results", rule$pipeline, rule$version),
              "```", "",
              "File paths in `samplesheet.csv` are copied from your seq-triage samplesheet; if they are relative, run from the same folder or make them absolute.")
  writeLines(readme, file.path(dir, "README.md"))
}
message(sprintf("Recommendations: %d of %d samples have a next pipeline (%s).", sum(rec$pipeline != "none"), nrow(rec),
                paste(unique(rec$pipeline[rec$pipeline != "none"]), collapse = ", ")))
