#!/usr/bin/env Rscript
#*******************************************************************************
# recommend.R: which nf-core pipeline each ready sample should go to next,
#              with a ready-to-run samplesheet and params file (CLAUDE.md D66)
#
# WHAT IT DOES  For every sample in the triage summary, either names ONE
#               downstream pipeline (from assets/downstream_pipelines.tsv,
#               pinned versions) or says why there is none. Then writes, per
#               rule (one folder each, e.g. bacass_short and bacass_long, so
#               every folder runs as it is), a samplesheet in that pipeline's
#               own format, a params.yaml and a README with the command.
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
# WRITES        <outdir>/recommendations.tsv and <outdir>/<rule>/...
#
# Usage:  recommend.R <triage_summary.tsv> <triage_metrics.tsv> <samplesheet.csv> <downstream_pipelines.tsv> <outdir>
#*******************************************************************************

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) stop("Usage: recommend.R <triage_summary.tsv> <triage_metrics.tsv> <samplesheet.csv> <downstream_pipelines.tsv> <outdir>")
rt   <- function(f) read.delim(f, colClasses = "character", na.strings = "", quote = "", check.names = FALSE)
summ <- rt(args[1]); mets <- rt(args[2]); rules <- rt(args[4]); out <- args[5]
sheet <- read.csv(args[3], colClasses = "character", na.strings = "", check.names = FALSE, strip.white = TRUE)
for (col in c("fastq_1", "fastq_2", "bam", "fastq_dir", "group")) if (!col %in% names(sheet)) sheet[[col]] <- NA_character_
sheet$sample <- trimws(sheet$sample)
dir.create(out, recursive = TRUE, showWarnings = FALSE)

# SARS-CoV-2 under every name it has in NCBI taxonomy versions: the virus,
# its old species ("SARS-related coronavirus") and the 2023 ICTV binomial,
# "Betacoronavirus pandemicum", which the Oct 2025 Kraken2 databases use. The
# first version matched only the virus name, so validation layer E's real
# SARS-CoV-2 amplicons would have been sent to ampliseq (found 2026-10-22).
SARSCOV2 <- c("Severe acute respiratory syndrome coronavirus 2", "Severe acute respiratory syndrome-related coronavirus",
              "Betacoronavirus pandemicum")
has <- function(x) !is.na(x) & nzchar(x)
`%|%` <- function(a, b) ifelse(is.na(a), b, a)
in_list <- function(value, list) value %in% strsplit(list, ",")[[1]]

# WHAT IT DOES  Decides for one sample. RETURNS list(rule = row of rules or
#               NULL, reason = text).
# row: the sample's samplesheet rows (several for a sample in several
# files, D75), in samplesheet order.
decide <- function(s, row) {
  none <- function(why) list(rule = NULL, reason = why)
  if (!s$verdict %in% c("PASS", "WARN"))
    return(none(sprintf("Verdict is %s: resolve the triage findings first.", s$verdict)))
  if (any(grepl("\\.(fa|fasta|fna)(\\.(gz|bz2|xz|zst))?$", c(row$fastq_1, row$fastq_2))))
    return(none("Reads given as FASTA (no quality values): the downstream pipelines take FASTQ; convert or obtain the FASTQ first."))
  if (any(has(row$fastq_dir)))
    return(none("Reads given as a folder (fastq_dir): list the files in the next pipeline's samplesheet, one row each, or merge them first."))
  if (any(has(row$bam)))
    return(none("BAM input: the downstream pipelines take FASTQ. Convert first (samtools collate | samtools fastq)."))
  flags <- if (has(s$detection_flags)) strsplit(s$detection_flags, ",")[[1]] else character(0)
  conflict <- intersect(flags, c("TYPE_CONTRADICTS", "BAM_HEADER_DISAGREES"))
  if (length(conflict))
    return(none(sprintf("Metadata conflict (%s): review the sample first.", paste(conflict, collapse = ", "))))
  platform <- switch(s$routed_by, detected = s$detected_platform, declared = s$declared_platform, NA)
  if (!has(platform)) return(none("Platform not confirmed by the reads."))
  type <- if (has(s$declared_type) && s$declared_type %in% c("isolate", "metagenome", "amplicon")) s$declared_type
          else if (has(s$detected_type) && s$type_confidence %in% "high") s$detected_type else NA
  # D78: the library's structure, read from the reads, decides the kind of
  # analysis before the isolate / metagenome / amplicon question does.
  if ("BARCODE_READ" %in% flags) type <- "single_cell"
  if ("BISULFITE" %in% flags)    type <- "methylation"
  if (!has(type)) return(none("Sample type not declared and not detected with high confidence: declare sample_type."))
  # Single-cell libraries are almost always RNA (D78); otherwise DNA is assumed.
  molecule <- if (has(s$declared_molecule)) s$declared_molecule else if (type == "single_cell") "rna" else "dna"
  dom <- mets$value[mets$sample == s$sample & mets$metric == "dominant_species"]
  vir <- suppressWarnings(as.numeric(sub("%", "", mets$value[mets$sample == s$sample & mets$metric == "viral_fraction"]))) / 100
  bac <- suppressWarnings(as.numeric(sub("%", "", mets$value[mets$sample == s$sample & mets$metric == "bacterial_fraction"]))) / 100
  # D76: a declared human target goes to sarek, never to a bacterial
  # assembler; bacass also needs most classified reads to be bacterial when
  # the screen ran (a yeast isolate is not a job for bacass).
  human <- has(s$declared_target) && s$declared_target %in% c("9606", "9605")
  bacterial <- !human && (!length(bac) || is.na(bac) || bac >= 0.5)
  # target: SARS-CoV-2 (by name, every taxonomy version); another virus (by
  # taxonomy: most classified reads under Viruses, D67); or a cellular
  # organism/community. A mumps amplicon must not go to ampliseq, which is
  # for marker genes (layer E, 2026-10-22).
  target <- if (!length(dom) || !has(dom)) NA
            else if (any(startsWith(dom, SARSCOV2))) "sarscov2"
            else if (length(vir) && !is.na(vir) && vir >= 0.5) "virus"
            else if (length(vir) && !is.na(vir)) "other" else NA
  fits <- which(vapply(seq_len(nrow(rules)), function(i) {
    r <- rules[i, ]
    in_list(platform, r$platform) && r$sample_type == type && in_list(molecule, r$molecule) &&
      (r$target == "any" || (r$target == "sarscov2" && target %in% "sarscov2") ||
       (r$target == "not_sarscov2" && target %in% "other") ||
       (r$target == "human" && human) || (r$target == "bacterial" && bacterial))
  }, TRUE))
  # ampliseq's samplesheet takes one row per sample (its schema says so,
  # whatever its code allows), so several files must be merged first (D75).
  fits <- fits[!(rules$pipeline[fits] == "ampliseq" & nrow(row) > 1)]
  if (length(fits) == 0 && nrow(row) > 1 && type == "amplicon" && !target %in% c("sarscov2", "virus"))
    return(none("Several files for one sample: nf-core/ampliseq takes one row per sample, so merge the files first (cat), then use nf-core/ampliseq."))
  # D78: the files as given cannot go to the next pipeline as they are.
  if (length(fits) && s$detected_layout %in% "interleaved")
    return(none(sprintf("Interleaved FASTQ (both mates in one file): split it into R1 and R2 first (e.g. seqtk seq -1 and seqtk seq -2, or reformat.sh in= out1= out2=), then use nf-core/%s %s.",
                        rules$pipeline[fits[1]], rules$version[fits[1]])))
  if ("INDEX_READ_AS_MATE" %in% flags)
    return(none("fastq_2 looks like an index read, not R2: fix the samplesheet first."))
  if (length(fits)) {
    # D80: depth of coverage, where the next step depends on it.
    dep <- suppressWarnings(as.numeric(sub("x$", "", mets$value[mets$sample == s$sample & mets$metric == "depth_of_coverage"])))
    dep_note <- if (!length(dep) || is.na(dep) || !rules$pipeline[fits[1]] %in% c("bacass", "sarek")) NULL
                else if (rules$pipeline[fits[1]] == "bacass" && dep < 30) sprintf("depth about %.0fx: below 30x, so the assembly will be fragmented; more data would help", dep)
                else if (rules$pipeline[fits[1]] == "sarek" && dep < 20) sprintf("depth about %.0fx: below 20x, thin for variant calling", dep)
                else sprintf("depth about %.0fx", dep)
    note <- c(dep_note, if (!has(s$declared_molecule)) sprintf("molecule not declared, so %s assumed", toupper(molecule)) else NULL,
              if (!has(s$declared_type)) sprintf("%s detected from the reads", c(single_cell = "single-cell library", methylation = "methylation library")[type] %|% sprintf("sample type %s", type)) else NULL)
    return(list(rule = rules[fits[1], ], reason = paste(c(rules$note[fits[1]], if (length(note)) sprintf("(%s)", paste(note, collapse = "; "))), collapse = " ")))
  }
  # Known gaps, said plainly.
  if (type == "isolate" && human)
    return(none(sprintf("Human %s on %s: no pipeline set up here for that combination.", molecule, platform)))
  if (type == "isolate" && molecule == "dna" && !bacterial)
    return(none(sprintf("Isolate whose classified reads are mostly not bacterial (bacterial: %.0f%%): no pipeline set up here for non-bacterial organisms (a eukaryote may also be missing from the Kraken2 database).", 100 * bac)))
  if (type == "amplicon" && platform != "illumina")
    return(none("Long-read amplicons: no pipeline set up here (nf-core/viralrecon expects the sequencer's barcode folders)."))
  if (type == "amplicon" && is.na(target))
    return(none("Amplicon target unknown without the contamination screen (--kraken2_db): marker gene or virus?"))
  if (type == "amplicon" && target == "virus")
    return(none(sprintf("Viral amplicon other than SARS-CoV-2 (%s): no pipeline set up here.", sub(" \\(.*$", "", dom))))
  if (type == "isolate" && platform == "pacbio")
    return(none("PacBio HiFi isolates: no pipeline set up here yet."))
  if (type == "metagenome" && platform == "element")
    return(none("Element metagenome: nf-core/taxprofiler does not list Element as a platform; nf-core/mag (assembly) does."))
  if (molecule == "rna" && type != "isolate")
    return(none(sprintf("RNA %s (e.g. metatranscriptome): no pipeline set up here yet.", type)))
  none(sprintf("No pipeline in assets/downstream_pipelines.tsv for %s / %s / %s.", platform, type, molecule))
}

rows <- lapply(seq_len(nrow(summ)), function(i) {
  s   <- summ[i, ]
  row <- sheet[sheet$sample == s$sample, , drop = FALSE]
  d   <- decide(s, row)
  list(sample = s$sample, row = row, platform = switch(s$routed_by, detected = s$detected_platform, declared = s$declared_platform, NA),
       rule = d$rule, reason = d$reason)
})

# D83: samples of one isolate sequenced on a short- and a long-read platform
# (same samplesheet `group`) are assembled together: one bacass hybrid row
# (R1, R2 and the long reads), which resolves plasmids and repeats. Only
# when both are ready for bacass on their own and look like the same
# organism (same most common species, when the screen ran).
group_of <- function(id) { g <- sheet$group[sheet$sample == id]; if (length(g) && has(g[1])) g[1] else NA_character_ }
hybrid_rule <- rules[rules$rule == "bacass_short", ][1, ]
hybrid_rule$rule <- "bacass_hybrid"; hybrid_rule$format <- "bacass_hybrid"; hybrid_rule$params <- "assembly_type=hybrid"
hybrid_rule$note <- "Hybrid assembly of short and long reads of one isolate: long reads span repeats and close plasmids, short reads correct errors."
hybrids <- list()
for (g in unique(na.omit(vapply(rows, function(r) group_of(r$sample), "")))) {
  idx   <- which(vapply(rows, function(r) group_of(r$sample) %in% g, TRUE))
  short <- idx[vapply(rows[idx], function(r) !is.null(r$rule) && r$rule$rule == "bacass_short", TRUE)]
  long  <- idx[vapply(rows[idx], function(r) !is.null(r$rule) && r$rule$rule == "bacass_long", TRUE)]
  if (length(short) != 1 || length(long) != 1) next
  # One samplesheet row each: bacass merges repeated rows of one ID, which
  # would repeat the long reads once per short-read lane.
  if (nrow(rows[[short]]$row) != 1 || nrow(rows[[long]]$row) != 1) next
  sp <- vapply(c(short, long), function(i) { v <- mets$value[mets$sample == rows[[i]]$sample & mets$metric == "dominant_species"]
                                              if (length(v) && has(v[1])) sub(" \\(.*$", "", v[1]) else NA_character_ }, "")
  if (all(!is.na(sp)) && sp[1] != sp[2]) {
    for (i in c(short, long)) rows[[i]]$reason <- paste(rows[[i]]$reason, sprintf("(group %s: the short- and long-read samples look like different organisms, %s and %s, so no hybrid assembly is advised)", g, sp[1], sp[2]))
    next
  }
  for (i in c(short, long)) {
    rows[[i]]$rule   <- hybrid_rule
    rows[[i]]$reason <- sprintf("%s (group %s: %s short reads + %s long reads)", hybrid_rule$note, g, rows[[short]]$sample, rows[[long]]$sample)
  }
  hybrids[[g]] <- list(short = rows[[short]], long = rows[[long]])
}

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
# One output row per samplesheet row: the pipelines that accept a sample
# in several files take repeated rows (bacass, viralrecon, rnaseq; taxprofiler
# with a different run_accession per row), checked in their code (D75).
sheet_row <- function(fmt, r) {
  f1 <- r$row$fastq_1; f2 <- r$row$fastq_2
  run <- if (nrow(r$row) > 1) sprintf("%s_%d", r$sample, seq_len(nrow(r$row))) else r$sample
  switch(fmt,
    bacass_short = data.frame(ID = r$sample, R1 = NAv(f1), R2 = NAv(f2), LongFastQ = "NA", Fast5 = "NA", GenomeSize = "NA"),
    bacass_long  = data.frame(ID = r$sample, R1 = "NA", R2 = "NA", LongFastQ = NAv(f1), Fast5 = "NA", GenomeSize = "NA"),
    taxprofiler  = data.frame(sample = r$sample, run_accession = run,
                              instrument_platform = c(illumina = "ILLUMINA", dnbseq = "DNBSEQ", ont = "OXFORD_NANOPORE",
                                                      pacbio = "PACBIO_SMRT")[[r$platform]],
                              fastq_1 = f1, fastq_2 = blank(f2), fasta = ""),
    ampliseq     = data.frame(sample = r$sample, fastq_1 = f1, fastq_2 = blank(f2)),
    viralrecon   = data.frame(sample = r$sample, fastq_1 = f1, fastq_2 = blank(f2)),
    rnaseq       = data.frame(sample = r$sample, fastq_1 = f1, fastq_2 = blank(f2), strandedness = "auto"),
    # sarek (D76): one row per lane, lanes distinct (its uniqueEntries:
    # patient + sample + lane); patient = the sample ID.
    scrnaseq     = data.frame(sample = r$sample, fastq_1 = f1, fastq_2 = blank(f2)),
    methylseq    = data.frame(sample = r$sample, fastq_1 = f1, fastq_2 = blank(f2)),
    sarek        = data.frame(patient = r$sample, sample = r$sample, lane = sprintf("L%d", seq_len(nrow(r$row))),
                              fastq_1 = f1, fastq_2 = blank(f2)),
    stop("unknown samplesheet format: ", fmt))
}

for (rl in unique(rec$rule[rec$pipeline != "none"])) {
  sel  <- Filter(function(r) !is.null(r$rule) && r$rule$rule == rl, rows)
  rule <- sel[[1]]$rule
  dir  <- file.path(out, rl); dir.create(dir, showWarnings = FALSE)
  ss <- if (rl == "bacass_hybrid")
          do.call(rbind, lapply(names(hybrids), function(g) data.frame(ID = g, R1 = NAv(hybrids[[g]]$short$row$fastq_1[1]),
                  R2 = NAv(hybrids[[g]]$short$row$fastq_2[1]), LongFastQ = NAv(hybrids[[g]]$long$row$fastq_1[1]), Fast5 = "NA", GenomeSize = "NA")))
        else do.call(rbind, lapply(sel, function(r) sheet_row(r$rule$format, r)))
  write.csv(ss, file.path(dir, "samplesheet.csv"), row.names = FALSE, quote = FALSE)

  fixed <- if (has(rule$params)) strsplit(rule$params, ";")[[1]] else character(0)
  if (rule$pipeline == "ampliseq" && all(vapply(sel, function(r) all(!has(r$row$fastq_2)), TRUE))) fixed <- c(fixed, "single_end=true")
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
              "Run from this folder:", "", "```bash",
              sprintf("nextflow run nf-core/%s -r %s -profile <docker|singularity|conda> -params-file params.yaml --outdir results", rule$pipeline, rule$version),
              "```", "",
              "File paths in `samplesheet.csv` are copied from your seq-triage samplesheet; if they are relative, run from the same folder or make them absolute.")
  writeLines(readme, file.path(dir, "README.md"))
}
message(sprintf("Recommendations: %d of %d samples have a next pipeline (%s).", sum(rec$pipeline != "none"), nrow(rec),
                paste(unique(rec$pipeline[rec$pipeline != "none"]), collapse = ", ")))


# Analysis menu (D82) ------------------------------------------------------------
# WHAT IT DOES  For every sample, the analyses its data can support next, each
#               with a status and the reason, so the user sees the options and
#               what is missing, not only one pipeline:
#                 ready     the data supports it now
#                 possible  it can be run, with the caveat given
#                 not now   something is missing (more data, more samples,
#                           a declaration); the reason says what
#               Decided from what the triage found (type, molecule, platform,
#               library structure, depth, rRNA share, target) and from the
#               whole batch (comparative analyses need several samples).
#               Limits are provisional and stated in each reason.
# WRITES        <outdir>/analyses.tsv and <outdir>/analyses.md
metric <- function(sample, name) {
  v <- mets$value[mets$sample == sample & mets$metric == name]
  if (length(v)) v[1] else NA_character_
}
num_of <- function(x) suppressWarnings(as.numeric(gsub("[,x%]", "", x)))
dom_species <- function(sample) { v <- metric(sample, "dominant_species"); if (has(v)) sub(" \\(.*$", "", v) else NA_character_ }
SHORT <- c("illumina", "dnbseq", "element"); LONG <- c("ont", "pacbio")
kind_of <- function(s) {
  fl <- if (has(s$detection_flags)) strsplit(s$detection_flags, ",")[[1]] else character(0)
  if ("BARCODE_READ" %in% fl) return("single_cell")
  if ("BISULFITE" %in% fl)    return("methylation")
  t <- if (has(s$declared_type) && s$declared_type %in% c("isolate", "metagenome", "amplicon")) s$declared_type
       # A medium "amplicon" means the reads cannot settle it (long reads,
       # RNA with uneven starts, D85), so it does not open the amplicon menu.
       else if (has(s$detected_type) && (s$type_confidence %in% "high" ||
                (s$type_confidence %in% "medium" && s$detected_type != "amplicon"))) s$detected_type else NA
  mol <- if (has(s$declared_molecule)) s$declared_molecule else "dna"
  if (is.na(t)) return(NA_character_)
  if (t == "isolate" && "declared_target" %in% names(s) && s$declared_target %in% c("9606", "9605")) return(if (mol == "rna") "rna_isolate" else "human")
  if (mol == "rna") return(if (t == "metagenome") "metatranscriptome" else if (t == "isolate") "rna_isolate" else "amplicon")
  t
}
ready_samples <- summ[summ$verdict %in% c("PASS", "WARN"), ]
kinds <- setNames(vapply(seq_len(nrow(ready_samples)), function(i) kind_of(ready_samples[i, ]) %|% "", ""), ready_samples$sample)
iso_species <- vapply(names(kinds)[kinds == "isolate"], dom_species, "")
row_ <- function(sample, analysis, status, pipeline, reason) data.frame(sample = sample, analysis = analysis, status = status, pipeline = pipeline, reason = reason)
menu_for <- function(s) {
  id <- s$sample
  if (!s$verdict %in% c("PASS", "WARN"))
    return(row_(id, "any", "not now", "", sprintf("Verdict is %s: resolve the triage findings first.", s$verdict)))
  k <- kind_of(s)
  platform <- switch(s$routed_by, detected = s$detected_platform, declared = s$declared_platform, s$detected_platform)
  long <- platform %in% LONG
  depth <- num_of(metric(id, "depth_of_coverage"))
  reads <- num_of(metric(id, if (long) "number_of_reads" else "total_reads"))
  rrna  <- num_of(metric(id, "rrna_fraction"))
  bac   <- num_of(metric(id, "bacterial_fraction"))
  n_of  <- function(kind) sum(kinds == kind)
  if (is.na(k))
    return(row_(id, "any", "not now", "", "Sample type not known (not declared, and not detected with confidence): declare sample_type to see the analyses it supports."))
  rrna_note <- if (!is.na(rrna) && rrna > 25) sprintf(" rRNA is %.0f%% of reads, so fewer reads describe genes; sequence deeper or improve depletion.", rrna) else ""
  switch(k,
    isolate = {
      sp <- dom_species(id); same <- if (has(sp)) sum(iso_species == sp, na.rm = TRUE) else 0
      non_bact <- !is.na(bac) && bac < 50
      asm_status <- if (non_bact) "not now" else if (is.na(depth)) "possible" else if (depth >= 30) "ready" else if (depth >= 10) "possible" else "not now"
      asm_why <- if (non_bact) "Most classified reads are not bacterial; no assembler for non-bacterial genomes is set up here."
                 else if (is.na(depth)) "Depth unknown: declare genome_size (or target_taxid) to check it is enough (30x or more)."
                 else if (depth >= 30) sprintf("Depth about %.0fx: enough for de novo assembly.", depth)
                 else if (depth >= 10) sprintf("Depth about %.0fx: an assembly will be fragmented (30x or more is better).", depth)
                 else sprintf("Depth about %.0fx: too thin for de novo assembly; sequence more (30x or more).", depth)
      rbind(row_(id, "Genome assembly and annotation", asm_status, sprintf("nf-core/bacass 2.6.1 (%s)", if (long) "long" else "short"), asm_why),
            row_(id, "Species confirmation and typing (ANI, MLST)", if (asm_status == "not now" && !non_bact && !is.na(depth) && depth >= 5) "possible" else asm_status, "",
                 "From the assembly (e.g. fastANI, mlst), or read-based typing at lower depth."),
            row_(id, "AMR and virulence genes", if (asm_status == "ready") "ready" else if (asm_status == "possible") "possible" else "not now", "nf-core/funcscan 4.0.0 (on the assembly)",
                 if (asm_status == "ready") "On the assembled contigs (AMRFinderPlus, RGI, ABRicate)." else "Needs an assembly; with thin data genes split across contigs can be missed."),
            {
              g <- group_of(id); hyb <- !is.na(g) && g %in% names(hybrids)
              row_(id, "Plasmids and mobile genetic elements", if (hyb || (long && asm_status == "ready")) "ready" else if (asm_status == "not now") "not now" else "possible",
                   if (hyb) "nf-core/bacass 2.6.1 (hybrid)" else "",
                   if (hyb) sprintf("Short and long reads of this isolate (group %s): a hybrid assembly closes plasmids and places mobile elements.", g)
                   else if (long) "Long reads can close plasmids and resolve where mobile elements sit."
                   else "Short reads fragment plasmids and repeats; long reads of the same isolate (give both the same group) allow a hybrid assembly that resolves them.")
            },
            row_(id, "Comparative genomics (pangenome, core-genome phylogeny)", if (same >= 3) "ready" else "not now", "",
                 if (same >= 3) sprintf("%d isolates of %s in this batch.", same, sp)
                 else if (!has(sp)) "Needs several isolates of one species; the species is not known here (turn on the contamination screen to name it)."
                 else sprintf("Needs several isolates of one species; this batch has %d of %s.", same, sp)))
    },
    human = rbind(
      row_(id, "Germline variant calling", if (is.na(depth)) "possible" else if (depth >= 20) "ready" else if (depth >= 10) "possible" else "not now", "nf-core/sarek 3.10.0",
           if (is.na(depth)) "Depth unknown." else sprintf("Depth about %.1fx (20x or more is the usual minimum for germline calls).", depth)),
      row_(id, "Somatic variant calling", "not now", "nf-core/sarek 3.10.0", "Needs a matched normal sample, which the reads cannot show; declare tumour/normal in sarek's samplesheet.")),
    metagenome = rbind(
      row_(id, "Taxonomic profiling", if (platform == "element") "not now" else "ready", "nf-core/taxprofiler 2.0.1",
           if (platform == "element") "taxprofiler does not list Element as a platform." else "Read-based profiles (Kraken2, MetaPhlAn and others)."),
      row_(id, "Assembly and binning (MAGs)", if (is.na(reads)) "possible" else if (reads >= 5e6) "ready" else if (reads >= 1e6) "possible" else "not now", "nf-core/mag 5.5.0",
           if (is.na(reads)) "Read count unknown." else sprintf("%s reads: binning needs millions (5 million or more is a usual minimum; deeper for diverse communities).", format(reads, big.mark = ","))),
      row_(id, "Resistome (antimicrobial resistance genes)", "possible", "nf-core/funcscan 4.0.0 (on MAG contigs)", "Contig-based after assembly, or read-based mapping to a resistance database; read-based is more sensitive at low depth."),
      row_(id, "Diversity across samples (alpha, beta, differential abundance)", if (n_of("metagenome") >= 3) "ready" else "not now", "",
           sprintf("Needs several metagenomes; this batch has %d.", n_of("metagenome")))),
    amplicon = {
      dom <- metric(id, "dominant_species"); vir <- num_of(metric(id, "viral_fraction"))
      sc2 <- has(dom) && any(startsWith(dom, SARSCOV2))
      if (sc2) row_(id, "Viral consensus, variants and lineage", if (platform %in% SHORT) "ready" else "possible", "nf-core/viralrecon 3.0.0", "SARS-CoV-2 amplicons: needs the primer scheme.")
      else if (!is.na(vir) && vir >= 50) row_(id, "Viral consensus and variants", "not now", "", "A viral amplicon other than SARS-CoV-2: no pipeline set up here.")
      else rbind(row_(id, "Marker-gene analysis (ASVs, taxonomy)", if (platform %in% SHORT) "ready" else "not now", "nf-core/ampliseq 2.18.0",
                      if (platform %in% SHORT) "Needs the primer sequences." else "Long-read amplicons: no pipeline set up here."),
                 row_(id, "Diversity across samples", if (n_of("amplicon") >= 3) "ready" else "not now", "", sprintf("Needs several amplicon samples; this batch has %d.", n_of("amplicon"))))
    },
    rna_isolate = rbind(
      row_(id, "Expression quantification", if (platform %in% SHORT) "ready" else "possible", "nf-core/rnaseq 3.27.0", paste0("Needs the reference genome and annotation.", rrna_note)),
      row_(id, "Differential expression", "not now", "nf-core/rnaseq 3.27.0, then DESeq2 / edgeR", "Needs biological replicates (2-3 or more per condition), which the reads cannot show.")),
    metatranscriptome = rbind(
      row_(id, "De novo metatranscriptome (genes, expression)", "possible", "nf-core/metatdenovo 1.4.1", paste0("Assembles and quantifies community transcripts.", rrna_note)),
      {
        g <- group_of(id)
        mates <- if (is.na(g)) character(0) else names(kinds)[kinds == "metagenome" & vapply(names(kinds), group_of, "") %in% g]
        row_(id, "Activity vs abundance (RNA/DNA ratio)", if (length(mates)) "ready" else if (n_of("metagenome") >= 1) "possible" else "not now", "",
             if (length(mates)) sprintf("Matched DNA metagenome in group %s: %s.", g, paste(mates, collapse = ", "))
             else if (n_of("metagenome") >= 1) "This batch also has DNA metagenomes: give each RNA sample and the DNA of the same sample one group to pair them."
             else "Needs a matched DNA metagenome of the same sample.")
      }),
    single_cell = row_(id, "Single-cell expression (cell barcodes, UMI counts)", "ready", "nf-core/scrnaseq 4.2.0", "Needs the reference genome and annotation."),
    methylation = row_(id, "DNA methylation calling", "ready", "nf-core/methylseq 4.2.0", "Needs the reference genome; say whether the library is WGBS, EM-seq or RRBS."),
    row_(id, "any", "not now", "", "No analyses set up for this kind of sample yet."))
}
menu <- do.call(rbind, lapply(seq_len(nrow(summ)), function(i) menu_for(summ[i, ])))
write.table(menu, file.path(out, "analyses.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
md <- c("# What can be done with each sample", "",
        "From seq-triage's findings (D82). **ready**: the data supports it now; **possible**: it can be run, with the caveat given; **not now**: something is missing, and the reason says what. Limits are provisional.", "")
for (id in unique(menu$sample)) {
  m <- menu[menu$sample == id, ]
  md <- c(md, sprintf("## %s", id), "", "| Analysis | Status | Pipeline | Why |", "|---|---|---|---|",
          sprintf("| %s | %s | %s | %s |", m$analysis, m$status, ifelse(nzchar(m$pipeline), m$pipeline, "-"), m$reason), "")
}
writeLines(md, file.path(out, "analyses.md"))
message(sprintf("Analysis menu: %d analyses ready, %d possible, %d not now, over %d samples.",
                sum(menu$status == "ready"), sum(menu$status == "possible"), sum(menu$status == "not now"), length(unique(menu$sample))))
