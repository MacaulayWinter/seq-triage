#!/usr/bin/env Rscript
#*******************************************************************************
# build_mixtures.R: contamination mixtures with known answers (D43, layer D)
#
# WHAT IT DOES  For every row of mixtures.tsv, on each platform, takes N reads
#               in total: round(fraction * N) from the contaminant and the
#               rest from the base, then shuffles their order (fixed seed per
#               mixture) so contaminant reads are spread through the file, as
#               they would be in a real sample. Illumina pairs are moved
#               together, so mates stay mates.
# WRITES        <out_dir>/fastq/*.fastq.gz, <out_dir>/samplesheet.csv and
#               <out_dir>/truth.tsv (the true contaminant fraction of each
#               sample, known because we put it there).
# N             Illumina 100,000 pairs (200,000 reads); nanopore 12,000 reads:
#               both comfortably above every read-count limit, so only
#               contamination can change a verdict.
#
# Usage:  Rscript validation/layer_D/build_mixtures.R <sources_dir> <mixtures.tsv> <out_dir>
#*******************************************************************************

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) stop("Usage: build_mixtures.R <sources_dir> <mixtures.tsv> <out_dir>")
script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(script_dir, "..", "lib_fastq.R"))
src_dir <- args[1]
mix     <- read.delim(args[2], colClasses = "character", na.strings = "")
out     <- normalizePath(args[3], mustWork = FALSE)
dir.create(file.path(out, "fastq"), recursive = TRUE, showWarnings = FALSE)

N <- c(ILM = 100000L, ONT = 12000L)

cache <- list()
reads_of <- function(src, platform) {
  key <- paste(src, platform)
  if (is.null(cache[[key]])) {
    files <- if (platform == "ILM") file.path(src_dir, paste0(src, "_ILM_R", 1:2, ".fastq.gz"))
             else file.path(src_dir, paste0(src, "_ONT.fastq.gz"))
    cache[[key]] <<- lapply(files, read_fq)
    message(sprintf("  read %s %s: %s reads", src, platform, format(length(cache[[key]][[1]]$h), big.mark = ",")))
  }
  cache[[key]]
}

sheet <- list(); truth <- list()
for (i in seq_len(nrow(mix))) {
  m <- mix[i, ]
  for (pl in names(N)) {
    n_total <- N[[pl]]
    k <- if (is.na(m$contaminant)) 0L else as.integer(round(as.numeric(m$fraction) * n_total))
    base <- reads_of(m$base, pl)
    if (length(base[[1]]$h) < n_total - k) stop(sprintf("%s %s: base has only %d reads", m$base, pl, length(base[[1]]$h)))
    parts <- lapply(seq_along(base), function(f) subset_fq(base[[f]], seq_len(n_total - k)))
    if (k > 0) {
      cont <- reads_of(m$contaminant, pl)
      if (length(cont[[1]]$h) < k) stop(sprintf("%s %s: contaminant has only %d reads", m$contaminant, pl, length(cont[[1]]$h)))
      parts <- lapply(seq_along(parts), function(f) {
        cf <- subset_fq(cont[[f]], seq_len(k))
        list(h = c(parts[[f]]$h, cf$h), s = c(parts[[f]]$s, cf$s), q = c(parts[[f]]$q, cf$q))
      })
    }
    # One shuffle per mixture, the same for both mates.
    set.seed(20261020 + i)
    ord <- sample(n_total)
    id <- sprintf("D_%s_%s", m$mixture, pl)
    files <- vapply(seq_along(parts), function(f) {
      p <- file.path(out, "fastq", if (length(parts) == 2) sprintf("%s_R%d.fastq.gz", id, f) else sprintf("%s.fastq.gz", id))
      write_fq(subset_fq(parts[[f]], ord), p)
      p
    }, "")
    sheet[[id]] <- data.frame(sample = id, fastq_1 = files[1], fastq_2 = if (length(files) == 2) files[2] else "",
                              platform = if (pl == "ILM") "illumina" else "ont", sample_type = m$sample_type, molecule = "dna")
    truth[[id]] <- data.frame(sample = id, mixture = m$mixture, platform = pl, base = m$base,
                              contaminant = ifelse(is.na(m$contaminant), "", m$contaminant),
                              contaminant_reads = k, total_reads = n_total, true_fraction = k / n_total,
                              metric = m$metric, expected_verdict = m$expected_verdict, score_verdict = m$score_verdict)
    message(sprintf("%-26s %s of %s reads from %s", id, format(k, big.mark = ","), format(n_total, big.mark = ","),
                    ifelse(is.na(m$contaminant), "-", m$contaminant)))
  }
}
write.csv(do.call(rbind, sheet), file.path(out, "samplesheet.csv"), row.names = FALSE, quote = FALSE)
write.table(do.call(rbind, truth), file.path(out, "truth.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
message(sprintf("\n%d samples written to %s", length(sheet), out))
