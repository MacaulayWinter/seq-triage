#!/usr/bin/env Rscript
#*******************************************************************************
# explore_estimators_D98.R: post hoc look at other ways to estimate off-target
#                           contamination, on layer D rounds 1-2 (CLAUDE.md D98)
#
# DATA ALREADY SEEN: this can only suggest a rule, never validate one.
# Recomputes off_target_fraction from the saved qc_values.tsv three ways and
# scores each with the layer D tolerance (max(20% of true, 0.2 pp); controls
# < 5%), strain mixtures left out (the expected limit):
#   bracken   the current rule: species units (look-alike groups, human out),
#             1% floor, every unit but the top one (Kraken2 for nanopore)
#   kraken2   the same from Kraken2's own species-level counts
#   split     Bracken, but the main genus's reads shared among its species
#             by Kraken2's ratio (aimed at sister species on short reads)
#
# Usage:  Rscript validation/layer_D/explore_estimators_D98.R   (from the repo root)
#*******************************************************************************

o <- "validation/outcomes"
groups <- read.delim("assets/species_groups.tsv", colClasses = "character")
sp_rows <- function(q, tool) {
  b <- q[q$tool == tool, ]; if (!nrow(b)) return(NULL)
  f <- do.call(rbind, strsplit(b$field, "|", fixed = TRUE))
  x <- data.frame(rank = f[, 1], taxid = f[, 2], name = f[, 3], n = as.numeric(b$value))
  x <- x[x$rank == "S" & x$n > 0 & x$taxid != "9606", ]
  x$unit <- ifelse(x$taxid %in% groups$taxid, groups$group[match(x$taxid, groups$taxid)], x$taxid)
  x$genus <- sub(" .*", "", x$name); x$genus[x$genus == "Shigella"] <- "Escherichia"
  x
}
shares <- function(x) { u <- tapply(x$n, x$unit, sum); u / sum(u) }
off <- function(s) { m <- names(sort(s, decreasing = TRUE))[1]; o <- s[names(s) != m]; sum(o[o >= 0.01]) }
estimate <- function(q, how) {
  br <- sp_rows(q, "Bracken"); k2 <- sp_rows(q, "Kraken2"); base <- if (is.null(br)) k2 else br
  if (how == "bracken") return(off(shares(base)))
  if (how == "kraken2") return(off(shares(k2)))
  s <- shares(base); gmain <- base$genus[base$unit == names(sort(s, decreasing = TRUE))[1]][1]
  ug <- tapply(base$genus, base$unit, `[`, 1); in_g <- names(s)[ug[names(s)] == gmain]
  ks <- shares(k2[k2$genus == gmain, ])
  s2 <- s; s2[in_g] <- 0; s2[names(ks)] <- sum(s[in_g]) * ks
  off(s2 / sum(s2))
}
for (r in c("layer_D_round1", "layer_D_round2_std8", "layer_D_round2_full")) {
  qc <- read.delim(file.path(o, r, "qc_values.tsv"), colClasses = "character", quote = "")
  tr <- read.delim(file.path(o, r, "truth.tsv"), colClasses = "character")
  tr <- tr[tr$metric == "off_target_fraction" & tr$score_verdict != "limit", ]
  tf <- as.numeric(tr$true_fraction); tol <- pmax(0.002, 0.2 * tf)
  res <- sapply(c("bracken", "kraken2", "split"), function(h) {
    e <- sapply(tr$sample, function(s) estimate(qc[qc$sample == s, ], h))
    ok <- ifelse(tf == 0, e < 0.05, e > 0 & abs(e - tf) <= tol)
    sprintf("%d/%d", sum(ok), length(ok))
  })
  cat(sprintf("%-22s within tolerance: %s\n", r, paste(names(res), res, sep = " ", collapse = ", ")))
}
