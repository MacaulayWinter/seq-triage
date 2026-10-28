#!/usr/bin/env Rscript
# =============================================================================
# score_fresh.R  (D96, pass mark P7)
# Scores the pathogen flags of the fresh-data run against ENA's species label.
# Independent of bin/triage_report.R: it counts species from qc_values.tsv with
# its own code, and looks families up with NCBI E-utilities.
#
#   label not confirmed  the run's main unit does not contain the ENA species
#                        (replaced by the next backup, as pre-registered)
#   P7a (listed arm)     the ENA species is flagged, right status, no caveat
#   P7b (unlisted arm)   the main unit is not flagged
#   P7c (both)           every other flag carries the caveat exactly when its
#                        family equals the main species' family
#   at least 8 confirmed runs per arm, else INCONCLUSIVE
#
# Usage:  score_fresh.R <panel.tsv> <qc_values.tsv> <pathogen_flags.tsv> <pathogen_list_dir> <out.tsv>
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) stop("Usage: score_fresh.R <panel.tsv> <qc_values.tsv> <pathogen_flags.tsv> <pathogen_list_dir> <out.tsv>")
rd <- function(f) read.delim(f, colClasses = "character", na.strings = "", quote = "")
panel <- rd(args[1]); qc <- rd(args[2]); flags <- rd(args[3]); pl <- rd(file.path(args[4], "pathogen_list.tsv"))
flags$caveat[is.na(flags$caveat)] <- ""            # read.delim turns blank cells into NA
repo   <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "..", ".."))
groups <- rd(file.path(repo, "assets/species_groups.tsv"))

fam_cache <- new.env()
efam <- function(ids) {
  unique(unlist(lapply(unique(ids), function(id) {
    if (is.null(fam_cache[[id]])) {
      Sys.sleep(0.4)                                 # NCBI: at most 3 requests a second without a key
      x <- paste(readLines(sprintf("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id=%s", id), warn = FALSE), collapse = "")
      f <- regmatches(x, regexpr("<TaxId>[0-9]+</TaxId>\\s*<ScientificName>[^<]*</ScientificName>\\s*<Rank>family</Rank>", x))
      fam_cache[[id]] <- if (length(f)) sub(".*<TaxId>([0-9]+)</TaxId>.*", "\\1", f) else NA_character_
    }
    fam_cache[[id]]
  })))
}

ran <- panel[panel$run_accession %in% qc$sample, ]
res <- do.call(rbind, lapply(seq_len(nrow(ran)), function(i) {
  r <- ran[i, ]; s <- r$run_accession
  q <- qc[qc$sample == s, ]
  b <- q[q$tool == "Bracken", ]
  sp <- if (nrow(b)) b else q[q$tool == "Kraken2", ]
  f <- do.call(rbind, strsplit(sp$field, "|", fixed = TRUE))
  sp <- data.frame(rank = f[, 1], taxid = f[, 2], n = as.numeric(sp$value))
  sp <- sp[sp$rank == "S" & sp$taxid != "9606" & sp$n > 0, ]
  if (sum(sp$n) < 10) return(data.frame(run = s, arm = r$arm, ena = r$scientific_name, main = "", confirmed = FALSE,
                                        p7a = NA, p7b = NA, p7c = NA, detail = "no species-level estimate"))
  sp$unit <- ifelse(sp$taxid %in% groups$taxid, groups$group[match(sp$taxid, groups$taxid)], paste0("taxid:", sp$taxid))
  tot <- tapply(sp$n, sp$unit, sum); share <- sort(tot / sum(tot), decreasing = TRUE)
  main_u <- names(share)[1]
  main_ids <- sp$taxid[sp$unit == main_u]
  # The ENA taxid may be newer than the database's: compare through the list's
  # merged-taxid rows too.
  truth_ids <- unique(c(r$tax_id, pl$match_taxid[pl$taxid == r$tax_id]))
  confirmed <- any(main_ids %in% truth_ids)
  fl <- flags[flags$sample == s & flags$result == "listed", ]
  unit_of <- function(k) unique(sp$unit[sp$taxid %in% pl$match_taxid[pl$taxid %in% strsplit(fl$ncbi_taxid[k], ",")[[1]]]])[1]
  fu <- vapply(seq_len(nrow(fl)), unit_of, "")
  main_flag <- which(fu == main_u)
  p7a <- if (r$arm == "listed") {
    st <- pl$status[pl$taxid == r$tax_id | pl$match_taxid == r$tax_id][1]
    length(main_flag) == 1 && fl$list_status[main_flag] == st && !nzchar(fl$caveat[main_flag])
  } else NA
  p7b <- if (r$arm == "unlisted") length(main_flag) == 0 else NA
  wrong <- character(0)
  for (k in setdiff(seq_len(nrow(fl)), main_flag)) {
    same <- length(intersect(efam(sp$taxid[sp$unit == fu[k]]), efam(main_ids))) > 0
    if (same != nzchar(fl$caveat[k])) wrong <- c(wrong, fl$species[k])
  }
  data.frame(run = s, arm = r$arm, ena = r$scientific_name, main = sprintf("%s %.1f%%", main_u, 100 * share[[1]]),
             confirmed = confirmed, p7a = p7a, p7b = p7b, p7c = !length(wrong),
             detail = paste(c(sprintf("%s (%s%s)", fl$species, fl$share, ifelse(nzchar(fl$caveat), ", caveat", "")),
                              if (length(wrong)) paste("caveat wrong:", paste(wrong, collapse = ", "))), collapse = "; "))
}))
write.table(res, args[5], sep = "\t", quote = FALSE, row.names = FALSE)
options(width = 250)
print(res, row.names = FALSE, right = FALSE)
c_ <- res[res$confirmed, ]
for (a in c("listed", "unlisted")) {
  x <- c_[c_$arm == a, ]
  mark <- if (a == "listed") x$p7a else x$p7b
  cat(sprintf("%s arm: %d confirmed (%d not confirmed); %s %d/%d; P7c %d/%d -> %s\n", a, nrow(x), sum(res$arm == a & !res$confirmed),
              if (a == "listed") "P7a" else "P7b", sum(mark), nrow(x), sum(x$p7c), nrow(x),
              if (nrow(x) < 8) "INCONCLUSIVE (fewer than 8 confirmed)" else if (all(mark) && all(x$p7c)) "PASS" else "FAIL"))
}
