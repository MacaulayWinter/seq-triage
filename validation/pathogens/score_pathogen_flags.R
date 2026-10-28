#!/usr/bin/env Rscript
# =============================================================================
# score_pathogen_flags.R  (D96, pass marks P2 and P3)
# Scores the pathogen flags of a re-run layer D report against the known
# composition of each mixture. Independent of bin/triage_report.R: it counts
# species from qc_values.tsv with its own code, and takes the truth from the
# mixture sources (organism -> NCBI taxid written here, not looked up).
#
# Expected flags: every counted unit (>= 1% of species-level reads, look-alike
# groups as one) that contains a species of the mixture whose taxid is on the
# list. PASS (amended P2/P3) when the report flags every expected unit with
# the right status, and every other flag carries the close-relative caveat.
# The caveat's correctness (same family as the main species) is checked with
# families fetched here from NCBI E-utilities, not the pipeline's table.
#
# Usage:  score_pathogen_flags.R <outcomes_dir> <sources.tsv> <pathogen_list.tsv> <triage_pathogens.tsv> <out.tsv>
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) stop("Usage: score_pathogen_flags.R <outcomes_dir> <sources.tsv> <pathogen_list.tsv> <triage_pathogens.tsv> <out.tsv>")
rd <- function(f) read.delim(f, colClasses = "character", na.strings = "", check.names = FALSE, quote = "")
run <- args[1]
truth   <- rd(file.path(run, "truth.tsv"))
qc      <- rd(file.path(run, "qc_values.tsv"))
sources <- rd(args[2])
plist   <- rd(args[3])
flags   <- rd(args[4])
repo    <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "..", ".."))
groups  <- rd(file.path(repo, "assets/species_groups.tsv"))

ORGANISM_TAXID <- c("Escherichia coli" = "562", "Klebsiella pneumoniae" = "573", "Klebsiella variicola" = "244366",
                    "Enterobacter hormaechei" = "158836", "Homo sapiens" = "9606")
if (!all(sources$organism %in% names(ORGANISM_TAXID))) stop("organism without a taxid here: ", paste(setdiff(sources$organism, names(ORGANISM_TAXID)), collapse = ", "))
src_taxid <- setNames(ORGANISM_TAXID[sources$organism], sources$source)

# NCBI family of each taxid, by E-utilities efetch (cached for the run).
fam_cache <- new.env()
efam <- function(ids) {
  unique(unlist(lapply(unique(ids), function(id) {
    if (is.null(fam_cache[[id]])) {
      Sys.sleep(0.4)
      x <- paste(readLines(sprintf("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id=%s", id), warn = FALSE), collapse = "")
      f <- regmatches(x, regexpr("<TaxId>[0-9]+</TaxId>\\s*<ScientificName>[^<]*</ScientificName>\\s*<Rank>family</Rank>", x))
      fam_cache[[id]] <- if (length(f)) sub(".*<TaxId>([0-9]+)</TaxId>.*", "\\1", f) else NA_character_
    }
    fam_cache[[id]]
  })))
}
# The unit a flag was reported as, in the scorer's own unit names.
unit_of <- function(fl, k, sp) {
  ids <- strsplit(fl$ncbi_taxid[k], ",")[[1]]
  plist_ids <- plist$match_taxid[plist$taxid %in% ids]
  unique(sp$unit[sp$taxid %in% plist_ids])[1]
}

res <- do.call(rbind, lapply(seq_len(nrow(truth)), function(i) {
  t <- truth[i, ]; s <- t$sample
  q <- qc[qc$sample == s, ]
  b <- q[q$tool == "Bracken", ]
  sp <- if (t$platform == "ILM" && nrow(b)) b else q[q$tool == "Kraken2", ]
  f <- do.call(rbind, strsplit(sp$field, "|", fixed = TRUE))
  sp <- data.frame(rank = f[, 1], taxid = f[, 2], n = as.numeric(sp$value))
  sp <- sp[sp$rank == "S" & sp$taxid != "9606", ]          # no target declared in layer D
  sp$unit <- ifelse(sp$taxid %in% groups$taxid, groups$group[match(sp$taxid, groups$taxid)], paste0("taxid:", sp$taxid))
  tot <- tapply(sp$n, sp$unit, sum); share <- tot / sum(tot)
  units <- names(share)[share >= 0.01]
  truth_taxa <- src_taxid[c(t$base, if (!is.na(t$contaminant)) t$contaminant)]
  listed <- function(u) unique(plist$taxid[plist$match_taxid %in% sp$taxid[sp$unit == u]])
  status <- function(u) { st <- plist$status[plist$taxid %in% listed(u)]; if (any(st == "established")) "established" else "putative" }
  expected <- units[vapply(units, function(u) length(listed(u)) > 0 && any(sp$taxid[sp$unit == u] %in% truth_taxa), TRUE)]
  any_listed <- units[vapply(units, function(u) length(listed(u)) > 0, TRUE)]
  fl <- flags[flags$sample == s & flags$result == "listed", ]
  fl$caveat[is.na(fl$caveat)] <- ""                  # read.delim turns blank cells into NA, and nzchar(NA) is TRUE
  # The report names a unit by its species name or group; the scorer by taxid
  # or group, so compare on the listed taxids each unit carries.
  key_report <- vapply(fl$ncbi_taxid, function(x) paste(sort(strsplit(x, ",")[[1]]), collapse = ","), "")
  key_exp    <- vapply(expected, function(u) paste(sort(listed(u)), collapse = ","), "")
  st_exp     <- vapply(expected, status, "")
  missing <- setdiff(key_exp, key_report)
  extra   <- setdiff(key_report, key_exp)
  main_u  <- names(share)[which.max(share)]
  caveat_wrong <- character(0)
  for (k in seq_len(nrow(fl))) {
    uk <- unit_of(fl, k, sp)
    if (identical(uk, main_u)) {                       # the main species carries no caveat
      if (nzchar(fl$caveat[k])) caveat_wrong <- c(caveat_wrong, key_report[k])
      next
    }
    fam_same <- length(intersect(efam(sp$taxid[sp$unit == uk]), efam(sp$taxid[sp$unit == main_u]))) > 0
    if (fam_same != nzchar(fl$caveat[k])) caveat_wrong <- c(caveat_wrong, key_report[k])
  }
  uncaveated_extra <- extra[!nzchar(fl$caveat[match(extra, key_report)])]
  wrong_status <- key_exp[key_exp %in% key_report & st_exp != fl$list_status[match(key_exp, key_report)]]
  data.frame(sample = s, expected = paste(key_exp, collapse = " | "), flagged = paste(key_report, collapse = " | "),
             missing = paste(missing, collapse = " | "), extra = paste(extra, collapse = " | "),
             extra_listed_not_truth = paste(setdiff(vapply(any_listed, function(u) paste(sort(listed(u)), collapse = ","), ""), key_exp), collapse = " | "),
             wrong_status = paste(wrong_status, collapse = " | "),
             uncaveated_extra = paste(uncaveated_extra, collapse = " | "),
             caveat_wrong = paste(caveat_wrong, collapse = " | "),
             result = if (!length(missing) && !length(uncaveated_extra) && !length(wrong_status) && !length(caveat_wrong)) "PASS" else "FAIL",
             stringsAsFactors = FALSE)
}))
write.table(res, args[5], sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("%s: %d of %d samples PASS\n", basename(run), sum(res$result == "PASS"), nrow(res)))
if (any(res$result == "FAIL")) print(res[res$result == "FAIL", c("sample", "missing", "uncaveated_extra", "wrong_status", "caveat_wrong")], row.names = FALSE)
