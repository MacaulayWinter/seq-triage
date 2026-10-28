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
# list. PASS when the report flags exactly those units, each with the right
# status, and nothing else.
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
  # The report names a unit by its species name or group; the scorer by taxid
  # or group, so compare on the listed taxids each unit carries.
  key_report <- vapply(fl$ncbi_taxid, function(x) paste(sort(strsplit(x, ",")[[1]]), collapse = ","), "")
  key_exp    <- vapply(expected, function(u) paste(sort(listed(u)), collapse = ","), "")
  st_exp     <- vapply(expected, status, "")
  missing <- setdiff(key_exp, key_report)
  extra   <- setdiff(key_report, key_exp)
  wrong_status <- key_exp[key_exp %in% key_report & st_exp != fl$list_status[match(key_exp, key_report)]]
  data.frame(sample = s, expected = paste(key_exp, collapse = " | "), flagged = paste(key_report, collapse = " | "),
             missing = paste(missing, collapse = " | "), extra = paste(extra, collapse = " | "),
             extra_listed_not_truth = paste(setdiff(vapply(any_listed, function(u) paste(sort(listed(u)), collapse = ","), ""), key_exp), collapse = " | "),
             wrong_status = paste(wrong_status, collapse = " | "),
             result = if (!length(missing) && !length(extra) && !length(wrong_status)) "PASS" else "FAIL",
             stringsAsFactors = FALSE)
}))
write.table(res, args[5], sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("%s: %d of %d samples PASS\n", basename(run), sum(res$result == "PASS"), nrow(res)))
if (any(res$result == "FAIL")) print(res[res$result == "FAIL", c("sample", "expected", "flagged", "missing", "extra")], row.names = FALSE)
