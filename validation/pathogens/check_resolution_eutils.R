#!/usr/bin/env Rscript
# =============================================================================
# check_resolution_eutils.R  (D96, pass mark P1)
# Checks pathogen_list.R's name -> taxid step by an independent route: NCBI
# E-utilities esearch on the taxonomy database, for 60 resolved names drawn
# with set.seed(96). Agreement = esearch returns exactly the same taxid.
#
# Usage:  check_resolution_eutils.R <pathogen_list.tsv> <out.tsv>
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("Usage: check_resolution_eutils.R <pathogen_list.tsv> <out.tsv>")
pl <- read.delim(args[1], colClasses = "character", quote = "")
pl <- pl[pl$match_taxid == pl$taxid, ]
names_one <- unique(data.frame(species = unlist(strsplit(pl$species, " / ")),
                               taxid = rep(pl$taxid, lengths(strsplit(pl$species, " / "))),
                               resolved_by = rep(pl$resolved_by, lengths(strsplit(pl$species, " / ")))))
set.seed(96)
pick <- names_one[sort(sample(nrow(names_one), 60)), ]
pick$esearch <- vapply(pick$species, function(n) {
  Sys.sleep(0.4)                                   # NCBI: at most 3 requests a second without a key
  u <- sprintf("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=taxonomy&term=%s", URLencode(n, reserved = TRUE))
  x <- paste(readLines(u, warn = FALSE), collapse = "")
  ids <- regmatches(x, gregexpr("<Id>[0-9]+</Id>", x))[[1]]
  paste(gsub("\\D", "", ids), collapse = ",")
}, "")
pick$agree <- pick$esearch == pick$taxid
write.table(pick, args[2], sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("P1: %d of %d agree\n", sum(pick$agree), nrow(pick)))
if (!all(pick$agree)) print(pick[!pick$agree, ], row.names = FALSE)
