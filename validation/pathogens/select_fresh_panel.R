#!/usr/bin/env Rscript
# =============================================================================
# select_fresh_panel.R  (D96, pass mark P7)
# Draws the fresh-data panel for the pathogen flags exactly as pre-registered
# in D96: ENA read runs (WGS, GENOMIC, Illumina, paired, under Bacteria,
# species-rank tax_id, 500,000-5,000,000 reads, first public in June 2024),
# one run per study (lowest accession), no study used in any earlier panel;
# split by whether tax_id is on the resolved list; 10 drawn from each arm with
# set.seed(97), the rest in drawn order as backups.
#
# Usage:  select_fresh_panel.R <pathogen_list_dir> <out.tsv>
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("Usage: select_fresh_panel.R <pathogen_list_dir> <out.tsv>")
rd <- function(f) read.delim(f, colClasses = "character", na.strings = "", quote = "")
pl  <- rd(file.path(args[1], "pathogen_list.tsv"))
fam <- rd(file.path(args[1], "pathogen_list_families.tsv"))   # every bacterial species-rank taxid

query <- paste('library_strategy="WGS"', 'library_source="GENOMIC"', 'instrument_platform="ILLUMINA"',
               'library_layout="PAIRED"', "tax_tree(2)", "read_count>=500000", "read_count<=5000000",
               "first_public>=2024-06-01", "first_public<=2024-06-30", sep = " AND ")
url <- paste0("https://www.ebi.ac.uk/ena/portal/api/search?result=read_run&format=tsv&limit=0",
              "&fields=run_accession,study_accession,tax_id,scientific_name,read_count,first_public,fastq_ftp",
              "&query=", URLencode(query, reserved = TRUE))
f <- tempfile(fileext = ".tsv")
download.file(url, f, quiet = TRUE)
pool <- rd(f)
n_all <- nrow(pool)

# Species-rank bacterial taxa only (NCBI taxonomy of the list build).
pool <- pool[pool$tax_id %in% fam$taxid & !is.na(pool$fastq_ftp), ]
# One run per study: the lowest accession.
pool <- pool[order(pool$study_accession, pool$run_accession), ]
pool <- pool[!duplicated(pool$study_accession), ]
# No study from an earlier panel.
here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
panels <- list.files(file.path(here, "..", "outcomes"), pattern = "^panel\\.tsv$", recursive = TRUE, full.names = TRUE)
used <- unique(unlist(lapply(panels, function(p) { x <- rd(p); x$study_accession })))
pool <- pool[!pool$study_accession %in% used, ]

pool$arm <- ifelse(pool$tax_id %in% pl$match_taxid, "listed", "unlisted")
set.seed(97)
out <- do.call(rbind, lapply(c("listed", "unlisted"), function(a) {
  p <- pool[pool$arm == a, ]
  p <- p[sample(nrow(p)), ]
  p$order <- seq_len(nrow(p))
  p$role  <- ifelse(p$order <= 10, "panel", "backup")
  p
}))
write.table(out, args[2], sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("ENA returned %d runs; %d studies eligible after the rules (%d listed, %d unlisted); %d earlier panel studies excluded\n",
            n_all, nrow(pool), sum(pool$arm == "listed"), sum(pool$arm == "unlisted"), length(used)))
