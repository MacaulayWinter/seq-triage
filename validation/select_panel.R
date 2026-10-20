#!/usr/bin/env Rscript
#*******************************************************************************
# select_panel.R: choose the validation panel from ENA, by rule (CLAUDE.md D37)
#
# WHAT IT DOES  For each category (stratum) below, asks ENA's public API for
#               candidate runs, removes anything excluded by the rules, and
#               picks runs at random with a fixed seed. Writes panel.tsv: one
#               row per run, with the TRUE platform and layout from ENA.
# RULES         1. No run from a study used to calibrate detection.
#               2. No human samples (tax_id 9606).
#               3. Runs must have downloadable FASTQ (1 file single, 2 paired).
#               4. Same inputs + same seed = same panel, for anyone.
# CAVEAT        ENA returns candidates in its own order, and only the first
#               CANDIDATES per stratum are drawn from, inside a stated date
#               window. The panel is random within that pool, not across all
#               of ENA. The pool definition is recorded in panel.tsv.
#
# Usage:  Rscript validation/select_panel.R <round> <out_panel.tsv>
#         round = 1 (40 runs) or 2 (~150 runs)
#*******************************************************************************

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("Usage: select_panel.R <round> <out_panel.tsv>")
round   <- as.integer(args[1])
out_tsv <- args[2]

SEED       <- 20261020
CANDIDATES <- 500     # runs fetched per stratum before the random pick
BACKUPS    <- 2       # extra runs per stratum, used only if a download fails

# Studies whose data tuned the detection limits or was inspected while
# choosing test data (CLAUDE.md "Test data survey", D31, D33). Looked up from
# each run accession on ENA (2026-10-20). Survey files without an ENA
# accession (bacass A1403KPN, nanoseq NA12878, simulated B. fragilis) cannot
# be excluded by study; they are not in ENA's index anyway.
CALIBRATION_STUDIES <- c(
  "PRJEB2510",    # ERR044595  HiSeq 2000
  "PRJEB2489",    # ERR064912  HiSeq 2000
  "PRJNA299729",  # SRR2838702 HiSeq 2500 + ERR3772599's Portiera
  "PRJEB44507",   # ERR5766176/80 MiSeq
  "PRJEB29152",   # ERR3201952 ONT
  "PRJEB35780",   # ERR3772599 ONT
  "PRJNA586863",  # SRR10382244 PacBio HiFi
  "PRJEB37886",   # ERR5069949 SARS-CoV-2 (nf-core test data)
  "PRJNA1287288", # ECE107 NovaSeq + R10.4.1
  "PRJNA607328",  # SRR30615993, inspected for R10 metadata
  "PRJNA1206708", # SRR30655751, inspected
  "PRJDB5929",    # DRR293631, inspected
  "PRJNA1165667"  # SRR30810011, inspected
)
NOT_HUMAN <- "NOT tax_eq(9606)"
BACTERIA  <- "tax_tree(2)"

# Strata: name, ENA query, date window, runs in round 1 and round 2, and the
# metadata the samplesheet will declare (platform is never declared: the
# reads must decide). truth = what ENA says the platform is.
strata <- read.table(header = TRUE, sep = "|", strip.white = TRUE, stringsAsFactors = FALSE, quote = "", text = '
stratum            | truth        | query                                                                                                         | from       | to         | r1 | r2 | sample_type | molecule
ILM_MISEQ          | illumina     | instrument_model="Illumina MiSeq" AND library_strategy="WGS" AND BACTERIA                                     | 2023-01-01 | 2024-12-31 | 2  | 8  | isolate     | dna
ILM_NEXTSEQ        | illumina     | instrument_model="NextSeq 2000" AND library_strategy="WGS" AND BACTERIA                                       | 2023-01-01 | 2024-12-31 | 2  | 8  | isolate     | dna
ILM_HISEQ          | illumina     | instrument_model="Illumina HiSeq 2500" AND library_strategy="WGS" AND BACTERIA                                | 2016-01-01 | 2019-12-31 | 2  | 8  | isolate     | dna
ILM_NOVASEQ6000    | illumina     | instrument_model="Illumina NovaSeq 6000" AND library_strategy="WGS" AND BACTERIA                              | 2022-01-01 | 2024-12-31 | 2  | 8  | isolate     | dna
ILM_NOVASEQX       | illumina     | instrument_model="Illumina NovaSeq X Plus" AND library_strategy="WGS" AND BACTERIA                            | 2024-01-01 | 2025-12-31 | 2  | 8  | isolate     | dna
ILM_METAGENOME     | illumina     | instrument_platform="ILLUMINA" AND library_source="METAGENOMIC" AND library_strategy="WGS"                    | 2023-01-01 | 2024-12-31 | 2  | 5  | metagenome  | dna
ILM_AMPLICON       | illumina     | instrument_platform="ILLUMINA" AND library_strategy="AMPLICON" AND library_source="METAGENOMIC"              | 2023-01-01 | 2024-12-31 | 2  | 5  | amplicon    | dna
ILM_RNASEQ         | illumina     | instrument_platform="ILLUMINA" AND library_strategy="RNA-Seq" AND library_source="TRANSCRIPTOMIC" AND BACTERIA | 2022-01-01 | 2024-12-31 | 2  | 5  | isolate     | rna
ILM_METATRANSCRIPT | illumina     | instrument_platform="ILLUMINA" AND library_source="METATRANSCRIPTOMIC"                                       | 2022-01-01 | 2024-12-31 | 2  | 5  | metagenome  | rna
ONT_OLD            | ont          | instrument_platform="OXFORD_NANOPORE" AND library_strategy="WGS" AND BACTERIA                                 | 2019-01-01 | 2021-12-31 | 4  | 10 | isolate     | dna
ONT_R10_ERA        | ont          | instrument_platform="OXFORD_NANOPORE" AND library_strategy="WGS" AND BACTERIA                                 | 2024-06-01 | 2025-12-31 | 4  | 10 | isolate     | dna
ONT_RNA            | ont          | instrument_platform="OXFORD_NANOPORE" AND library_source="TRANSCRIPTOMIC"                                    | 2021-01-01 | 2025-12-31 | 2  | 5  |             | rna
ONT_METAGENOME     | ont          | instrument_platform="OXFORD_NANOPORE" AND library_source="METAGENOMIC" AND library_strategy="WGS"            | 2022-01-01 | 2025-12-31 | 2  | 5  | metagenome  | dna
PB_HIFI            | pacbio       | instrument_model="Revio"                                                                                      | 2023-06-01 | 2025-12-31 | 2  | 5  |             |
PB_CLR             | pacbio       | instrument_model="PacBio RS II"                                                                               | 2015-01-01 | 2019-12-31 | 2  | 5  |             |
DNBSEQ             | dnbseq       | instrument_platform="BGISEQ" AND library_strategy="WGS" AND BACTERIA                                          | 2021-01-01 | 2025-12-31 | 2  | 5  | isolate     | dna
ELEMENT            | element      | instrument_platform="ELEMENT"                                                                                 | 2023-01-01 | 2025-12-31 | 2  | 5  |             |
ION_TORRENT        | ion_torrent  | instrument_platform="ION_TORRENT" AND BACTERIA                                                                | 2018-01-01 | 2025-12-31 | 2  | 5  | isolate     | dna
')

FIELDS <- c("run_accession", "study_accession", "sample_accession", "instrument_platform", "instrument_model",
            "library_strategy", "library_source", "library_layout", "tax_id", "scientific_name",
            "first_public", "read_count", "fastq_ftp", "fastq_bytes")

ena_search <- function(query, n) {
  url <- paste0("https://www.ebi.ac.uk/ena/portal/api/search?result=read_run",
                "&query=", URLencode(query, reserved = TRUE),
                "&fields=", paste(FIELDS, collapse = ","),
                "&limit=", n, "&format=tsv")
  read.delim(url(url), colClasses = "character", na.strings = "", quote = "")
}

n_col <- if (round == 1) "r1" else "r2"
panel <- list()
for (k in seq_len(nrow(strata))) {
  s <- strata[k, ]
  q <- sprintf("%s AND %s AND first_public>=%s AND first_public<=%s",
               gsub("BACTERIA", BACTERIA, s$query, fixed = TRUE), NOT_HUMAN, s$from, s$to)
  cand <- ena_search(q, CANDIDATES)
  n_files <- lengths(strsplit(ifelse(is.na(cand$fastq_ftp), "", cand$fastq_ftp), ";"))
  keep <- !is.na(cand$fastq_ftp) & !cand$study_accession %in% CALIBRATION_STUDIES &
          cand$tax_id != "9606" &
          ((cand$library_layout == "SINGLE" & n_files == 1) | (cand$library_layout == "PAIRED" & n_files %in% 2:3))
  cand <- cand[keep, ]
  # One run per study, so one large study cannot dominate a stratum.
  cand <- cand[!duplicated(cand$study_accession), ]
  want <- s[[n_col]] + BACKUPS
  set.seed(SEED + k)
  pick <- cand[sample(nrow(cand), min(want, nrow(cand))), ]
  if (nrow(pick) < s[[n_col]]) warning(sprintf("%s: only %d eligible runs", s$stratum, nrow(pick)))
  pick$stratum     <- s$stratum
  pick$truth       <- s$truth
  pick$role        <- c(rep("panel", min(s[[n_col]], nrow(pick))), rep("backup", max(0, nrow(pick) - s[[n_col]])))
  pick$sample_type <- s$sample_type
  pick$molecule    <- s$molecule
  pick$pool        <- sprintf("%s; first_public %s..%s; first %d candidates; %d eligible after exclusions",
                              s$query, s$from, s$to, CANDIDATES, nrow(cand))
  panel[[k]] <- pick
  message(sprintf("%-19s %3d eligible -> %d picked", s$stratum, nrow(cand), nrow(pick)))
}
panel <- do.call(rbind, panel)
# Paired runs on ENA sometimes carry a third file of unpaired reads: keep _1 and _2 only.
panel$fastq_ftp <- vapply(strsplit(panel$fastq_ftp, ";"), function(f) {
  if (length(f) == 3) paste(grep("_[12]\\.fastq\\.gz$", f, value = TRUE), collapse = ";") else paste(f, collapse = ";")
}, "")
cols <- c("stratum", "role", "run_accession", "truth", "library_layout", "instrument_platform", "instrument_model",
          "library_strategy", "library_source", "scientific_name", "tax_id", "study_accession",
          "first_public", "read_count", "sample_type", "molecule", "fastq_ftp", "pool")
write.table(panel[, cols], out_tsv, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
message(sprintf("\n%d panel runs + %d backups written to %s", sum(panel$role == "panel"),
                sum(panel$role == "backup"), out_tsv))
