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
# CAVEAT        Round 1: ENA returns candidates in its own order, and only
#               the first CANDIDATES per stratum are drawn from, inside a
#               stated date window. One run per study then left pools as
#               small as 4 runs, and 19 of 40 picks were DDBJ runs.
# ROUND 2       (D51) Draws from the WHOLE pool: every matching run's
#               accession and study is fetched, studies are put in a seeded
#               random order, and one random run per study is checked for
#               eligibility until the stratum is full. Every study in the
#               round 1 panel file is excluded too: the D40 read-name markers
#               were designed after looking at those runs. The pool
#               definition and size are recorded in panel.tsv.
#
# Usage:  Rscript validation/select_panel.R <round> <out_panel.tsv>
#         round = 1 (40 runs), 2 (115 runs) or 3 (54 runs, D61: tests D58/D60;
#         reuses round 2's saved pools when copied next to the output)
#*******************************************************************************

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("Usage: select_panel.R <round> <out_panel.tsv>")
round   <- as.integer(args[1])
out_tsv <- args[2]

options(timeout = 900)   # whole-pool queries return up to ~360,000 rows; R's default is 60 s
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
stratum            | truth        | query                                                                                                         | from       | to         | r1 | r2 | r3 | sample_type | molecule
ILM_MISEQ          | illumina     | instrument_model="Illumina MiSeq" AND library_strategy="WGS" AND BACTERIA                                     | 2023-01-01 | 2024-12-31 | 2  | 8  | 3  | isolate     | dna
ILM_NEXTSEQ        | illumina     | instrument_model="NextSeq 2000" AND library_strategy="WGS" AND BACTERIA                                       | 2023-01-01 | 2024-12-31 | 2  | 8  | 3  | isolate     | dna
ILM_HISEQ          | illumina     | instrument_model="Illumina HiSeq 2500" AND library_strategy="WGS" AND BACTERIA                                | 2016-01-01 | 2019-12-31 | 2  | 8  | 3  | isolate     | dna
ILM_NOVASEQ6000    | illumina     | instrument_model="Illumina NovaSeq 6000" AND library_strategy="WGS" AND BACTERIA                              | 2022-01-01 | 2024-12-31 | 2  | 8  | 6  | isolate     | dna
ILM_NOVASEQX       | illumina     | instrument_model="Illumina NovaSeq X Plus" AND library_strategy="WGS" AND BACTERIA                            | 2024-01-01 | 2025-12-31 | 2  | 8  | 4  | isolate     | dna
ILM_METAGENOME     | illumina     | instrument_platform="ILLUMINA" AND library_source="METAGENOMIC" AND library_strategy="WGS"                    | 2023-01-01 | 2024-12-31 | 2  | 5  | 4  | metagenome  | dna
ILM_AMPLICON       | illumina     | instrument_platform="ILLUMINA" AND library_strategy="AMPLICON" AND library_source="METAGENOMIC"              | 2023-01-01 | 2024-12-31 | 2  | 5  | 0  | amplicon    | dna
ILM_RNASEQ         | illumina     | instrument_platform="ILLUMINA" AND library_strategy="RNA-Seq" AND library_source="TRANSCRIPTOMIC" AND BACTERIA | 2022-01-01 | 2024-12-31 | 2  | 5  | 4  | isolate     | rna
ILM_METATRANSCRIPT | illumina     | instrument_platform="ILLUMINA" AND library_source="METATRANSCRIPTOMIC"                                       | 2022-01-01 | 2024-12-31 | 2  | 5  | 4  | metagenome  | rna
ONT_OLD            | ont          | instrument_platform="OXFORD_NANOPORE" AND library_strategy="WGS" AND BACTERIA                                 | 2019-01-01 | 2021-12-31 | 4  | 10 | 2  | isolate     | dna
ONT_R10_ERA        | ont          | instrument_platform="OXFORD_NANOPORE" AND library_strategy="WGS" AND BACTERIA                                 | 2024-06-01 | 2025-12-31 | 4  | 10 | 2  | isolate     | dna
ONT_RNA            | ont          | instrument_platform="OXFORD_NANOPORE" AND library_source="TRANSCRIPTOMIC"                                    | 2021-01-01 | 2025-12-31 | 2  | 5  | 0  |             | rna
ONT_METAGENOME     | ont          | instrument_platform="OXFORD_NANOPORE" AND library_source="METAGENOMIC" AND library_strategy="WGS"            | 2022-01-01 | 2025-12-31 | 2  | 5  | 0  | metagenome  | dna
PB_HIFI            | pacbio       | instrument_model="Revio"                                                                                      | 2023-06-01 | 2025-12-31 | 2  | 5  | 2  |             |
PB_CLR             | pacbio       | instrument_model="PacBio RS II"                                                                               | 2015-01-01 | 2019-12-31 | 2  | 5  | 2  |             |
DNBSEQ             | dnbseq       | instrument_platform="BGISEQ" AND library_strategy="WGS" AND BACTERIA                                          | 2021-01-01 | 2025-12-31 | 2  | 5  | 8  | isolate     | dna
ELEMENT            | element      | instrument_platform="ELEMENT"                                                                                 | 2023-01-01 | 2025-12-31 | 2  | 5  | 4  |             |
ION_TORRENT        | ion_torrent  | instrument_platform="ION_TORRENT" AND BACTERIA                                                                | 2018-01-01 | 2025-12-31 | 2  | 5  | 3  | isolate     | dna
')

FIELDS <- c("run_accession", "study_accession", "sample_accession", "instrument_platform", "instrument_model",
            "library_strategy", "library_source", "library_layout", "tax_id", "scientific_name",
            "first_public", "read_count", "fastq_ftp", "fastq_bytes")

ena_search <- function(query, n, fields = FIELDS) {
  url <- paste0("https://www.ebi.ac.uk/ena/portal/api/search?result=read_run",
                "&query=", URLencode(query, reserved = TRUE),
                "&fields=", paste(fields, collapse = ","),
                "&limit=", n, "&format=tsv")
  # ENA's API sometimes answers 500 to a large query; try again before giving up.
  for (attempt in 1:4) {
    res <- tryCatch(read.delim(url(url), colClasses = "character", na.strings = "", quote = ""),
                    error = function(e) e)
    if (!inherits(res, "error")) return(res)
    message(sprintf("  ENA query failed (attempt %d): %s", attempt, conditionMessage(res)))
    Sys.sleep(10 * attempt)
  }
  stop("ENA query failed 4 times: ", query)
}
# A run is usable if its FASTQ is downloadable as 1 file (single) or 2 (paired;
# a third file of unpaired reads is allowed and dropped later).
eligible <- function(cand) {
  n_files <- lengths(strsplit(ifelse(is.na(cand$fastq_ftp), "", cand$fastq_ftp), ";"))
  # %in% rather than ==, so a blank field in ENA's record (NA) means "not
  # eligible" instead of an NA that stops the script.
  !is.na(cand$fastq_ftp) & !(cand$tax_id %in% "9606") &
    ((cand$library_layout %in% "SINGLE" & n_files == 1) | (cand$library_layout %in% "PAIRED" & n_files %in% 2:3))
}

# Round 2 and later: exclude every study in every EARLIER round's panel file,
# picked or backup (round 3, D61: rounds 1 and 2).
if (round >= 2) {
  here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
  for (r in seq_len(round - 1)) {
    prev <- read.delim(file.path(here, "outcomes", paste0("round", r), "panel.tsv"), colClasses = "character")
    st   <- unique(prev$study_accession)
    message(sprintf("Excluding %d round %d studies", length(st), r))
    CALIBRATION_STUDIES <- c(CALIBRATION_STUDIES, st)
  }
  CALIBRATION_STUDIES <- unique(CALIBRATION_STUDIES)
}

n_col <- paste0("r", round)
panel <- list()
for (k in seq_len(nrow(strata))) {
  s <- strata[k, ]
  q <- sprintf("%s AND %s AND first_public>=%s AND first_public<=%s",
               gsub("BACTERIA", BACTERIA, s$query, fixed = TRUE), NOT_HUMAN, s$from, s$to)
  if (s[[n_col]] == 0) next          # stratum not in this round
  want <- s[[n_col]] + BACKUPS
  if (round == 1) {
    cand <- ena_search(q, CANDIDATES)
    keep <- eligible(cand) & !cand$study_accession %in% CALIBRATION_STUDIES
    cand <- cand[keep, ]
    # One run per study, so one large study cannot dominate a stratum.
    cand <- cand[!duplicated(cand$study_accession), ]
    set.seed(SEED + k)
    pick <- cand[sample(nrow(cand), min(want, nrow(cand))), ]
    pool_text <- sprintf("first %d candidates; %d eligible after exclusions", CANDIDATES, nrow(cand))
  } else {
    # Whole pool, two columns only, fetched one calendar year at a time
    # (smaller requests, which ENA's API answers more reliably; together they
    # are the same pool). Each pool is saved next to the panel, so a re-run
    # reuses it and the panel records exactly what it was drawn from.
    pool_dir <- paste0(out_tsv, ".pools")
    dir.create(pool_dir, showWarnings = FALSE)
    pool_file <- file.path(pool_dir, paste0(s$stratum, ".tsv"))
    if (file.exists(pool_file)) {
      all_runs <- read.delim(pool_file, colClasses = "character")
    } else {
      years <- as.integer(substr(s$from, 1, 4)):as.integer(substr(s$to, 1, 4))
      all_runs <- do.call(rbind, lapply(years, function(y) {
        lo <- max(s$from, sprintf("%d-01-01", y)); hi <- min(s$to, sprintf("%d-12-31", y))
        qy <- sprintf("%s AND %s AND first_public>=%s AND first_public<=%s",
                      gsub("BACTERIA", BACTERIA, s$query, fixed = TRUE), NOT_HUMAN, lo, hi)
        ena_search(qy, 0, c("run_accession", "study_accession"))
      }))
      write.table(all_runs, pool_file, sep = "\t", quote = FALSE, row.names = FALSE)
    }
    all_runs <- all_runs[!all_runs$study_accession %in% CALIBRATION_STUDIES, ]
    # Sorted first, so the seeded shuffle depends only on the pool's contents,
    # not on the order ENA happened to return it in.
    all_runs <- all_runs[order(all_runs$study_accession, all_runs$run_accession), ]
    studies  <- sort(unique(all_runs$study_accession))
    set.seed(SEED + 1000 * round + k)
    studies  <- studies[sample(length(studies))]
    pick <- list(); tried <- 0L
    for (st in studies) {
      if (length(pick) == want) break
      runs <- all_runs$run_accession[all_runs$study_accession == st]
      # Up to 3 random runs per study: some runs lack FASTQ or have odd files.
      for (r in runs[sample(length(runs), min(3L, length(runs)))]) {
        tried <- tried + 1L
        d <- ena_search(sprintf('run_accession="%s"', r), 1)
        if (nrow(d) == 1 && eligible(d)) { pick[[length(pick) + 1]] <- d; break }
      }
    }
    pick <- if (length(pick)) do.call(rbind, pick)
            else as.data.frame(setNames(replicate(length(FIELDS), character(0), simplify = FALSE), FIELDS))
    pool_text <- sprintf("whole pool: %d runs in %d studies after exclusions; one random run per study, %d runs checked",
                         nrow(all_runs), length(studies), tried)
  }
  if (nrow(pick) < s[[n_col]]) warning(sprintf("%s: only %d eligible runs", s$stratum, nrow(pick)))
  pick$stratum     <- s$stratum
  pick$truth       <- s$truth
  pick$role        <- c(rep("panel", min(s[[n_col]], nrow(pick))), rep("backup", max(0, nrow(pick) - s[[n_col]])))
  pick$sample_type <- s$sample_type
  pick$molecule    <- s$molecule
  pick$pool        <- sprintf("%s; first_public %s..%s; %s", s$query, s$from, s$to, pool_text)
  panel[[k]] <- pick
  message(sprintf("%-19s %s -> %d picked", s$stratum, pool_text, nrow(pick)))
}
panel <- do.call(rbind, panel)
# Paired runs on ENA sometimes carry a third file of unpaired reads: keep _1 and _2 only.
panel$fastq_ftp <- vapply(strsplit(panel$fastq_ftp, ";"), function(f) {
  if (length(f) == 3) paste(grep("_[12]\\.fastq\\.gz$", f, value = TRUE), collapse = ";") else paste(f, collapse = ";")
}, "")
# Archive the run was submitted to (DRR = DDBJ, ERR = ENA, SRR = NCBI SRA):
# DDBJ-served files can lose their original read names (open question A1).
panel$archive <- c(D = "DDBJ", E = "ENA", S = "SRA")[substr(panel$run_accession, 1, 1)]
cols <- c("stratum", "role", "run_accession", "archive", "truth", "library_layout", "instrument_platform", "instrument_model",
          "library_strategy", "library_source", "scientific_name", "tax_id", "study_accession",
          "first_public", "read_count", "sample_type", "molecule", "fastq_ftp", "pool")
write.table(panel[, cols], out_tsv, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
message(sprintf("\n%d panel runs + %d backups written to %s", sum(panel$role == "panel"),
                sum(panel$role == "backup"), out_tsv))
