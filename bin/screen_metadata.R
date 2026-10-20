#!/usr/bin/env Rscript
#*******************************************************************************
# screen_metadata.R: seq-triage metadata screen
#
# WHAT IT DOES  Reads the raw samplesheet and flags inconsistencies that a
#               per-row schema cannot see: IDs that collide once normalised,
#               IDs that look like file names, files that don't match their
#               sample, shared files, and platform values that needed
#               normalising. Writes one row per sample.
# WHY           The schema (assets/schema_input.json) already rejected
#               everything unambiguously wrong. What reaches this script is
#               valid but possibly suspicious, so every finding here is a
#               WARNING carried into the final report, never an error.
# RULE          Never rewrite a value silently. Each normalised field keeps
#               its raw value alongside, with the tier that explains how it
#               was matched (CLAUDE.md D6, D7).
# PROVENANCE    blank_na(), normalise_identifier() and the file-suffix regex
#               are adapted from the author's R_Bible (Parts 01 and 23).
#               Base R only, so the container needs nothing but R.
#
# Usage:  screen_metadata.R <samplesheet.csv> <output.tsv>
#*******************************************************************************

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: screen_metadata.R <samplesheet.csv> <output.tsv>")
}
in_csv  <- args[1]
out_tsv <- args[2]


# Helpers ----------------------------------------------------------------------

# WHAT IT DOES  TRUE where a value is one of the many ways people write "empty".
# PROVENANCE    R_Bible 1.7 blank_na(). The token list matches the schema's
#               accepted empty tokens, so anything the schema let through as
#               "empty" is recognised here.
EMPTY_TOKENS <- c("na", "n/a", "nan", "null", "none", "missing", ".", "-", "--")
is_empty_token <- function(x) !is.na(x) & tolower(trimws(x)) %in% EMPTY_TOKENS

# WHAT IT DOES  Builds a join key: upper case, whitespace and _ ( ) . -
#               removed. "S_1", "s-1" and "S 1" all become "S1".
# PROVENANCE    R_Bible 23 normalise_identifier(), without prefix stripping.
# CAUTION       A key, never a display value. Lossy by design, which is
#               exactly what makes collisions detectable.
normalise_identifier <- function(x) {
  gsub("[[:space:]_().-]+", "", toupper(as.character(x)))
}

# Sequencing-file noise: Illumina's _S<n>, _L<lane>, _R<mate>, _001 and the
# FASTQ extension. PROVENANCE: R_Bible 1.10 diagnose_id_mismatch().
FILE_NOISE <- "(_S[0-9]+)?(_L[0-9]{3})?(_R[12])?(_001)?(\\.(fastq|fq)(\\.gz)?)?$"

# WHAT IT DOES  File name without directory, extension or Illumina suffixes,
#               e.g. ".../ISO_01_S3_L001_R1_001.fastq.gz" -> "ISO_01".
strip_file_noise <- function(path) sub(FILE_NOISE, "", basename(path))

# WHAT IT DOES  Which mate a file name claims to be: "R1", "R2" or NA.
#               Matches _R1 / _R2 / _1 / _2 (optionally _001) right before
#               the FASTQ extension, as in "X_R1_001.fastq.gz" or "SRR1_2.fq.gz".
mate_from_name <- function(path) {
  pat <- "^.*[_.]R?([12])(_001)?\\.f(ast)?q(\\.gz)?$"
  b   <- basename(path)
  ifelse(!is.na(b) & grepl(pat, b), paste0("R", sub(pat, "\\1", b)), NA_character_)
}

# WHAT IT DOES  Normalises a controlled-vocabulary field and records HOW each
#               value was matched, in the spirit of R_Bible 23's match tiers.
# RETURNS       data.frame(value, tier) where tier is one of
#                 blank        cell was empty
#                 empty_token  "N/A", "none", "." ... treated as not declared
#                 exact        already the canonical form
#                 normalised   mapped from a variant (case or synonym)
normalise_vocab <- function(raw, synonyms) {
  key   <- tolower(trimws(raw))
  value <- unname(synonyms[key])
  tier  <- ifelse(is.na(raw) | trimws(raw) == "", "blank",
           ifelse(is_empty_token(raw),          "empty_token",
           ifelse(!is.na(value) & raw == value,  "exact", "normalised")))
  value[tier %in% c("blank", "empty_token")] <- NA_character_
  data.frame(value = value, tier = tier, stringsAsFactors = FALSE)
}

PLATFORM_SYNONYMS <- c(illumina = "illumina", ont = "ont", nanopore = "ont")
TYPE_SYNONYMS     <- c(isolate = "isolate", metagenome = "metagenome",
                       amplicon = "amplicon", unknown = "unknown")
MOLECULE_SYNONYMS <- c(dna = "dna", rna = "rna")


# Read ---------------------------------------------------------------------------

# Everything as character, and na.strings = "" only: a literal "NA" must stay
# visible as text so it can be flagged, instead of silently becoming missing.
sheet <- read.csv(in_csv, colClasses = "character", na.strings = "",
                  check.names = FALSE, strip.white = FALSE)
for (col in c("fastq_2", "platform", "sample_type", "molecule")) {
  if (!col %in% names(sheet)) sheet[[col]] <- NA_character_
}
n <- nrow(sheet)


# Screen -----------------------------------------------------------------------

# Each check appends "CODE: message" strings to that sample's list of flags.
flags <- vector("list", n)
add_flag <- function(rows, code, msg) {
  for (k in seq_along(rows)) {
    flags[[rows[k]]] <<- c(flags[[rows[k]]], paste0(code, ": ", msg[k]))
  }
}

plat <- normalise_vocab(sheet$platform,    PLATFORM_SYNONYMS)
type <- normalise_vocab(sheet$sample_type, TYPE_SYNONYMS)
mol  <- normalise_vocab(sheet$molecule,    MOLECULE_SYNONYMS)

# Values that needed interpreting: report them, never hide them.
for (f in list(list(name = "platform",    raw = sheet$platform,    norm = plat),
               list(name = "sample_type", raw = sheet$sample_type, norm = type),
               list(name = "molecule",    raw = sheet$molecule,    norm = mol))) {
  r <- which(f$norm$tier == "empty_token")
  add_flag(r, "EMPTY_TOKEN", sprintf("%s '%s' treated as not declared", f$name, f$raw[r]))
  r <- which(f$norm$tier == "normalised")
  add_flag(r, "NORMALISED", sprintf("%s '%s' read as '%s'", f$name, f$raw[r], f$norm$value[r]))
}

# Stray whitespace in any field. nf-schema trims it silently before
# validation (tested 2026-10-20), so this flag is the only record that the
# value was changed. It is also a sign the sheet was hand-edited or exported
# from a spreadsheet.
for (col in c("sample", "fastq_1", "fastq_2", "platform", "sample_type", "molecule")) {
  v <- sheet[[col]]
  r <- which(!is.na(v) & v != trimws(v))
  add_flag(r, "WHITESPACE", sprintf("%s has leading or trailing spaces", rep(col, length(r))))
}

# ID collisions: distinct IDs that become identical once normalised. The
# schema already rejected exact duplicates, so any hit here is S1 vs s1.
id_key <- normalise_identifier(sheet$sample)
for (k in unique(id_key[duplicated(id_key)])) {
  r <- which(id_key == k)
  for (i in r) {
    add_flag(i, "ID_COLLISION",
             sprintf("'%s' collides with %s after normalising case and punctuation",
                     sheet$sample[i], paste0("'", sheet$sample[setdiff(r, i)], "'", collapse = ", ")))
  }
}

# IDs carrying sequencing-file noise were probably copied from a file name.
r <- which(sheet$sample != strip_file_noise(sheet$sample))
add_flag(r, "ID_LOOKS_LIKE_FILENAME",
         sprintf("'%s' ends like a file name; did you mean '%s'?",
                 sheet$sample[r], strip_file_noise(sheet$sample[r])))

# File names that don't contain their sample ID. A warning, not an error:
# public and facility data often use accession numbers as file names.
for (col in c("fastq_1", "fastq_2")) {
  v   <- sheet[[col]]
  has <- !is.na(v)
  # Strip file noise from BOTH sides, so an ID like "X_R1_001" (already
  # flagged above) does not also raise a second flag for the same mistake.
  id_core <- normalise_identifier(strip_file_noise(sheet$sample))
  f_core  <- normalise_identifier(ifelse(has, strip_file_noise(v), ""))
  hit <- has & !mapply(grepl, id_core, f_core, fixed = TRUE)
  r <- which(hit)
  add_flag(r, "BASENAME_MISMATCH",
           sprintf("%s '%s' does not contain the sample ID", rep(col, length(r)), basename(v[r])))
}

# Mate naming. Paired rows: mates named the wrong way round, using either
# _R2 or bare _2. Single-file rows: "swapped" means nothing, so flag only an
# EXPLICIT _R2 (a lone second mate suggests R1 is missing). A bare _2 is too
# often just part of a name ("test_2.fastq.gz"): the first end-to-end run
# showed exactly that false positive on real ONT test data.
m1       <- mate_from_name(sheet$fastq_1)
m2       <- mate_from_name(sheet$fastq_2)
paired   <- !is.na(sheet$fastq_2)
explicit <- grepl("[_.]R[12](_001)?\\.f(ast)?q(\\.gz)?$", basename(sheet$fastq_1))
r <- which(paired & ((m1 %in% "R2") | (m2 %in% "R1")))
add_flag(r, "MATES_POSSIBLY_SWAPPED", rep("fastq_1 is named like R2 or fastq_2 like R1", length(r)))
r <- which(!paired & m1 %in% "R2" & explicit)
add_flag(r, "SINGLE_FILE_NAMED_R2", rep("only one file given and it is named like a second mate (_R2); check R1 is not missing", length(r)))

# The same file used by more than one sample.
files  <- data.frame(row = c(seq_len(n), seq_len(n)),
                     path = c(sheet$fastq_1, sheet$fastq_2), stringsAsFactors = FALSE)
files  <- files[!is.na(files$path), ]
shared <- unique(files$path[duplicated(files$path)])
for (p in shared) {
  r <- unique(files$row[files$path == p])
  for (i in r) {
    add_flag(i, "FASTQ_SHARED",
             sprintf("'%s' is also used by %s", basename(p),
                     paste0("'", sheet$sample[setdiff(r, i)], "'", collapse = ", ")))
  }
}

# Declared ONT with two files: nanopore reads are not paired-end.
r <- which(plat$value %in% "ont" & !is.na(sheet$fastq_2))
add_flag(r, "LAYOUT_CONFLICT", rep("platform declared as ONT but two FASTQ files given", length(r)))


# Write ------------------------------------------------------------------------

# `sample` is the key the pipeline joins on. nf-schema trims whitespace before
# building meta.id (tested), so the key is trimmed here to match. Otherwise a
# sample written as " B" would silently drop out of the join. The WHITESPACE
# flag above still records that the raw value had spaces.
out <- data.frame(
  sample                 = trimws(sheet$sample),
  declared_platform_raw  = sheet$platform,
  declared_platform      = plat$value,
  platform_tier          = plat$tier,
  declared_type_raw      = sheet$sample_type,
  declared_type          = type$value,
  type_tier              = type$tier,
  declared_molecule_raw  = sheet$molecule,
  declared_molecule      = mol$value,
  molecule_tier          = mol$tier,
  n_flags                = lengths(flags),
  flags                  = vapply(flags, function(f) paste(f, collapse = " | "), character(1)),
  stringsAsFactors       = FALSE
)
write.table(out, out_tsv, sep = "\t", quote = FALSE, row.names = FALSE, na = "")

message(sprintf("Screened %d sample(s): %d flag(s) on %d sample(s).",
                n, sum(out$n_flags), sum(out$n_flags > 0)))
