#!/usr/bin/env Rscript
#*******************************************************************************
# study_metadata.R: seq-triage study metadata and its column guide (CLAUDE.md D95)
#
# WHAT IT DOES  Collects the study's own columns (every samplesheet column the
#               pipeline does not use itself, and optionally a separate
#               metadata table) into one row per sample, and checks them
#               against the column guide: a file the user writes, one row per
#               column, saying what the column does (role) and how to handle
#               it (type). Writes what it found in plain language, the design
#               as read, and a draft guide to fill in.
# THE RULE      Nothing about a column is guessed. Role and type come only from
#               the guide; the data are used to CHECK the guide, never to fill
#               it in. Anything that looks like a typo stops the run, naming
#               it, until the user confirms or fixes it. A column the guide
#               does not cover is never analysed downstream.
# STOPS         (every case listed in one message; the outputs, including the
#               draft guide, are written first):
#               a sample with two different values of one variable; a
#               --metadata table with no sample ID column, or with columns
#               named like pipeline columns; a guide with a missing, unknown or
#               duplicated entry; a value that contradicts its type or its
#               listed levels; values that look like typos of each other
#               (differing only in case, spaces or punctuation, or by one
#               keystroke) in a column whose levels are not listed; in a
#               sample_control column (D97), an empty cell, a value not in
#               CONTROLS, or a second such column.
# WRITES        study_metadata.tsv, metadata_check.tsv, metadata_check.md,
#               column_guide_draft.csv
#
# Usage:  study_metadata.R <samplesheet> <schema_input.json> [--metadata <file>] [--column_guide <file>]
#*******************************************************************************

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: study_metadata.R <samplesheet> <schema_input.json> [--metadata <file>] [--column_guide <file>]")
opt <- function(flag) { i <- match(flag, args); if (is.na(i) || i == length(args)) NA_character_ else args[i + 1] }
sheet_file <- args[1]; schema_file <- args[2]
meta_file <- opt("--metadata"); guide_file <- opt("--column_guide")

# The two closed lists (README, Study metadata). Each ends in a catch-all, so
# every column fits one: a column the pipeline cannot model is `text` or
# `ignore`, described but never analysed.
ROLES <- c(sample_id  = "the sample ID (only needed when a --metadata table does not call it `sample`)",
           subject_id = "the unit sampled more than once (animal, patient, reactor)",
           predictor  = "a variable that may explain differences (treatment, temperature, diet)",
           outcome    = "a measurement on the sample to be explained (MIC, weight)",
           time       = "when the sample was taken (day, visit, date)",
           block      = "a grouping within which comparisons are made (site, litter)",
           batch      = "a technical grouping (extraction run, sequencing lane)",
           sample_control = "which samples are controls (D97): `sample`, or the kind of control",
           ignore     = "kept in the outputs, never analysed")
TYPES <- c(binary      = "two values (yes/no, 0/1)",
           categorical = "named groups with no order",
           ordinal     = "groups with an order; `levels` lists them in order",
           count       = "whole numbers 0 or more",
           continuous  = "any number",
           proportion  = "a number from 0 to 1",
           date        = "a date written YYYY-MM-DD",
           text        = "free text or an identifier: described, never modelled")
# The values of a sample_control column (D97), a closed list. `sample` is a
# real sample; `none` cannot be used, as it is read as an empty cell (EMPTY).
CONTROLS <- c(sample              = "a real sample, not a control",
              negative_extraction = "an extraction blank",
              negative_pcr        = "a PCR or library blank",
              negative_other      = "any other negative control",
              positive_mock       = "a mock community of known make-up",
              positive_other      = "any other positive control")
GUIDE_COLS <- c("column", "role", "type", "levels", "missing_codes", "unit", "nested_in", "derived_from", "description", "seen_in_data")
# Values conventionally used to mean "not recorded" in survey and lab data. A
# numeric column holding one must say whether it is a code (missing_codes) or
# a real value (missing_codes: none); it is never assumed either way.
CONVENTIONAL_CODES <- c("99", "999", "9999", "-9", "-99", "-999", "-9999")
EMPTY <- c("", "na", "n/a", "nan", "null", "none", "missing", ".", "-", "--")
NUM <- "^[-+]?([0-9]+\\.?[0-9]*|\\.[0-9]+)([eE][-+]?[0-9]+)?$"
`%||%` <- function(a, b) if (length(a) && !is.na(a[1])) a[1] else b
q <- function(x) paste(sprintf("`%s`", x), collapse = ", ")

read_table <- function(f) {
  if (!length(readLines(f, n = 1, warn = FALSE)) || !nzchar(trimws(readLines(f, n = 1, warn = FALSE))))
    stop(sprintf("%s is empty: it needs a first row of column names, then one row per sample (or per column, for a column guide).", basename(f)), call. = FALSE)
  sep <- if (grepl("\\.(tsv|txt)$", tolower(f))) "\t" else ","
  x <- read.table(f, sep = sep, header = TRUE, colClasses = "character", na.strings = character(0), quote = "\"",
                  comment.char = "", check.names = FALSE, strip.white = TRUE, fill = TRUE, fileEncoding = "UTF-8-BOM")
  # CRLF (Windows) files leave \r on the last column; removed, never kept.
  names(x) <- trimws(gsub("\r", "", names(x)))
  x[] <- lapply(x, function(v) { v <- trimws(gsub("\r", "", v)); ifelse(tolower(v) %in% EMPTY, NA_character_, v) })
  x
}
# Property names of the samplesheet schema: the columns the pipeline reads.
known_columns <- function(f) {
  s <- readLines(f, warn = FALSE)
  i <- grep('"properties"\\s*:', s)[1]
  ind <- nchar(sub("\\S.*$", "", s[i])) + 4
  hits <- regmatches(s, regexec(sprintf('^ {%d}"([^"]+)"\\s*:', ind), s))
  unique(vapply(Filter(length, hits), `[`, "", 2))
}
stops <- character(0); stop_if <- function(cond, msg) if (isTRUE(cond)) stops <<- c(stops, msg)
findings <- list(); add <- function(column, status, message) findings[[length(findings) + 1]] <<- data.frame(column = column, status = status, message = message)

# What the data show for certain, per column (facts, not choices) ---------------------
is_num <- function(v) grepl(NUM, v)
seen <- function(v) {
  v <- na.omit(v); u <- unique(v)
  if (!length(v)) return("no values")
  if (all(is_num(v))) {
    x <- as.numeric(v); whole <- all(x == round(x))
    # Values listed whenever they repeat, so an odd one (a 99) is visible.
    return(if (length(u) < length(v)) sprintf("%s: %s", if (whole) "whole numbers" else "numbers", paste(format(sort(unique(x)), trim = TRUE), collapse = "|"))
           else sprintf("%s, %s to %s, all different", if (whole) "whole numbers" else "numbers", format(min(x)), format(max(x))))
  }
  if (all(!is.na(as.Date(v, optional = TRUE, tryFormats = c("%Y-%m-%d"))))) return(sprintf("dates, %s to %s", min(v), max(v)))
  if (length(u) == length(v)) return(sprintf("a different value for every sample (%d)", length(u)))
  sprintf("%d values: %s", length(u), paste(sort(u), collapse = "|"))
}
# Column guide ---------------------------------------------------------------------
guide <- NULL
if (!is.na(guide_file)) {
  guide <- read_table(guide_file)
  names(guide) <- tolower(names(guide))
  miss <- setdiff(c("column", "role", "type"), names(guide))
  if (length(miss)) stop(sprintf("--column_guide %s needs the columns column, role and type; missing: %s. Start from assets/column_guide_template.csv.", basename(guide_file), q(miss)))
  odd <- setdiff(names(guide), GUIDE_COLS)
  if (length(odd)) stop(sprintf("--column_guide %s has column(s) the guide does not have: %s (allowed: %s). Nothing has been assumed: check the spelling.", basename(guide_file), q(odd), q(GUIDE_COLS)))
  for (c in setdiff(GUIDE_COLS, names(guide))) guide[[c]] <- NA_character_
  guide$role <- tolower(guide$role); guide$type <- tolower(guide$type)
  guide <- guide[!is.na(guide$column), , drop = FALSE]
  dup <- unique(guide$column[duplicated(guide$column)])
  if (length(dup)) stop(sprintf("--column_guide lists column(s) more than once: %s. Give each column one row.", q(dup)))
}
g <- function(c, field) if (is.null(guide)) NA_character_ else guide[[field]][match(c, guide$column)]

# Read the metadata ----------------------------------------------------------------------
sheet <- read_table(sheet_file)
known <- known_columns(schema_file)
study <- setdiff(names(sheet), known)
sheet$sample <- trimws(sheet$sample)
meta <- sheet[, c("sample", study), drop = FALSE]
source_of <- setNames(rep("samplesheet", length(study)), study)
id_col_meta <- NA_character_

if (!is.na(meta_file)) {
  extra <- read_table(meta_file)
  # The sample ID: a column called `sample`, or the one the guide gives role
  # sample_id. A similar name is never assumed to be it.
  id_col_meta <- if ("sample" %in% names(extra)) "sample" else if (!is.null(guide)) guide$column[guide$role %in% "sample_id"][1] else NA
  if (is.na(id_col_meta) || !id_col_meta %in% names(extra)) {
    like <- names(extra)[grepl("sample|^id$|_id$", tolower(names(extra)))]
    # A draft of every column first, so the user can fill in the whole guide at once.
    write.csv(data.frame(column = names(extra), role = "", type = "", levels = "", missing_codes = "", unit = "", nested_in = "", derived_from = "", description = "",
                         seen_in_data = vapply(names(extra), function(c) seen(extra[[c]]), ""), check.names = FALSE),
              "column_guide_draft.csv", row.names = FALSE, na = "")
    stop(sprintf("--metadata %s has no `sample` column, so its rows cannot be matched to the samplesheet.%s Rename it to `sample`, or give it role sample_id in --column_guide. A draft guide of all its columns is in %s/column_guide_draft.csv.",
                 basename(meta_file), if (length(like)) sprintf(" Column(s) that may be the sample ID: %s.", q(like)) else "", getwd()))
  }
  names(extra)[names(extra) == id_col_meta] <- "sample"
  clash <- intersect(setdiff(names(extra), "sample"), known)
  if (length(clash))
    stop(sprintf("--metadata %s has column(s) named like pipeline columns: %s. They would be mistaken for the pipeline's own declarations. Nothing has been assumed: rename them (for example `sample_type` holding tissue names -> `tissue`); declarations belong in the samplesheet.",
                 basename(meta_file), paste(sprintf("`%s` (e.g. '%s')", clash, vapply(clash, function(c) na.omit(extra[[c]])[1] %||% "", "")), collapse = ", ")))
  new_cols <- setdiff(names(extra), "sample")
  meta <- merge(meta, extra[, c("sample", new_cols), drop = FALSE], by = "sample", all.x = TRUE, all.y = TRUE, suffixes = c("", ".__meta"))
  for (c in intersect(new_cols, study)) {
    a <- meta[[c]]; b <- meta[[paste0(c, ".__meta")]]
    bad <- !is.na(a) & !is.na(b) & a != b
    if (any(bad)) stop(sprintf("Column `%s` differs between the samplesheet and --metadata for sample(s) %s (e.g. '%s' vs '%s'). Make them agree.",
                               c, paste(unique(meta$sample[bad]), collapse = ", "), a[bad][1], b[bad][1]))
    meta[[c]] <- ifelse(is.na(a), b, a); meta[[paste0(c, ".__meta")]] <- NULL
  }
  source_of[setdiff(new_cols, study)] <- "metadata"
  study <- union(study, new_cols)
}

# One row per sample; the rows of one sample (lanes, D75) must agree -----------------
in_sheet <- unique(sheet$sample)
conflicts <- character(0)
per <- lapply(split(meta, meta$sample), function(d) {
  for (c in study) {
    u <- unique(na.omit(d[[c]]))
    if (length(u) > 1) conflicts <<- c(conflicts, sprintf("sample %s, `%s`: %s", d$sample[1], c, paste(sprintf("'%s'", u), collapse = " vs ")))
  }
  data.frame(c(list(sample = d$sample[1]), lapply(d[study], function(v) unique(na.omit(v))[1] %||% NA_character_)),
             check.names = FALSE, stringsAsFactors = FALSE)
})
if (length(conflicts)) stop("Study metadata disagrees within a sample (one value per variable per sample): ", paste(conflicts, collapse = "; "))
md <- if (length(per)) do.call(rbind, per) else data.frame(sample = character(0))
md <- md[order(match(md$sample, in_sheet)), , drop = FALSE]
not_in_sheet <- setdiff(md$sample, in_sheet)
md <- md[md$sample %in% in_sheet, , drop = FALSE]
rownames(md) <- NULL
# Declared missing-value codes become missing, counted, before any check.
code_n <- list()
if (!is.null(guide)) for (c in intersect(guide$column, names(md))) {
  mc <- g(c, "missing_codes"); if (is.na(mc) || tolower(mc) == "none") next
  codes <- trimws(strsplit(mc, "|", fixed = TRUE)[[1]])
  hit <- md[[c]] %in% codes; code_n[[c]] <- sum(hit); md[[c]][hit] <- NA_character_
}
write.table(md, "study_metadata.tsv", sep = "\t", quote = FALSE, row.names = FALSE, na = "")
analysed <- setdiff(study, id_col_meta)   # the ID column is the join key, not a variable

info <- data.frame(column = analysed, source = unname(source_of[analysed]), stringsAsFactors = FALSE)
info$seen <- vapply(analysed, function(c) seen(md[[c]]), "")
info$n_missing <- vapply(analysed, function(c) sum(is.na(md[[c]])), 0L)
info$n_values <- vapply(analysed, function(c) length(unique(na.omit(md[[c]]))), 0L)
info$role <- vapply(analysed, function(c) g(c, "role"), ""); info$type <- vapply(analysed, function(c) g(c, "type"), "")

# Typos: values that look like misspellings of each other -----------------------------
typo_pairs <- function(v) {
  u <- unique(na.omit(v)); u <- u[!is_num(u)]
  if (length(u) < 2) return(character(0))
  out <- character(0)
  key <- tolower(gsub("[[:space:][:punct:]]", "", u))
  for (k in unique(key[duplicated(key)])) out <- c(out, sprintf("%s (differ only in case, spaces or punctuation)", paste(sprintf("'%s'", u[key == k]), collapse = " / ")))
  # One keystroke apart, where neither value is a single character (any two
  # single characters are one keystroke apart) and the values do not differ
  # only in their digits (A1 / A2 are a numbered series, not a typo).
  d <- adist(u)
  for (i in seq_along(u)) for (j in seq_along(u)) if (i < j && d[i, j] == 1 && nchar(u[i]) > 1 && nchar(u[j]) > 1 &&
      gsub("[0-9]", "", u[i]) != gsub("[0-9]", "", u[j]) && key[i] != key[j])
    out <- c(out, sprintf("'%s' / '%s' (one keystroke apart)", u[i], u[j]))
  out
}

# Check each column ------------------------------------------------------------------------
# D97: the real samples. A control's empty cells (a blank has no diet) are expected, not missing.
ctl_col <- if (is.null(guide)) NA_character_ else intersect(guide$column[guide$role %in% "sample_control"], analysed)[1]
real <- if (is.na(ctl_col)) rep(TRUE, nrow(md)) else md[[ctl_col]] %in% "sample"
for (c in analysed) {
  v <- na.omit(md[[c]]); r <- g(c, "role"); t <- g(c, "type")
  lv <- g(c, "levels"); lv <- if (is.na(lv)) character(0) else trimws(strsplit(lv, "|", fixed = TRUE)[[1]])
  if (g(c, "role") %in% c("ignore", "sample_control")) { }   # never analysed, or checked below (D97): no notes about its spread
  else if (!length(v)) add(c, "note", "No values for any sample in this run.")
  else if (length(unique(v)) == 1) add(c, "note", sprintf("The same value ('%s') for every sample: it cannot explain differences in this run.", v[1]))
  n_miss_real <- sum(is.na(md[[c]]) & real)
  if (!g(c, "role") %in% c("ignore", "sample_control") && n_miss_real > 0 && length(v))
    add(c, "note", sprintf("Missing for %d of %d %s.", n_miss_real, sum(real), if (is.na(ctl_col)) "sample(s)" else "real sample(s) (controls not counted)"))
  if (!is.null(code_n[[c]]) && code_n[[c]] > 0) add(c, "note", sprintf("%d value(s) were missing-value codes (%s) and are treated as missing.", code_n[[c]], g(c, "missing_codes")))
  codes_here <- intersect(unique(v), CONVENTIONAL_CODES)
  if (is.null(guide)) {
    if (length(codes_here) && all(is_num(v)))
      add(c, "note", sprintf("Contains %s, often a code for 'not recorded'. In the column guide, list it in missing_codes, or write none there if it is a real value.", paste(sprintf("'%s'", codes_here), collapse = ", ")))
    stop_if(any(is_num(v)) && any(!is_num(v)), sprintf("`%s` mixes numbers and text (e.g. '%s'): a unit or comment typed in a cell? Fix the values, or declare the column categorical or text in --column_guide.",
                                                      c, v[!is_num(v)][1]))
    for (p in typo_pairs(v)) stop_if(TRUE, sprintf("`%s`: %s. Fix the spelling, or confirm both are intended by listing them in `levels` in --column_guide.", c, p))
    next
  }
  if (!c %in% guide$column) { stop_if(TRUE, sprintf("`%s` is not in --column_guide: give it a role and type (role ignore keeps it without analysing it).", c)); next }
  stop_if(is.na(r) || !r %in% names(ROLES), sprintf("`%s`: role '%s' is not one of: %s.", c, r %||% "", paste(names(ROLES), collapse = ", ")))
  stop_if(is.na(t) || !t %in% names(TYPES), sprintf("`%s`: type '%s' is not one of: %s.", c, t %||% "", paste(names(TYPES), collapse = ", ")))
  if (r %in% "sample_control") {
    # D97: every sample says what it is; nothing outside the closed list.
    stop_if(!t %in% "categorical", sprintf("`%s`: a sample_control column must be categorical (it is %s).", c, t %||% ""))
    n_miss <- info$n_missing[match(c, info$column)]
    stop_if(n_miss > 0, sprintf("`%s`: empty for %d sample(s) (e.g. %s). Nothing has been assumed: write `sample` for a real sample, or the kind of control (%s).",
                                c, n_miss, paste(head(md$sample[is.na(md[[c]])], 3), collapse = ", "), paste(names(CONTROLS)[-1], collapse = ", ")))
    out <- setdiff(unique(v), names(CONTROLS))
    near <- vapply(out, function(x) { d <- adist(tolower(x), names(CONTROLS)); if (min(d) <= 3) sprintf(" (did you mean '%s'?)", names(CONTROLS)[which.min(d)]) else "" }, "")
    stop_if(length(out) > 0, sprintf("`%s`: value(s) not among the sample_control values: %s. Use only: %s.", c,
                                     paste(sprintf("'%s'%s", out, near), collapse = ", "), paste(names(CONTROLS), collapse = ", ")))
    stop_if(length(setdiff(lv, names(CONTROLS))) > 0, sprintf("`%s`: `levels` names value(s) that are not sample_control values: %s.", c, paste(sprintf("'%s'", setdiff(lv, names(CONTROLS))), collapse = ", ")))
    next
  }
  if (!t %in% names(TYPES) || !length(v)) next
  stop_if(length(codes_here) > 0 && t %in% c("count", "continuous", "proportion", "ordinal") && is.na(g(c, "missing_codes")),
          sprintf("`%s`: contains %s, often a code for 'not recorded'. Nothing has been assumed: list it in this column's missing_codes, or write none there if it is a real value.",
                  c, paste(sprintf("'%s'", codes_here), collapse = ", ")))
  bad <- switch(t,
    binary      = if (length(unique(v)) > 2) sprintf("declared binary, but has %d values: %s", length(unique(v)), paste(sprintf("'%s'", unique(v)), collapse = ", ")),
    count       = { x <- suppressWarnings(as.numeric(v)); if (any(!is_num(v)) || any(x < 0) || any(x != round(x))) sprintf("declared count, but '%s' is not a whole number 0 or more", v[!is_num(v) | is.na(x) | x < 0 | x != round(x)][1]) },
    continuous  = if (any(!is_num(v))) sprintf("declared continuous, but '%s' is not a number%s", v[!is_num(v)][1],
                                               if (grepl("^[-+]?[0-9]+,[0-9]+$", v[!is_num(v)][1])) " (a decimal comma? write 12.5, not 12,5)" else " (units go in the guide's `unit`, not in the cells)"),
    proportion  = { x <- suppressWarnings(as.numeric(v)); if (any(!is_num(v)) || any(x < 0 | x > 1)) sprintf("declared proportion, but '%s' is not a number from 0 to 1", v[!is_num(v) | is.na(x) | x < 0 | x > 1][1]) },
    date        = if (any(is.na(as.Date(v, optional = TRUE, tryFormats = "%Y-%m-%d")))) sprintf("declared date, but '%s' is not written YYYY-MM-DD", v[is.na(as.Date(v, optional = TRUE, tryFormats = "%Y-%m-%d"))][1]),
    NULL)
  stop_if(!is.null(bad), sprintf("`%s`: %s.", c, bad))
  stop_if(t == "ordinal" && !length(lv), sprintf("`%s`: declared ordinal, so `levels` must list the values in order, separated by |.", c))
  if (length(lv)) {
    out <- setdiff(unique(v), lv)
    stop_if(length(out) > 0, sprintf("`%s`: value(s) not among its listed levels: %s (check spelling and case).", c, paste(sprintf("'%s'", out), collapse = ", ")))
    unused <- setdiff(lv, unique(v))
    if (length(unused)) add(c, "note", sprintf("Listed level(s) not seen in this run: %s.", paste(unused, collapse = ", ")))
  } else if (t %in% c("binary", "categorical") && !r %in% "ignore" && is.na(g(c, "derived_from"))) {
    # Not for a derived column (its values are checked against its sources
    # below, which is stronger) or one never analysed (role ignore).
    for (p in typo_pairs(v)) stop_if(TRUE, sprintf("`%s`: %s. Fix the spelling, or confirm both are intended by listing the column's values in `levels`.", c, p))
  }
  stop_if(r %in% "time" && !t %in% c("ordinal", "count", "continuous", "date"), sprintf("`%s`: a time column must be ordinal, count, continuous or date, so samples can be put in order (it is %s).", c, t))
  stop_if(r %in% "outcome" && t %in% "text", sprintf("`%s`: an outcome cannot be text; give its type, or role ignore.", c))
  for (rel in c("nested_in", "derived_from")) {
    p <- g(c, rel); if (is.na(p)) next
    parents <- trimws(strsplit(p, "|", fixed = TRUE)[[1]])
    miss <- setdiff(parents, analysed)
    stop_if(length(miss) > 0, sprintf("`%s`: %s names column(s) not in the data: %s.", c, rel, q(miss)))
    if (length(miss)) next
    key <- do.call(paste, c(lapply(parents, function(x) md[[x]]), sep = "\r"))
    ok <- !is.na(md[[c]]) & !is.na(key)
    if (rel == "derived_from") {
      n <- tapply(md[[c]][ok], key[ok], function(z) length(unique(z)))
      stop_if(any(n > 1), sprintf("`%s`: declared derived from %s, but one combination of those has %d different values of it.", c, q(parents), max(n)))
    } else {
      n <- tapply(key[ok], md[[c]][ok], function(z) length(unique(z)))
      stop_if(any(n > 1), sprintf("`%s`: declared nested in %s, but value '%s' occurs under %d different values of it.", c, q(parents), names(n)[n > 1][1], max(n)))
    }
  }
}
if (!is.null(guide)) {
  extra_rows <- setdiff(guide$column, c(analysed, id_col_meta, "sample"))
  for (c in extra_rows) stop_if(TRUE, sprintf("`%s` is in --column_guide but not in the data (check the spelling of the column name).", c))
  ctl_cols <- guide$column[guide$role %in% "sample_control"]
  stop_if(length(ctl_cols) > 1, sprintf("--column_guide gives role sample_control to more than one column: %s. Put every sample's status in one column.", q(ctl_cols)))
  stop_if(sum(guide$role %in% "sample_id") > 1, sprintf("--column_guide gives role sample_id to more than one column: %s.", q(guide$column[guide$role %in% "sample_id"])))
}
if (length(not_in_sheet)) add("sample", "note", sprintf("%d sample(s) in --metadata are not in the samplesheet and were left out (e.g. %s).", length(not_in_sheet), paste(head(not_in_sheet, 3), collapse = ", ")))
no_meta <- if (length(analysed)) setdiff(in_sheet, md$sample[rowSums(!is.na(md[, analysed, drop = FALSE])) > 0]) else character(0)
if (length(no_meta)) add("sample", "note", sprintf("%d sample(s) have no study metadata (e.g. %s).", length(no_meta), paste(head(no_meta, 3), collapse = ", ")))

# Columns that always move together: a fact the guide should account for ----------------
cand <- info$column[info$n_values > 1 & !grepl("^a different value", info$seen)]
if (length(cand) >= 2) for (a in seq_along(cand)) for (b in seq_along(cand)) if (a < b) {
  x <- md[[cand[a]]]; y <- md[[cand[b]]]; ok <- !is.na(x) & !is.na(y)
  if (sum(ok) > 1 && all(tapply(y[ok], x[ok], function(z) length(unique(z))) == 1) && all(tapply(x[ok], y[ok], function(z) length(unique(z))) == 1)) {
    declared <- (g(cand[a], "derived_from") %in% cand[b]) || (g(cand[b], "derived_from") %in% cand[a]) ||
                any(c(g(cand[a], "role"), g(cand[b], "role")) %in% "ignore")
    if (!declared)
      add("design", "warning", sprintf("`%s` and `%s` always move together: each value of one comes with one value of the other, so their effects cannot be separated. If one is derived from the other, say so with derived_from (or give it role ignore).", cand[a], cand[b]))
  }
}

# Design as read (only from the guide) ---------------------------------------------------------
design <- character(0)
if (!is.null(guide) && nrow(md) && length(analysed)) {
  role_of <- function(r) info$column[info$role %in% r]
  # D97: controls are counted on their own line, and the design is read from
  # real samples only, so a blank or a mock is never a group or a subject.
  ctl <- role_of("sample_control")[1]
  md_all <- md
  if (!is.na(ctl)) {
    kinds <- table(factor(md[[ctl]], levels = names(CONTROLS)))
    k <- kinds[-1][kinds[-1] > 0]
    design <- c(design, sprintf("- `%s` (sample_control): %d real sample(s); %s", ctl, kinds[["sample"]],
                                if (length(k)) sprintf("controls: %s. Controls are left out of the lines below and never modelled as biology.",
                                                       paste(sprintf("%s (n = %d)", names(k), as.vector(k)), collapse = ", "))
                                else "no controls."))
    md <- md[md[[ctl]] %in% "sample", , drop = FALSE]
    if (!length(k)) add("design", "warning", sprintf("Every sample is a real sample (`%s` is `sample` throughout): this run has no negative or positive controls, so contamination from reagents and the kit cannot be measured.", ctl))
  } else add("design", "warning", "No column has role sample_control, so controls (extraction or PCR blanks, mock communities) cannot be told from real samples: any in this run would be analysed as biology. If there are controls, add a column with role sample_control (README, Study metadata).")
  pred <- role_of("predictor"); subj <- role_of("subject_id"); tim <- role_of("time")
  grouping <- pred[info$type[match(pred, info$column)] %in% c("binary", "categorical", "ordinal")]
  for (c in c(pred, role_of("block"))) {
    t <- info$type[match(c, info$column)]
    design <- c(design, if (t %in% c("count", "continuous", "proportion")) {
      x <- suppressWarnings(as.numeric(md[[c]])); sprintf("- `%s` (%s, %s): %s to %s, median %s; a covariate, not groups.", c, info$role[match(c, info$column)], t,
                                         format(min(x, na.rm = TRUE)), format(max(x, na.rm = TRUE)), format(median(x, na.rm = TRUE)))
    } else sprintf("- `%s` (%s, %s): %s", c, info$role[match(c, info$column)], t,
                   paste(sprintf("%s (n = %d)", names(table(md[[c]])), as.vector(table(md[[c]]))), collapse = ", ")))
  }
  long <- FALSE
  if (length(subj)) {
    n_per <- table(md[[subj[1]]])
    design <- c(design, sprintf("- `%s` (subject_id): %d subjects, %s sample(s) each.", subj[1], length(n_per),
                                if (length(unique(as.vector(n_per))) == 1) n_per[1] else paste(range(n_per), collapse = " to ")))
    if (length(tim)) {
      tt <- tapply(md[[tim[1]]], md[[subj[1]]], function(x) length(unique(na.omit(x))))
      long <- any(tt >= 2)
      design <- c(design, sprintf("- `%s` (time): %d timepoint(s); %d of %d subjects sampled at 2 or more.", tim[1], length(unique(na.omit(md[[tim[1]]]))), sum(tt >= 2), length(tt)))
    }
  } else if (length(tim)) design <- c(design, sprintf("- `%s` (time): %d timepoint(s), but no subject_id column, so repeated samples of one subject cannot be paired.", tim[1], length(unique(na.omit(md[[tim[1]]])))))
  design <- c(design, sprintf("- Layout: **%s**.", if (long) "change over time (the same subjects sampled repeatedly)"
                                                   else if (length(tim) && length(unique(na.omit(md[[tim[1]]]))) > 1) "several timepoints, not paired by subject"
                                                   else "a snapshot (one time)"))
  if (length(grouping)) {
    by <- c(grouping, tim[seq_len(min(1, length(tim)))])
    cells <- table(do.call(paste, c(lapply(by, function(x) md[[x]]), sep = " / ")))
    design <- c(design, sprintf("- Samples per combination of %s: %s.", q(by),
                                if (length(unique(as.vector(cells))) == 1) sprintf("%d in each of %d", cells[1], length(cells))
                                else sprintf("%d to %d over %d combinations", min(cells), max(cells), length(cells))))
    if (any(cells == 1)) add("design", "note", sprintf("%d combination(s) have a single sample: no variation can be estimated within them.", sum(cells == 1)))
  }
}

if (exists("md_all")) md <- md_all

# Draft guide: every column, what the data show, choices left blank ---------------------------
if (length(analysed) || !is.na(id_col_meta)) {
  cols <- c(if (!is.na(id_col_meta) && id_col_meta != "sample") id_col_meta, analysed)
  draft <- data.frame(column = cols,
                      role = vapply(cols, function(c) if (identical(c, id_col_meta)) "sample_id" else g(c, "role") %||% "", ""),
                      type = vapply(cols, function(c) g(c, "type") %||% "", ""),
                      levels = vapply(cols, function(c) g(c, "levels") %||% "", ""),
                      missing_codes = vapply(cols, function(c) g(c, "missing_codes") %||% "", ""),
                      unit = vapply(cols, function(c) g(c, "unit") %||% "", ""),
                      nested_in = vapply(cols, function(c) g(c, "nested_in") %||% "", ""),
                      derived_from = vapply(cols, function(c) g(c, "derived_from") %||% "", ""),
                      description = vapply(cols, function(c) g(c, "description") %||% "", ""),
                      seen_in_data = vapply(cols, function(c) if (c %in% info$column) info$seen[match(c, info$column)] else "the sample ID", ""),
                      check.names = FALSE)
  write.csv(draft, "column_guide_draft.csv", row.names = FALSE, na = "")
}

# Write ----------------------------------------------------------------------------------------
for (s in stops) add(sub("^`([^`]+)`.*$", "\\1", s), "stop", s)
check <- if (length(findings)) do.call(rbind, findings) else data.frame(column = character(0), status = character(0), message = character(0))
write.table(check, "metadata_check.tsv", sep = "\t", quote = FALSE, row.names = FALSE, na = "")
rest <- check[check$status != "stop", , drop = FALSE]
lines <- c("# Study metadata", "",
  if (!length(analysed)) "No study metadata: the samplesheet has only the columns the pipeline uses. Add your study's variables as extra columns (or `--metadata`) and describe each in a column guide (`--column_guide`, template `assets/column_guide_template.csv`) to see the design here and to analyse them downstream."
  else c(sprintf("%d sample(s), %d study column(s)%s.", nrow(md), length(analysed), if (!is.na(meta_file)) sprintf(" (samplesheet and %s)", basename(meta_file)) else ""),
         if (is.null(guide)) c("", "**No column guide given**, so no role or type is known: nothing here will be analysed downstream, and no design is described. `column_guide_draft.csv` lists every column with what the data show; fill in `role` and `type` (README, Study metadata) and give it with `--column_guide`.")
         else sprintf("Column guide: %s.", basename(guide_file)), "",
         "| Column | Role | Type | Seen in the data | Missing |", "|---|---|---|---|---|",
         sprintf("| `%s` | %s | %s | %s | %d |", info$column, ifelse(is.na(info$role), "-", info$role), ifelse(is.na(info$type), "-", info$type), info$seen, info$n_missing), "",
         if (length(design)) c("## Design as read", "", design, "") else NULL),
  if (length(stops)) c("## Stopped: fix these, then run again", "", sprintf("- %s", stops), "") else NULL,
  "## Findings", "",
  if (nrow(rest)) sprintf("- **%s** `%s`: %s", rest$status, rest$column, rest$message) else "- None.",
  "", "Study metadata is not used for QC; it is carried to `study_metadata.tsv` for the analysis pipelines.")
writeLines(lines, "metadata_check.md")
if (length(stops))
  stop(sprintf("Study metadata: %d problem(s) to fix before the run can continue. Nothing has been assumed.\n  - %s\nThe full check and a draft column guide are in %s (metadata_check.md, column_guide_draft.csv).",
               length(stops), paste(stops, collapse = "\n  - "), getwd()), call. = FALSE)
message(sprintf("Study metadata: %d column(s); %d warning(s), %d note(s).", length(analysed), sum(check$status == "warning"), sum(check$status == "note")))
