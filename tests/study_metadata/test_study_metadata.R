#!/usr/bin/env Rscript
# Test for bin/study_metadata.R (CLAUDE.md D95; controls D97)
#
# One case per rule: what must pass, and what must stop (and that the stop
# names the problem). Then a fuzz test: random, messy metadata sheets must
# either pass or stop with a message, never crash with an R error. Every
# case writes its own small files in a temporary folder.
#
# Run from the repository root:
#   Rscript tests/study_metadata/test_study_metadata.R

script <- normalizePath("bin/study_metadata.R"); schema <- normalizePath("assets/schema_input.json")
tmp <- tempfile("sm_"); dir.create(tmp)
csv <- function(x, f) { write.csv(x, file.path(tmp, f), row.names = FALSE, na = "", quote = TRUE); file.path(tmp, f) }
run <- function(sheet, meta = NULL, guide = NULL) {
  force(sheet); force(meta); force(guide)   # paths are worked out before changing folder
  d <- tempfile("run_", tmp); dir.create(d); old <- setwd(d); on.exit(setwd(old))
  a <- c("--vanilla", script, sheet, schema, if (!is.null(meta)) c("--metadata", meta), if (!is.null(guide)) c("--column_guide", guide))
  out <- suppressWarnings(system2("Rscript", a, stdout = TRUE, stderr = TRUE))
  list(status = attr(out, "status") %||% 0L, text = paste(out, collapse = "\n"), dir = d,
       md = if (file.exists("metadata_check.md")) paste(readLines("metadata_check.md"), collapse = "\n") else "",
       draft = if (file.exists("column_guide_draft.csv")) read.csv("column_guide_draft.csv", colClasses = "character") else NULL,
       meta = if (file.exists("study_metadata.tsv")) read.delim("study_metadata.tsv", colClasses = "character", check.names = FALSE) else NULL)
}
`%||%` <- function(a, b) if (is.null(a)) b else a
options(width = 300)
res <- data.frame(case = character(0), got = character(0), expected = character(0))
expect <- function(case, got, expected) res <<- rbind(res, data.frame(case = case, got = as.character(got), expected = as.character(expected)))
reads <- function(ids) data.frame(sample = ids, fastq_1 = paste0(ids, "_R1.fastq.gz"), fastq_2 = paste0(ids, "_R2.fastq.gz"))
G <- function(...) { x <- rbind(...); colnames(x) <- c("column", "role", "type", "levels", "missing_codes", "derived_from", "nested_in"); as.data.frame(x) }
row <- function(column, role, type, levels = "", missing = "", derived = "", nested = "") c(column, role, type, levels, missing, derived, nested)
ids <- sprintf("S%02d", 1:8)
base <- cbind(reads(ids), trt = rep(c("control", "drug"), 4), temp = rep(c(26, 26, 30, 30), 2), day = rep(c(0, 7), each = 4),
              animal = rep(c("A", "B", "C", "D"), 2), weight = c(20.1, 19.8, 21.0, 20.4, 21.2, 20.0, 21.9, 20.8))
guide_ok <- G(row("trt", "predictor", "binary", "control|drug"), row("temp", "predictor", "categorical", "26|30"),
              row("day", "time", "count"), row("animal", "subject_id", "text"), row("weight", "outcome", "continuous", missing = "none"))

# Passes, and reads the design --------------------------------------------------------
r <- run(csv(base, "s.csv"), guide = csv(guide_ok, "g.csv"))
expect("OK_EXIT", r$status, 0)
expect("OK_LONGITUDINAL", grepl("change over time", r$md), TRUE)
expect("OK_CONTINUOUS_IS_COVARIATE", grepl("`weight`", r$md) && !grepl("Samples per combination of .*weight", r$md), TRUE)
expect("OK_DISCRETE_IS_GROUPS", grepl("`temp` \\(predictor, categorical\\): 26 \\(n = 4\\), 30 \\(n = 4\\)", r$md), TRUE)
# The same temperature declared continuous is a covariate: no threshold decides it.
r <- run(csv(base, "s.csv"), guide = csv(within(guide_ok, type[column == "temp"] <- "continuous"), "g2.csv"))
expect("TEMP_DECLARED_CONTINUOUS", grepl("`temp` \\(predictor, continuous\\): 26 to 30", r$md), TRUE)
# No guide: nothing about roles is guessed, and a draft with blank choices is written.
r <- run(csv(base, "s.csv"))
expect("NOGUIDE_EXIT", r$status, 0)
expect("NOGUIDE_NO_DESIGN", grepl("Design as read", r$md), FALSE)
expect("NOGUIDE_DRAFT_BLANK_CHOICES", !is.null(r$draft) && all(r$draft$role == "" & r$draft$type == ""), TRUE)
expect("NOGUIDE_DRAFT_SHOWS_VALUES", r$draft$seen_in_data[r$draft$column == "trt"], "2 values: control|drug")

# Stops: every one names the problem --------------------------------------------------
stops <- function(case, sheet, pattern, meta = NULL, guide = NULL) {
  r <- run(csv(sheet, paste0(case, "_s.csv")), if (!is.null(meta)) csv(meta, paste0(case, "_m.csv")), if (!is.null(guide)) csv(guide, paste0(case, "_g.csv")))
  expect(case, if (r$status != 0 && grepl(pattern, r$text)) "stopped, named" else sprintf("status %s: %s", r$status, substr(r$text, 1, 160)), "stopped, named")
  invisible(r)
}
stops("STOP_LANES_DISAGREE", rbind(base, transform(base[1, ], fastq_1 = "x_R1.fastq.gz", fastq_2 = "x_R2.fastq.gz", trt = "drug")), "sample S01, `trt`")
stops("STOP_CASE_TYPO", within(base, trt[1] <- "Control"), "'control' / 'Control'|'Control' / 'control'")
stops("STOP_ONE_KEYSTROKE", within(base, trt[1] <- "contol"), "one keystroke apart")
stops("STOP_NUMBERS_AND_TEXT", within(base, weight[2] <- "19.8 g"), "mixes numbers and text")
stops("STOP_NOT_IN_GUIDE", cbind(base, site = "X"), "`site` is not in --column_guide", guide = guide_ok)
stops("STOP_GUIDE_ROW_NOT_IN_DATA", base, "`weigth` is in --column_guide but not in the data", guide = rbind(guide_ok, row("weigth", "outcome", "continuous")))
stops("STOP_UNKNOWN_ROLE", base, "role 'predicter' is not one of", guide = within(guide_ok, role[column == "trt"] <- "predicter"))
stops("STOP_UNKNOWN_TYPE", base, "type 'factor' is not one of", guide = within(guide_ok, type[column == "trt"] <- "factor"))
stops("STOP_OUTSIDE_LEVELS", within(base, trt[3] <- "Drug"), "not among its listed levels: 'Drug'", guide = guide_ok)
stops("STOP_UNIT_IN_CELL", within(base, weight[2] <- "19.8g"), "'19.8g' is not a number", guide = guide_ok)
stops("STOP_DECIMAL_COMMA", within(base, weight[2] <- "19,8"), "decimal comma", guide = guide_ok)
stops("STOP_BINARY_THREE", within(base, trt[3] <- "placebo"), "declared binary, but has 3 values", guide = within(guide_ok, levels[column == "trt"] <- ""))
stops("STOP_COUNT_FRACTION", within(base, day[2] <- "0.5"), "declared count, but '0.5'", guide = guide_ok)
stops("STOP_ORDINAL_NO_LEVELS", base, "declared ordinal, so `levels`", guide = within(guide_ok, type[column == "day"] <- "ordinal"))
stops("STOP_TIME_NOT_ORDERED", base, "a time column must be", guide = within(guide_ok, { type[column == "day"] <- "categorical" }))
stops("STOP_MISSING_CODE_UNDECLARED", within(base, weight[4] <- "99"), "contains '99', often a code", guide = within(guide_ok, missing_codes[column == "weight"] <- ""))
stops("STOP_DERIVED_INCONSISTENT", cbind(base, grp = c("C", "D", "C", "D", "C", "D", "C", "C")), "declared derived from `trt`",
      guide = rbind(guide_ok, row("grp", "ignore", "categorical", derived = "trt")))
stops("STOP_NESTED_INCONSISTENT", cbind(base, cage = c("K1", "K1", "K2", "K2", "K1", "K2", "K2", "K2")), "declared nested in `cage`",
      guide = rbind(guide_ok[guide_ok$column != "animal", ], row("animal", "subject_id", "text", nested = "cage"), row("cage", "block", "categorical")))
stops("STOP_META_NO_SAMPLE", reads(ids), "no `sample` column", meta = data.frame(sample_id = ids, trt = base$trt))
stops("STOP_META_PIPELINE_NAME", reads(ids), "`sample_type` \\(e.g. 'plaque'\\)", meta = data.frame(sample = ids, sample_type = "plaque"))
stops("STOP_GUIDE_TYPO_HEADER", base, "column\\(s\\) the guide does not have: `lvels`", guide = setNames(guide_ok, sub("levels", "lvels", names(guide_ok))))

empty <- file.path(tmp, "empty.csv"); file.create(empty)
r <- run(csv(reads(ids), "r_e.csv"), meta = empty)
expect("STOP_EMPTY_METADATA", r$status != 0 && grepl("empty.csv is empty", r$text), TRUE)
hdr <- file.path(tmp, "header_only.csv"); writeLines("sample,trt", hdr)
r <- run(csv(reads(ids), "r_h.csv"), meta = hdr)
expect("HEADER_ONLY_METADATA_NO_CRASH", r$status == 0 || grepl("Study metadata|--metadata", r$text), TRUE)
r <- run(csv(base, "s_eg.csv"), guide = empty)
expect("STOP_EMPTY_GUIDE", r$status != 0 && grepl("is empty", r$text), TRUE)

# Stops that must NOT happen ------------------------------------------------------------
r <- run(csv(transform(base, animal = rep(c("A1", "A2", "A3", "A4"), 2)), "series.csv"))
expect("NO_STOP_NUMBERED_SERIES", r$status, 0)
r <- run(csv(transform(base, tooth = rep(c("H", "E", "D", "H"), 2)), "single.csv"))
expect("NO_STOP_SINGLE_CHARACTERS", r$status, 0)
r <- run(csv(within(base, trt[1] <- "Control"), "confirmed.csv"), guide = csv(within(guide_ok, { type[column == "trt"] <- "categorical"; levels[column == "trt"] <- "control|Control|drug" }), "gc.csv"))
expect("NO_STOP_CONFIRMED_BY_LEVELS", r$status, 0)
r <- run(csv(within(base, weight[4] <- "99"), "code.csv"), guide = csv(within(guide_ok, missing_codes[column == "weight"] <- "99"), "gcode.csv"))
expect("MISSING_CODE_BECOMES_MISSING", if (!is.null(r$meta)) is.na(r$meta$weight[4]) || r$meta$weight[4] == "" else "no output", TRUE)
r <- run(csv(reads(ids), "r.csv"), meta = csv(data.frame(sid = ids, trt = base$trt), "m_sid.csv"), guide = csv(G(row("sid", "sample_id", "text"), row("trt", "predictor", "binary")), "g_sid.csv"))
expect("META_ID_FROM_GUIDE", r$status, 0)
# CRLF (Windows) files read the same as Unix ones.
crlf <- file.path(tmp, "crlf.csv"); writeLines(gsub("\n", "\r\n", paste(readLines(csv(base, "u.csv")), collapse = "\n")), crlf, sep = "\r\n")
r <- run(crlf, guide = csv(guide_ok, "g_crlf.csv"))
expect("CRLF_READ", r$status == 0 && !any(grepl("\r", unlist(r$meta))), TRUE)
# Confounding is reported (a fact about the design), not a stop.
r <- run(csv(cbind(base, batch = rep(c("B1", "B2"), 4)), "conf.csv"), guide = csv(rbind(guide_ok, row("batch", "batch", "categorical")), "g_conf.csv"))
expect("CONFOUNDED_WARNS_NOT_STOPS", r$status == 0 && grepl("`trt` and `batch` always move together", r$md), TRUE)
# The shipped example passes cleanly.
r <- run(normalizePath("assets/examples/samplesheet_example.csv"), meta = normalizePath("assets/examples/metadata_example.csv"),
         guide = normalizePath("assets/examples/column_guide_example.csv"))
expect("EXAMPLE_PASSES_NO_WARNING", r$status == 0 && !grepl("\\*\\*warning\\*\\*", r$md), TRUE)
expect("TEMPLATE_HEADER", readLines("assets/column_guide_template.csv")[1], "column,role,type,levels,missing_codes,unit,nested_in,derived_from,description")

# Controls: the sample_control role (D97) -----------------------------------------------
ctl <- rbind(cbind(base, kind = "sample"),
             data.frame(reads(c("NEG1", "MOCK1")), trt = "", temp = "", day = "", animal = "", weight = "", kind = c("negative_extraction", "positive_mock")))
guide_ctl <- rbind(guide_ok, row("kind", "sample_control", "categorical"))
r <- run(csv(ctl, "ctl.csv"), guide = csv(guide_ctl, "g_ctl.csv"))
expect("CTL_OK_EXIT", r$status, 0)
expect("CTL_COUNTED_BY_KIND", grepl("`kind` \\(sample_control\\): 8 real sample\\(s\\); controls: negative_extraction \\(n = 1\\), positive_mock \\(n = 1\\)", r$md), TRUE)
expect("CTL_NOT_IN_GROUPS", grepl("`trt` \\(predictor, binary\\): control \\(n = 4\\), drug \\(n = 4\\)\n", r$md) &&
                            grepl("Samples per combination of `trt`, `temp`, `day`: 1 in each of 8\\.", r$md), TRUE)
expect("CTL_NOT_SUBJECTS", grepl("`animal` \\(subject_id\\): 4 subjects, 2 sample\\(s\\) each", r$md), TRUE)
expect("CTL_BLANKS_NOT_MISSING", grepl("Missing for", r$md), FALSE)
expect("CTL_NO_WARNING", grepl("\\*\\*warning\\*\\*", r$md), FALSE)
expect("CTL_KEPT_IN_OUTPUT", if (!is.null(r$meta)) paste(r$meta$kind[r$meta$sample %in% c("NEG1", "MOCK1")], collapse = ",") else "no output", "negative_extraction,positive_mock")
stops("STOP_CTL_CASE", within(ctl, kind[9] <- "negative_PCR"), "'negative_PCR' \\(did you mean 'negative_pcr'\\?\\)", guide = guide_ctl)
stops("STOP_CTL_UNKNOWN", within(ctl, kind[9] <- "blank"), "not among the sample_control values: 'blank'", guide = guide_ctl)
stops("STOP_CTL_EMPTY", within(ctl, kind[3] <- ""), "`kind`: empty for 1 sample\\(s\\) \\(e.g. S03\\)", guide = guide_ctl)
stops("STOP_CTL_NONE_IS_EMPTY", within(ctl, kind[3] <- "none"), "write `sample` for a real sample", guide = guide_ctl)
stops("STOP_CTL_TWO_COLUMNS", cbind(ctl, kind2 = ctl$kind), "role sample_control to more than one column", guide = rbind(guide_ctl, row("kind2", "sample_control", "categorical")))
stops("STOP_CTL_TYPE", ctl, "a sample_control column must be categorical \\(it is text\\)", guide = within(guide_ctl, type[column == "kind"] <- "text"))
stops("STOP_CTL_LEVELS", ctl, "`levels` names value\\(s\\) that are not sample_control values: 'blank'", guide = within(guide_ctl, levels[column == "kind"] <- "sample|blank"))
r <- run(csv(base, "s_noctl.csv"), guide = csv(guide_ok, "g_noctl.csv"))
expect("WARN_NO_CTL_COLUMN", r$status == 0 && grepl("\\*\\*warning\\*\\* `design`: No column has role sample_control", r$md), TRUE)
r <- run(csv(cbind(base, kind = "sample"), "s_allreal.csv"), guide = csv(guide_ctl, "g_allreal.csv"))
expect("WARN_ALL_REAL", r$status == 0 && grepl("Every sample is a real sample", r$md), TRUE)

# Fuzz: messy random sheets either pass or stop with a message, never crash --------------
set.seed(95)
crash <- 0L
for (k in 1:300) {
  n <- sample(2:20, 1); id <- sprintf("X%03d", seq_len(n)); s <- reads(id)
  for (j in seq_len(sample(1:12, 1))) {
    nm <- paste0(sample(c("Trt", "day ", "élan", "a b", "x.1", "ID2", "dosé (mg)", "#n", "time_h", "grp"), 1), j)
    s[[nm]] <- switch(sample(5, 1), sample(c("a", "b", "A", " b", "", "NA"), n, TRUE), sample(c(1, 2.5, -9, 99, NA), n, TRUE),
                      sample(c("2026-01-0", "2026-01-02", "x"), n, TRUE), paste0("v", seq_len(n)), rep("same", n))
  }
  r <- run(csv(s, sprintf("f%03d.csv", k)))
  if (r$status != 0 && (grepl("Error in |object '.*' not found|subscript out of bounds|argument is of length zero", r$text) || !grepl("Study metadata|--metadata|column", r$text))) {
    crash <- crash + 1L; if (crash == 1) message("first crash: ", substr(r$text, 1, 400))
  }
}
expect("FUZZ_300_NO_CRASH", crash, 0)

res$result <- ifelse(res$got == res$expected, "PASS", "FAIL")
print(res, row.names = FALSE, right = FALSE)
n_fail <- sum(res$result == "FAIL")
cat(sprintf("\n%d / %d cases as expected\n", nrow(res) - n_fail, nrow(res)))
unlink(tmp, recursive = TRUE)
if (n_fail > 0) quit(status = 1)
