#!/usr/bin/env Rscript
# Test for bin/pathogen_list.R and the known-pathogen section of
# bin/triage_report.R (D96)
#
# 1. Builds the list from SYNTHETIC fixtures (tests/pathogens/): a small .xlsx
#    shaped like Data S1, inside a Zenodo-style archive, and a small NCBI
#    taxonomy dump. Each list name exercises one way of resolving a name, or
#    one reason it cannot be resolved.
# 2. Runs the report on synthetic nanopore samples (seqkit, NanoPlot, Kraken2)
#    with that list, without one, and without the contamination screen, and
#    checks each sample's flags and caveats.
#
# Run from the repository root:
#   Rscript tests/pathogens/test_pathogen_flags.R
# Exits non-zero if any check differs, so CI can use it.

options(width = 300)
tmp <- tempfile("pathogen_test_")
dir.create(tmp)
write_tsv <- function(x, f) write.table(x, f, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
rd <- function(f) read.delim(f, colClasses = "character", na.strings = "", quote = "")
res <- list()
check <- function(case, got, expected) {
  res[[length(res) + 1]] <<- data.frame(case = case, got = got, expected = expected,
                                        result = if (identical(got, expected)) "PASS" else "FAIL")
}

# 1. The list builder ------------------------------------------------------------
plist <- file.path(tmp, "pathogen_list"); dir.create(plist)
status <- system2("Rscript", c("--vanilla", "bin/pathogen_list.R", "tests/pathogens/fixture_record.zip",
                               "tests/pathogens/fixture_taxdump.tar.gz", "v9.9", "10.0000/test",
                               file.path(plist, "pathogen_list")), stdout = FALSE, stderr = FALSE)
if (status != 0) stop("pathogen_list.R exited with status ", status)
pl <- rd(file.path(plist, "pathogen_list.tsv"))
un <- rd(file.path(plist, "pathogen_list_unresolved.tsv"))
fm <- rd(file.path(plist, "pathogen_list_families.tsv"))
taxid_of <- function(name) { t <- unique(pl$taxid[vapply(strsplit(pl$species, " / "), function(x) name %in% x, TRUE)]); if (length(t)) t else "unresolved" }
how_of   <- function(name) { r <- pl[vapply(strsplit(pl$species, " / "), function(x) name %in% x, TRUE) & pl$match_taxid == pl$taxid, ]; if (nrow(r)) r$resolved_by[1] else "unresolved" }
check("LIST_SCIENTIFIC_NAME",       taxid_of("Escherichia coli"), "562")          # a eukaryote of the same name is ignored
check("LIST_SYNONYM_WITH_AUTHOR",   taxid_of("Mycobacterium smegmatis"), "1772")
check("LIST_SYNONYM_CORRIG",        taxid_of("Mycobacterium chelonae"), "1774")
check("LIST_AUTHORITY_NAME",        taxid_of("Lactobacillus rhamnosus"), "47715")  # also cleaned of a stray "¬†"
check("LIST_STRAIN_NAMES_ONLY",     paste(taxid_of("Mycoplasma pneumoniae"), how_of("Mycoplasma pneumoniae")), "2104 strain_names")
check("LIST_TWO_NAMES_ONE_TAXID",   paste(sum(pl$match_taxid == "43767"), pl$status[pl$match_taxid == "43767"]), "1 established")
check("LIST_MERGED_OLD_TAXID",      pl$taxid[pl$match_taxid == "1001"], "562")
check("LIST_SEROVAR_UNRESOLVED",    grepl("below or above species", un$reason[un$species == "Salmonella typhi"]), TRUE)
check("LIST_AMBIGUOUS_UNRESOLVED",  grepl("several NCBI species", un$reason[un$species == "Ochrobactrum anthropi"]), TRUE)
check("LIST_ABSENT_UNRESOLVED",     grepl("no bacterial species", un$reason[un$species == "Madeupia nonexistens"]), TRUE)
check("LIST_COUNTS",                paste(pl$list_resolved[1], pl$list_unresolved[1]), "7 3")
check("LIST_FAMILIES",              fm$family_taxid[match(c("562", "28901", "43767", "1001"), fm$taxid)], c("543", "543", "85025", "543"))

# 2. The report ----------------------------------------------------------------------
SEQKIT <- c(num_seqs = 250000, sum_len = 75000000, sum_n = 1000)
NANO   <- c(number_of_reads = 20000, number_of_bases = 1e8, n50 = 10000, mean_qual = 15, median_read_length = 5000)
# Species-level Kraken2 counts per case ("rank|taxid|name" = reads), and the
# flags expected: "reported unit = species listed [caveat main species]".
cases <- list(
  list(id = "ECOLI_GROUP", k = c("S|562|Escherichia coli" = 900000, "S|623|Shigella flexneri" = 30000),
       expect = "Escherichia coli / Shigella=Escherichia coli"),
  list(id = "OLD_TAXID", k = c("S|1001|Escherichia coli" = 900000),
       expect = "Escherichia coli=Escherichia coli"),
  list(id = "SAME_FAMILY_CAVEAT", k = c("S|28901|Salmonella enterica" = 900000, "S|562|Escherichia coli" = 50000),
       expect = "Escherichia coli / Shigella=Escherichia coli [Salmonella enterica]"),
  list(id = "OTHER_FAMILY_NO_CAVEAT", k = c("S|562|Escherichia coli" = 900000, "S|43767|Prescottella equi" = 50000),
       expect = c("Escherichia coli / Shigella=Escherichia coli", "Prescottella equi=Rhodococcus equi / Prescottella equi")),
  list(id = "UNDER_1PCT_NOT_FLAGGED", k = c("S|562|Escherichia coli" = 900000, "S|43767|Prescottella equi" = 5000),
       expect = "Escherichia coli / Shigella=Escherichia coli"),
  list(id = "NONE_LISTED", k = c("S|28901|Salmonella enterica" = 900000), expect = "none_listed"),
  # The declared target is the main species even when another has more reads.
  list(id = "TARGET_IS_MAIN", target = "28901", k = c("S|562|Escherichia coli" = 600000, "S|28901|Salmonella enterica" = 400000),
       expect = "Escherichia coli / Shigella=Escherichia coli [Salmonella enterica]")
)
build <- function(cases, dir, screen_on = TRUE) {
  sc <- sa <- qc <- list()
  for (cs in cases) {
    sc[[cs$id]] <- data.frame(sample = cs$id, declared_platform_raw = "ont", declared_platform = "ont", declared_type = "isolate",
                              declared_molecule = NA, declared_library = NA, n_flags = "0", flags = NA,
                              declared_target = if (is.null(cs$target)) NA else cs$target, declared_genome_size = NA)
    sa[[cs$id]] <- data.frame(sample = cs$id, single_end = "true", n_files = "1", n_reads = "20000", detected_platform = "ont",
                              detected_layout = "single", confidence = "high", detection_flags = NA, reasons = "synthetic",
                              route = "ont", routed_by = "detected", declared_platform = "ont", declared_from = "samplesheet",
                              start_share = NA, probe_reads = NA, gc_excess = "0.01", kmer_repeated = "0.01",
                              kmer_bases = "7500000", kmer_genome_mb = "5.0")
    v <- function(tool, x) data.frame(sample = cs$id, tool = tool, file = "f", field = names(x), value = format(x, scientific = FALSE, trim = TRUE))
    qc[[paste0(cs$id, "s")]] <- v("seqkit", SEQKIT)
    qc[[paste0(cs$id, "n")]] <- v("NanoPlot", NANO)
    if (screen_on) {
      k <- c("U|0|unclassified" = 10000, "R|1|root" = sum(cs$k), "S|9606|Homo sapiens" = 0, cs$k)
      qc[[paste0(cs$id, "k")]] <- v("Kraken2", k)
    }
  }
  write_tsv(do.call(rbind, sc), file.path(dir, "screen.tsv"))
  write_tsv(do.call(rbind, sa), file.path(dir, "samples.tsv"))
  write_tsv(do.call(rbind, qc), file.path(dir, "qc.tsv"))
}
run_report <- function(dir, list_dir, note) {
  status <- system2("Rscript", c("--vanilla", "bin/triage_report.R", file.path(dir, "screen.tsv"), file.path(dir, "samples.tsv"),
                                 "assets/metric_guide.tsv", file.path(dir, "qc.tsv"), "assets/species_groups.tsv",
                                 file.path(dir, "out"), "assets/species_genome_size.tsv", list_dir, shQuote(note)),
                    stdout = FALSE, stderr = FALSE)
  if (status != 0) stop("triage_report.R exited with status ", status, " in ", dir)
  p <- rd(file.path(dir, "out_pathogens.tsv"))
  p$caveat[is.na(p$caveat)] <- ""
  list(p = p, md = readLines(file.path(dir, "out_summary.md")))
}

a <- file.path(tmp, "with_list"); dir.create(a); build(cases, a)
ra <- run_report(a, plist, "a test copy")
for (cs in cases) {
  r <- ra$p[ra$p$sample == cs$id, ]
  got <- if (all(r$result == "listed"))
    sprintf("%s=%s%s", r$reported_as, r$species, ifelse(nzchar(r$caveat), sprintf(" [%s]", sub("^may be misread DNA of (.*) \\(same family\\).*$", "\\1", r$caveat)), ""))
  else r$result
  check(cs$id, sort(got), sort(cs$expect))
}
check("UNRESOLVED_PUBLISHED", file.exists(file.path(a, "pathogen_list_unresolved.tsv")), TRUE)
check("REPORT_SECTION_AND_CITATION", any(grepl("Bartlett, Padfield, Lear, Bendall & Vos", ra$md)) &&
        any(grepl("version v9.9 \\(doi:10.0000/test\\), a test copy", ra$md)), TRUE)
check("REPORT_SAYS_INFORMATION_ONLY", any(grepl("^\\*\\*Information only.\\*\\*", ra$md)), TRUE)

# Flags change no verdict: the same samples give the same summary without a list.
b <- file.path(tmp, "no_list"); dir.create(b); build(cases, b)
rb <- run_report(b, "assets/no_pathogen_list", "turned off with --skip_pathogen_flags")
check("NO_LIST_NOT_CHECKED", unique(paste(rb$p$result, rb$p$note)), "not_checked turned off with --skip_pathogen_flags")
check("NO_LIST_SAME_VERDICTS", identical(readLines(file.path(a, "out_summary.tsv")), readLines(file.path(b, "out_summary.tsv"))) &&
        identical(readLines(file.path(a, "out_metrics.tsv")), readLines(file.path(b, "out_metrics.tsv"))), TRUE)

c_ <- file.path(tmp, "no_screen"); dir.create(c_); build(cases[1], c_, screen_on = FALSE)
rc <- run_report(c_, plist, "a test copy")
check("NO_SCREEN_NOT_CHECKED", paste(rc$p$result, rc$p$note), "not_checked contamination screen did not run")

res <- do.call(rbind, lapply(res, function(r) { r$got <- paste(r$got, collapse = " | "); r$expected <- paste(r$expected, collapse = " | "); r[1, ] }))
print(res, row.names = FALSE, right = FALSE)
n_fail <- sum(res$result == "FAIL")
cat(sprintf("\n%d / %d cases as expected\n", nrow(res) - n_fail, nrow(res)))
quit(status = if (n_fail) 1 else 0)
