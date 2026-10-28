#!/usr/bin/env Rscript
# =============================================================================
# pathogen_list.R  (D96)
# Turns the published list of bacterial human pathogens (Bartlett et al. 2022,
# Data S1, from Zenodo) into a table of NCBI taxids the contamination screen
# can be matched against. Names change between taxonomy versions, so seq-triage
# matches by taxid (D41); this script does the one name -> taxid step, with
# NCBI's own taxonomy dump, and writes down every name it could not resolve.
#
# Usage:  pathogen_list.R <list.xlsx | zenodo.zip> <taxdump.tar.gz> <version> <doi> <out_prefix>
# Writes: <out_prefix>.tsv             one row per (species, taxid it can be
#                                      reported under): the current taxid plus
#                                      every older taxid merged into it
#         <out_prefix>_unresolved.tsv  names with no single bacterial species
#                                      match, and why
#         <out_prefix>_families.tsv    NCBI family of every bacterial species
#                                      taxid (and older taxids merged into it),
#                                      for the close-relative caveat
# Base R only (D42).
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop("Usage: pathogen_list.R <list.xlsx | zenodo.zip> <taxdump.tar.gz> <version> <doi> <out_prefix>")
}
src <- args[1]; taxdump <- args[2]; version <- args[3]; doi <- args[4]; out_prefix <- args[5]

SHEET        <- "Tab 6 Full List"
LIST_FILE    <- "bacteria_human_pathogens.xlsx"   # Data S1, inside the Zenodo archive
BACTERIA     <- "2"
# NCBI name classes that name the same organism. "authority" entries are the
# taxon's own names with their author citation; "includes" and "in-part" name
# other organisms, so they are never used.
SYNONYM_CLASSES <- c("synonym", "equivalent name", "authority")

# WHAT IT DOES  Reads one sheet of an .xlsx file with base R: an .xlsx is a zip
#               of XML (workbook, shared strings, one file per sheet).
# RETURNS       data.frame of character columns, named by the first row
read_xlsx_sheet <- function(path, sheet) {
  d <- tempfile(); dir.create(d)
  utils::unzip(path, exdir = d)
  rd <- function(f) paste(readLines(file.path(d, f), warn = FALSE, encoding = "UTF-8"), collapse = "")
  unesc <- function(x) {
    x <- gsub("&lt;", "<", x, fixed = TRUE); x <- gsub("&gt;", ">", x, fixed = TRUE)
    x <- gsub("&quot;", "\"", x, fixed = TRUE); x <- gsub("&apos;", "'", x, fixed = TRUE)
    gsub("&amp;", "&", x, fixed = TRUE)
  }
  # The sheet's file, through the workbook's relationship ids.
  wb   <- rd("xl/workbook.xml")
  tags <- regmatches(wb, gregexpr("<sheet [^>]*>", wb))[[1]]
  nm   <- unesc(sub('.*name="([^"]*)".*', "\\1", tags))
  rid  <- sub('.*r:id="([^"]*)".*', "\\1", tags)
  if (!sheet %in% nm) stop(sprintf("Sheet '%s' not found in %s (sheets: %s).", sheet, basename(path), paste(nm, collapse = ", ")))
  rels <- rd("xl/_rels/workbook.xml.rels")
  rtag <- regmatches(rels, gregexpr("<Relationship [^>]*>", rels))[[1]]
  target <- sub('.*Target="([^"]*)".*', "\\1", rtag[grepl(sprintf('Id="%s"', rid[nm == sheet]), rtag, fixed = TRUE)])
  target <- sub("^/?(xl/)?", "xl/", target)
  # Shared strings: the text of each <si>, rich-text runs joined.
  ss <- if (file.exists(file.path(d, "xl/sharedStrings.xml"))) {
    s  <- rd("xl/sharedStrings.xml")
    si <- regmatches(s, gregexpr("<si>.*?</si>", s, perl = TRUE))[[1]]
    vapply(si, function(x) {
      t <- regmatches(x, gregexpr("<t[^>]*>.*?</t>|<t[^>]*/>", x, perl = TRUE))[[1]]
      unesc(paste(sub("<t[^>]*>(.*)</t>", "\\1", sub("<t[^>]*/>", "", t)), collapse = ""))
    }, "", USE.NAMES = FALSE)
  } else character(0)
  x     <- rd(target)
  cells <- regmatches(x, gregexpr("<c [^>]*/>|<c [^>]*>.*?</c>", x, perl = TRUE))[[1]]
  ref   <- sub('.*\\br="([A-Z]+)([0-9]+)".*', "\\1 \\2", sub(">.*", ">", cells))
  col   <- sub(" .*", "", ref); row <- as.integer(sub(".* ", "", ref))
  type  <- ifelse(grepl('^<c [^>]*\\bt="([a-zA-Z]+)"', cells), sub('^<c [^>]*\\bt="([a-zA-Z]+)".*', "\\1", cells), "n")
  val   <- ifelse(grepl("<v>", cells), sub(".*<v>(.*)</v>.*", "\\1", cells), "")
  val[type == "inlineStr"] <- sub(".*<t[^>]*>(.*)</t>.*", "\\1", cells[type == "inlineStr"])
  val[type == "s"] <- ss[as.integer(val[type == "s"]) + 1]
  val <- unesc(val)
  colnum <- vapply(strsplit(col, ""), function(ch) sum((match(ch, LETTERS)) * 26^(rev(seq_along(ch)) - 1)), 0)
  m <- matrix(NA_character_, max(row), max(colnum))
  m[cbind(row, colnum)] <- val
  out <- as.data.frame(m[-1, , drop = FALSE], stringsAsFactors = FALSE)
  names(out) <- m[1, ]
  out[rowSums(!is.na(out) & out != "") > 0, , drop = FALSE]
}

# Whitespace of any kind (27 names in v1.0 end in a non-breaking space)
# collapsed, so it cannot stop a name matching.
clean <- function(x) {
  x <- gsub("\xc2\xa0", " ", enc2utf8(x), fixed = TRUE, useBytes = TRUE)   # UTF-8 non-breaking space, in any locale
  trimws(gsub("[[:space:]]+", " ", x, useBytes = TRUE))
}

# The list ------------------------------------------------------------------
xlsx <- src
if (grepl("\\.zip$", src, ignore.case = TRUE)) {
  inside <- utils::unzip(src, list = TRUE)$Name
  hit <- inside[basename(inside) == LIST_FILE]
  if (length(hit) != 1) stop(sprintf("%s not found (once) in %s.", LIST_FILE, basename(src)))
  d <- tempfile(); utils::unzip(src, files = hit, exdir = d, junkpaths = TRUE)
  xlsx <- file.path(d, LIST_FILE)
}
tab <- read_xlsx_sheet(xlsx, SHEET)
need <- c("genus", "species", "status")
if (!all(need %in% names(tab))) stop(sprintf("Sheet '%s' lacks column(s): %s.", SHEET, paste(setdiff(need, names(tab)), collapse = ", ")))
# Genus and epithet are plain Latin, so anything but letters and hyphens is
# noise: v1.0 has "¬†" (a Mac-encoded non-breaking space) after 27 epithets.
latin <- function(x) gsub("[^A-Za-z-]", "", x, useBytes = TRUE)
raw_species <- clean(paste(tab$genus, tab$species))
lst <- data.frame(species = paste(latin(tab$genus), latin(tab$species)), status = clean(tab$status), stringsAsFactors = FALSE)
n_cleaned <- sum(lst$species != raw_species)
if (!all(lst$status %in% c("established", "putative"))) {
  stop(sprintf("Unexpected status value(s): %s.", paste(unique(setdiff(lst$status, c("established", "putative"))), collapse = ", ")))
}
list_md5 <- unname(tools::md5sum(xlsx))

# NCBI taxonomy --------------------------------------------------------------
td <- tempfile(); dir.create(td)
utils::untar(taxdump, files = c("names.dmp", "nodes.dmp", "merged.dmp"), exdir = td)
dmp <- function(f, cols) {
  x <- readLines(file.path(td, f), warn = FALSE, encoding = "UTF-8")
  p <- strsplit(sub("\t\\|$", "", x), "\t|\t", fixed = TRUE)
  as.data.frame(do.call(rbind, lapply(p, `[`, cols)), stringsAsFactors = FALSE)
}
nodes <- dmp("nodes.dmp", 1:3); names(nodes) <- c("taxid", "parent", "rank")
parent <- setNames(nodes$parent, nodes$taxid); rank <- setNames(nodes$rank, nodes$taxid)
nm_lines <- readLines(file.path(td, "names.dmp"), warn = FALSE, encoding = "UTF-8")
nm_lines <- nm_lines[grepl(paste0("\t\\|\t(", paste(c("scientific name", SYNONYM_CLASSES), collapse = "|"), ")\t\\|$"), nm_lines)]
np <- do.call(rbind, strsplit(sub("\t\\|$", "", nm_lines), "\t|\t", fixed = TRUE))
nms <- data.frame(taxid = np[, 1], name = clean(np[, 2]), class = np[, 4], stringsAsFactors = FALSE)
# NCBI keeps most old names only with their author citation ("Mycobacterium
# smegmatis (Trevisan 1889) Lehmann and Neumann 1899"). The citation is cut off
# only when what remains is exactly "Genus species" and the rest starts like a
# citation ("(", a capital or "["); "subsp." and the like never match.
cite <- nms$class != "scientific name"
nms$name[cite] <- sub("^([A-Z][a-z]+ [a-z][a-z-]+) ([([A-Z]|corrig\\.|emend\\.|ex ).*$", "\\1", nms$name[cite])
# Strain names ("Mycoplasma pneumoniae M129", filed under Mycoplasmoides
# pneumoniae) by their first two words, used only as the last route (resolve()).
st_name <- sub("^([A-Z][a-z]+ [a-z][a-z-]+) .*$", "\\1", nms$name)
strains <- unique(data.frame(taxid = nms$taxid, name = st_name)[st_name %in% lst$species & !nms$name %in% lst$species, ])
nms <- unique(nms[nms$name %in% lst$species, ])
merged <- dmp("merged.dmp", 1:2); names(merged) <- c("old", "new")

# Family of every bacterial species (the close-relative caveat): climb all
# nodes towards the root together, keeping the first family met.
sci <- np[np[, 4] == "scientific name", ]
sci_name <- setNames(sci[, 2], sci[, 1])
family <- setNames(rep(NA_character_, nrow(nodes)), nodes$taxid)
bact <- setNames(rep(FALSE, nrow(nodes)), nodes$taxid)
cur <- nodes$taxid
for (step in 1:100) {
  hit <- is.na(family) & rank[cur] == "family"
  family[hit] <- cur[hit]
  bact[cur == BACTERIA] <- TRUE
  nxt <- unname(parent[cur])
  if (all(nxt == cur)) break
  cur <- nxt
}
sp_bact <- nodes$taxid[nodes$rank == "species" & bact & !is.na(family)]
fam <- data.frame(taxid = sp_bact, family_taxid = family[sp_bact], family = unname(sci_name[family[sp_bact]]), stringsAsFactors = FALSE)
old_ids <- merged[merged$new %in% sp_bact, ]
fam <- rbind(fam, data.frame(taxid = old_ids$old, family_taxid = family[old_ids$new], family = unname(sci_name[family[old_ids$new]]),
                             stringsAsFactors = FALSE))

in_bacteria <- function(id) {
  for (i in 1:100) {
    if (id == BACTERIA) return(TRUE)
    p <- parent[[id]]
    if (is.null(p) || p == id) return(FALSE)
    id <- p
  }
  FALSE
}

# WHAT IT DOES  One list name -> one NCBI bacterial species taxid: scientific
#               name first, then synonyms; anything other than exactly one
#               match is unresolved, with the reason.
resolve <- function(name) {
  for (cls in list("scientific name", SYNONYM_CLASSES)) {
    ids <- unique(nms$taxid[nms$name == name & nms$class %in% cls])
    ids <- ids[rank[ids] %in% "species" & vapply(ids, in_bacteria, TRUE)]
    if (length(ids) == 1) return(c(ids, if (identical(cls, "scientific name")) "scientific_name" else "synonym", ""))
    if (length(ids) > 1) return(c(NA, NA, sprintf("several NCBI species: %s", paste(ids, collapse = ","))))
  }
  # Named in NCBI, but not as a species (a serovar such as Salmonella typhi,
  # a variant such as Mycobacterium bovis): species-level counts cannot tell it
  # apart, so it is not matched. Strains alone are different: when NCBI uses
  # the name only for strains, all under one species, that species is the one
  # the name means (Mycoplasma pneumoniae, filed under Mycoplasmoides).
  below <- unique(nms$taxid[nms$name == name])
  below <- below[!rank[below] %in% "species" & vapply(below, in_bacteria, TRUE)]
  if (length(below) && !all(rank[below] %in% "strain")) {
    return(c(NA, NA, sprintf("NCBI names it below or above species (%s), which species-level counts cannot tell apart",
                             paste(sprintf("%s %s", rank[below], below), collapse = ", "))))
  }
  st <- unique(c(below, strains$taxid[strains$name == name]))
  st <- st[rank[st] %in% "strain" & vapply(st, in_bacteria, TRUE)]
  if (length(st)) {
    sp <- unique(vapply(st, species_of, ""))
    if (length(sp) == 1 && !is.na(sp)) return(c(sp, "strain_names", ""))
    return(c(NA, NA, sprintf("NCBI uses it only for strains, under several species or none: %s", paste(sp, collapse = ","))))
  }
  c(NA, NA, "no bacterial species of this name or synonym in NCBI taxonomy")
}

# The species a strain sits under (NA if none within 100 steps).
species_of <- function(id) {
  for (i in 1:100) {
    if (rank[[id]] == "species") return(id)
    p <- parent[[id]]
    if (is.null(p) || p == id) return(NA_character_)
    id <- p
  }
  NA_character_
}
r <- t(vapply(lst$species, resolve, character(3)))
lst$taxid <- r[, 1]; lst$resolved_by <- r[, 2]; lst$reason <- r[, 3]

ok <- lst[!is.na(lst$taxid), ]
rows <- do.call(rbind, lapply(seq_len(nrow(ok)), function(i) {
  ids <- c(ok$taxid[i], merged$old[merged$new == ok$taxid[i]])
  data.frame(match_taxid = ids, taxid = ok$taxid[i], species = ok$species[i], status = ok$status[i],
             resolved_by = ok$resolved_by[i], stringsAsFactors = FALSE)
}))
# Two list names on one taxid (e.g. Rhodococcus equi and its new name
# Prescottia equi, both listed) become one row naming both; status is
# established if either is.
rows <- do.call(rbind, lapply(split(rows, rows$match_taxid), function(g) {
  data.frame(match_taxid = g$match_taxid[1], taxid = g$taxid[1], species = paste(unique(g$species), collapse = " / "),
             status = if (any(g$status == "established")) "established" else "putative",
             resolved_by = paste(unique(g$resolved_by), collapse = " / "), stringsAsFactors = FALSE)
}))
rows$list_version <- version; rows$list_doi <- doi; rows$list_md5 <- list_md5
rows$taxdump_date <- format(file.mtime(taxdump), "%Y-%m-%d")
rows$list_resolved <- nrow(ok); rows$list_unresolved <- sum(is.na(lst$taxid))
write.table(rows, paste0(out_prefix, ".tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(fam, paste0(out_prefix, "_families.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(lst[is.na(lst$taxid), c("species", "status", "reason")], paste0(out_prefix, "_unresolved.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("%d of %d list names resolved (%d by synonym, %d by strain names; %d names cleaned of stray characters); %d unresolved; %d taxids to match\n",
            nrow(ok), nrow(lst), sum(ok$resolved_by == "synonym"), sum(ok$resolved_by == "strain_names"), n_cleaned, sum(is.na(lst$taxid)), nrow(rows)))
