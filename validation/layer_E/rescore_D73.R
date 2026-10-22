# Re-score a layer E result folder with D73: metagenome HIGH -> medium when the
# top genus (Kraken2 rank G, genus Homo excluded) holds >= 50% of genus-level reads.
a <- commandArgs(TRUE); setwd(a[1])
s <- read.delim("scores.tsv", colClasses = "character")
q <- read.delim("results/qc_values.tsv", colClasses = "character", quote = "")
k <- q[q$tool == "Kraken2" & startsWith(q$field, "G|"), ]
k$taxid <- sapply(strsplit(k$field, "|", fixed = TRUE), `[`, 2); k$n <- as.numeric(k$value)
k <- k[k$taxid != "9605" & k$n > 0, ]; k$run <- sub("^.*_", "", k$sample)
gs <- tapply(k$n, k$run, function(n) max(n) / sum(n))
s$top_genus_share <- round(gs[s$run], 3)
hit <- s$detected_type == "metagenome" & s$type_confidence == "high" & !is.na(s$top_genus_share) & s$top_genus_share >= 0.5
s$type_confidence[hit] <- "medium"
s$outcome <- ifelse(s$detected_type == "uncertain", "uncertain",
             ifelse(s$detected_type == s$truth, paste("correct,", s$type_confidence), paste("WRONG,", s$type_confidence)))
s$wrong_high <- s$outcome == "WRONG, high"; s$changed_by_D73 <- hit
write.table(s, a[2], sep = "\t", quote = FALSE, row.names = FALSE, na = "")
cat(sprintf("%s: %d changed by D73 (%s); wrong at high %d; correct at high %d of %d\n", basename(a[1]), sum(hit),
            paste(s$run[hit], sub(".*_", "", s$outcome[hit]), collapse = ", "), sum(s$wrong_high), sum(s$outcome == "correct, high"), nrow(s)))
