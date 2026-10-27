V <- path.expand("~/seqtriage-scratch/validation/")
a <- read.delim(paste0(V, "gc_spread_layerE_all.tsv"), quote = "")
k <- do.call(rbind, lapply(c("", "_r2", "_r3"), function(r) { x <- read.delim(paste0(V, "kmer_layerE", r, ".tsv"), quote = ""); x$round <- c("r1", "r2", "r3")[match(r, c("", "_r2", "_r3"))]; x }))
k$run <- regmatches(k$file, regexpr("(SRR|ERR|DRR)[0-9]+", k$file))
m <- merge(a, k[, c("run", "round", "bases", "frac_repeated", "est_genome_mb")], by = c("run", "round"))
cat("joined:", nrow(m), "\n\n== frac_repeated by TRUE type and read length\n")
for (r in c("short", "long")) for (t in c("isolate", "metagenome", "amplicon")) { v <- m$frac_repeated[m$truth == t & m$reads == r]
  if (length(v)) cat(sprintf("  %-5s %-10s n=%2d  min %.3f  q10 %.3f  median %.3f  q90 %.3f  max %.3f\n", r, t, length(v), min(v), quantile(v, .1), median(v), quantile(v, .9), max(v))) }
options(width = 220)
mg <- m[m$detected_type == "metagenome" & m$reads == "short", ]
cat("\n== called metagenome, short reads (truth isolate = an error), by redundancy\n")
print(mg[order(-mg$frac_repeated), c("round", "run", "scientific_name", "truth", "type_confidence", "frac_repeated", "excess_sd", "bases")], row.names = FALSE)
tr <- mg[mg$round %in% c("r1", "r2"), ]; te <- mg[mg$round == "r3", ]
right_tr <- tr$frac_repeated[tr$truth == "metagenome" & tr$type_confidence == "high"]
wrong_tr <- tr$frac_repeated[tr$truth != "metagenome"]
cat(sprintf("\nTRAIN r1+r2: correct high metagenomes max %.3f (n=%d); wrongly called metagenome min %.3f (n=%d)\n", max(right_tr), length(right_tr), min(wrong_tr), length(wrong_tr)))
thr <- round((max(right_tr) + min(wrong_tr[wrong_tr > max(right_tr)])) / 2, 3)
cat(sprintf("Threshold from TRAIN (midpoint of the gap): metagenome high only if frac_repeated < %.3f\n", thr))
hi <- te$type_confidence == "high"
cat(sprintf("TEST r3: wrong at high %d -> %d; correct high metagenomes %d -> kept %d\n",
  sum(hi & te$truth != "metagenome"), sum(hi & te$truth != "metagenome" & te$frac_repeated < thr),
  sum(hi & te$truth == "metagenome"), sum(hi & te$truth == "metagenome" & te$frac_repeated < thr)))
write.table(m, paste0(V, "kmer_gc_layerE_all.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
