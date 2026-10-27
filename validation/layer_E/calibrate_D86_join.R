rd <- function(f) read.delim(f, colClasses = "character", quote = "")
V <- "~/seqtriage-scratch/validation/"; O <- "validation/outcomes/"
sets <- list(r1 = c("gc_spread_layerE.tsv", "layer_E/scores_D73.tsv"),
             r2 = c("gc_spread_layerE_r2.tsv", "layer_E_round2/scores_D73.tsv"),
             r3 = c("gc_spread_layerE_r3.tsv", "layer_E_round3/scores.tsv"))
all <- do.call(rbind, lapply(names(sets), function(k) {
  g <- rd(path.expand(paste0(V, sets[[k]][1]))); s <- rd(paste0(O, sets[[k]][2]))
  g$run <- regmatches(g$file, regexpr("(SRR|ERR|DRR)[0-9]+", g$file))
  m <- merge(g, s[, c("run", "stratum", "truth", "detected_type", "type_confidence", "outcome", "scientific_name")], by = "run")
  m$round <- k; m }))
for (c in c("excess_sd", "mean_len")) all[[c]] <- as.numeric(all[[c]])
all$reads <- ifelse(all$mean_len >= 1000, "long", "short")
cat("runs joined:", nrow(all), "\n")
cat("\n== short reads: excess_sd by TRUE type\n")
for (t in c("isolate", "metagenome", "amplicon")) { v <- all$excess_sd[all$truth == t & all$reads == "short"]
  if (length(v)) cat(sprintf("  %-10s n=%2d  min %.3f  q10 %.3f  median %.3f  q90 %.3f  max %.3f\n", t, length(v), min(v), quantile(v, .1), median(v), quantile(v, .9), max(v))) }
cat("\n== long reads: excess_sd by TRUE type\n")
for (t in c("isolate", "metagenome", "amplicon")) { v <- all$excess_sd[all$truth == t & all$reads == "long"]
  if (length(v)) cat(sprintf("  %-10s n=%2d  min %.3f  median %.3f  max %.3f\n", t, length(v), min(v), median(v), max(v))) }
options(width = 220)
cat("\n== called METAGENOME (any confidence), short reads, by truth; sorted by spread\n")
x <- all[all$detected_type == "metagenome" & all$reads == "short", c("round", "run", "stratum", "scientific_name", "truth", "type_confidence", "outcome", "excess_sd")]
print(x[order(x$excess_sd), ], row.names = FALSE)
write.table(all, path.expand(paste0(V, "gc_spread_layerE_all.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
