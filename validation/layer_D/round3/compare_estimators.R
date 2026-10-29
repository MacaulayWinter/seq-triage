# compare_estimators.R: layer D round 3 (D99), sylph vs the Kraken2/Bracken
# estimate the report notes beside it, on the same reads, with the layer D
# tolerance. Run from an unpacked results folder (results/, mixtures/).

m <- read.delim("results/triage_metrics.tsv", colClasses = "character", na.strings = "")
t <- read.delim("mixtures/truth.tsv", colClasses = "character"); t <- t[t$score_verdict != "limit", ]
r <- do.call(rbind, lapply(seq_len(nrow(t)), function(i) {
  s <- t$sample[i]; tf <- as.numeric(t$true_fraction[i]); tol <- max(0.002, 0.2 * tf)
  v <- m[m$sample == s & m$metric == t$metric[i], ][1, ]; sy <- as.numeric(sub("%", "", v$value)) / 100
  k2 <- if (t$metric[i] == "off_target_fraction") as.numeric(sub(".*Bracken: ([0-9.]+)%.*", "\\1", v$note)) / 100 else sy
  ok <- function(e) if (tf == 0) e < 0.05 else e > 0 && abs(e - tf) <= tol
  data.frame(sample = s, true = 100 * tf, sylph = 100 * sy, kraken_bracken = 100 * k2, sylph_ok = ok(sy), kb_ok = ok(k2))
}))
cat(sprintf("within tolerance (24 scored): sylph %d, Kraken2/Bracken %d\n", sum(r$sylph_ok), sum(r$kb_ok)))
print(r[!(r$sylph_ok & r$kb_ok), ], row.names = FALSE)
