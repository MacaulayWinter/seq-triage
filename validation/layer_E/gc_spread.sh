#!/usr/bin/env bash
#*******************************************************************************
# gc_spread.sh: per-read GC spread of the first 2,000 reads of each FASTQ
#               (calibration for an isolate-vs-community signal, D86)
#
# WHAT IT DOES  For each file: mean read length, mean GC, the SD of per-read
#               GC, the SD expected from chance alone for reads that long
#               (binomial), and the EXCESS spread beyond chance:
#               sqrt(SD^2 - binomial^2). One genome's reads share a narrow GC
#               range; a community's spread wide. Reads under 50 bp skipped.
# NEEDS         bash, gzip, awk only (runs on an HPC login node without R).
# Usage:  bash gc_spread.sh <out.tsv> <reads.fastq.gz> [more files ...]
#*******************************************************************************
set -euo pipefail
out=$1; shift
printf 'file\tn\tmean_len\tmean_gc\tsd_gc\tbinom_sd\texcess_sd\n' > "$out"
for f in "$@"; do
    gzip -dc "$f" 2>/dev/null | head -n 8000 | awk -v name="$(basename "$f")" '
        NR % 4 == 2 && length($0) >= 50 {
            L = length($0); g = gsub(/[GCgc]/, ""); p = g / L
            n++; sl += L; s += p; ss += p * p; b += p * (1 - p) / L }
        END {
            if (n < 200) exit
            m = s / n; v = (ss - n * m * m) / (n - 1); bv = b / n; ex = v - bv
            printf "%s\t%d\t%.0f\t%.3f\t%.4f\t%.4f\t%.4f\n", name, n, sl / n, m, sqrt(v), sqrt(bv), (ex > 0 ? sqrt(ex) : 0) }' >> "$out" || true
done
