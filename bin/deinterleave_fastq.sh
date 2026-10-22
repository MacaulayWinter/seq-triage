#!/usr/bin/env bash
#*******************************************************************************
# deinterleave_fastq.sh: split an interleaved FASTQ into R1 and R2 (D78)
#
# WHAT IT DOES  Records 1, 3, 5, ... go to <sample>_1.fastq.gz and 2, 4, 6, ...
#               to <sample>_2.fastq.gz. Detection already checked that each
#               pair of records shares a read name; this checks every pair
#               does, and stops with a message naming the first that does
#               not, rather than writing mates that are out of step.
# Usage:  deinterleave_fastq.sh <sample> <threads> <interleaved.fastq.gz>
#*******************************************************************************
set -euo pipefail
sample=$1; threads=$2; in=$3
gzip -dc "$in" | paste - - - - - - - - | awk -F'\t' -v r1="${sample}_1.fastq" -v r2="${sample}_2.fastq" '
    { a = $1; b = $5; sub(/[ \t].*$/, "", a); sub(/[ \t].*$/, "", b); sub(/\/[12]$/, "", a); sub(/\/[12]$/, "", b)
      if (a != b) { printf "record pair %d: read names differ (%s, %s), so the file is not interleaved mates\n", NR, a, b > "/dev/stderr"; exit 1 }
      print $1 "\n" $2 "\n" $3 "\n" $4 > r1; print $5 "\n" $6 "\n" $7 "\n" $8 > r2 }'
bgzip -@ "$threads" "${sample}_1.fastq"
bgzip -@ "$threads" "${sample}_2.fastq"
