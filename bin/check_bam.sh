#!/usr/bin/env bash
#*******************************************************************************
# check_bam.sh: is a BAM file intact, from header to end? (CLAUDE.md D55)
#
# WHAT IT DOES  Two checks, cheapest first:
#               1. samtools quickcheck -u: a BAM header and the end-of-file
#                  marker are present (-u: an unaligned BAM, with no
#                  reference sequences in its header, is fine);
#               2. samtools view -c: every block decompresses and every
#                  record parses, so damage part-way through is found too.
#                  quickcheck alone misses it: a test BAM with 64 bytes
#                  flipped in the middle passed quickcheck.
# WRITES        <sample>.check.tsv: file, status (ok | unreadable), problem,
#               the same columns as check_fastq.sh, so the workflow reads
#               both the same way.
# EXIT          Always 0: a damaged file is a result (a FAIL row), not an
#               error that stops the batch (D45).
#
# Usage:  check_bam.sh <sample> <file.bam>
#*******************************************************************************
set -u
sample=$1
bam=$2
out="${sample}.check.tsv"
printf 'file\tstatus\tproblem\n' > "$out"

if ! samtools quickcheck -u "$bam" 2>/dev/null; then
    printf '%s\tunreadable\t%s\n' "$(basename "$bam")" \
        "not a complete BAM file (no valid header or no end-of-file marker: truncated or not BAM format)" >> "$out"
    exit 0
fi
if ! err=$(samtools view -c "$bam" 2>&1 >/dev/null); then
    first=$(printf '%s\n' "$err" | grep -m1 -E 'E::|ERROR|error' | sed 's/\t/ /g')
    printf '%s\tunreadable\t%s\n' "$(basename "$bam")" \
        "BAM is damaged part-way through: ${first:-samtools could not read every record}" >> "$out"
    exit 0
fi
printf '%s\tok\t\n' "$(basename "$bam")" >> "$out"
