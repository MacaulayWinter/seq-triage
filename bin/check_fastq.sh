#!/usr/bin/env bash
#*******************************************************************************
# check_fastq.sh: is each FASTQ file readable and well formed? (CLAUDE.md D45)
#
# One pass per file checks both layers:
#   gzip   decompression must complete (catches truncated or corrupted files)
#   FASTQ  every record is 4 lines: header starting "@", sequence, "+" line,
#          and a quality string as long as its sequence
# Needed because seqkit exits 0 on both kinds of damage and simply leaves the
# file out of its table (found in validation round 1 and on nf-core's
# test_1_broken.fastq.gz).
#
# Writes <sample>.check.tsv: sample, file, status (ok | unreadable), problem.
# Always exits 0: a bad file is a RESULT for the report (a FAIL row), and
# the other samples carry on.
#
# Usage:  check_fastq.sh <sample_id> <reads_1> [reads_2]
#*******************************************************************************
set -uo pipefail
sample=$1; shift
out="${sample}.check.tsv"
printf 'sample\tfile\tstatus\tproblem\n' > "$out"

for f in "$@"; do
    gz_err=$(mktemp)
    problem=$(gzip -dc "$f" 2> "$gz_err" | awk '
        NR % 4 == 1 && substr($0, 1, 1) != "@" { printf "record %d: header does not start with @", (NR + 3) / 4; bad = 1; exit }
        NR % 4 == 2 { seq_len = length($0) }
        NR % 4 == 3 && substr($0, 1, 1) != "+" { printf "record %d: third line does not start with +", (NR + 1) / 4; bad = 1; exit }
        NR % 4 == 0 && length($0) != seq_len   { printf "record %d: quality has %d characters but the sequence has %d", NR / 4, length($0), seq_len; bad = 1; exit }
        END { if (!bad && NR % 4 != 0) printf "file ends part-way through record %d", int(NR / 4) + 1
              if (!bad && NR == 0)     printf "file contains no reads" }')
    gz_status=${PIPESTATUS[0]}
    # awk stopping early (a bad record) closes the pipe, and gzip then exits
    # on SIGPIPE, which is not a gzip fault. So gzip's own error is reported
    # only when awk read to the end: then it explains a cut-off last record
    # (a truncated download) better than the record count does.
    if [ "$gz_status" -ne 0 ] && { [ -z "$problem" ] || [[ "$problem" == "file ends part-way"* ]]; }; then
        problem="gzip data is truncated or corrupt ($(head -c 200 "$gz_err" | tr '\n' ' ' | sed 's/ *$//'))"
    fi
    rm -f "$gz_err"
    # An empty file is not unreadable: the pipeline already reports it as
    # "No reads in the file" (FAIL), so it passes this check.
    if [ -z "$problem" ] || [ "$problem" = "file contains no reads" ]; then
        printf '%s\t%s\tok\t\n' "$sample" "$(basename "$f")" >> "$out"
    else
        printf '%s\t%s\tunreadable\t%s\n' "$sample" "$(basename "$f")" "$problem" >> "$out"
    fi
done
