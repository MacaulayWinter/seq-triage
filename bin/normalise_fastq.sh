#!/usr/bin/env bash
#*******************************************************************************
# normalise_fastq.sh: make every read file gzipped FASTQ (CLAUDE.md D69)
#
# WHAT IT DOES  For each file: .gz is passed through untouched (linked); plain
#               FASTQ is compressed; .bz2, .xz and .zst are decompressed and
#               recompressed as gzip (bgzip, which any gzip reader reads).
#               Everything after this step sees .fastq.gz only.
# INTEGRITY     The decompressor's own exit status is checked (pipefail): a
#               truncated or corrupt .bz2/.xz/.zst becomes an "unreadable" row,
#               like a truncated gzip (D45), instead of a silently short file.
# WRITES        <sample>.normalise.tsv (file, status, problem: the same columns
#               as check_fastq.sh) and <sample>_<n>.fastq.gz for each input.
# EXIT          Always 0: a damaged file is a result, not an error (D45).
#
# Usage:  normalise_fastq.sh <sample> <threads> <file1> [file2]
#*******************************************************************************
set -uo pipefail
sample=$1; threads=$2; shift 2
out="${sample}.normalise.tsv"
printf 'file\tstatus\tproblem\n' > "$out"
n=0
for f in "$@"; do
    n=$((n + 1))
    dest="${sample}_${n}.fastq.gz"
    case $f in
        *.gz)  ln -s "$f" "$dest"; printf '%s\tok\t\n' "$f" >> "$out"; continue ;;
        *.bz2) dec=(bzip2 -dc) ;;
        *.xz)  dec=(xz -dc) ;;
        *.zst) dec=(zstd -dcq) ;;
        *)     dec=(cat) ;;
    esac
    if err=$( { "${dec[@]}" "$f" | bgzip -c -@ "$threads" > "$dest"; } 2>&1 ); then
        printf '%s\tok\t\n' "$f" >> "$out"
    else
        first=$(printf '%s\n' "$err" | grep -m1 -v '^$' | tr '\t' ' ')
        printf '%s\tunreadable\t%s\n' "$f" "could not be decompressed (${dec[0]}): ${first:-the file is truncated or corrupt}" >> "$out"
        rm -f -- "$dest"
    fi
done
