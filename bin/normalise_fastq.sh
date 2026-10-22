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
# FASTA READS   (D79) .fa/.fasta/.fna (any of the compressions) are turned
#               into FASTQ: multi-line sequences joined, and every base given
#               the same placeholder quality ("I"). Detection then sees a
#               single quality value and flags QUALITY_UNINFORMATIVE, so
#               quality is never judged and the sample cannot PASS (D26, D40).
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
        *.gz)  dec=(gzip -dc) ;;
        *.bz2) dec=(bzip2 -dc) ;;
        *.xz)  dec=(xz -dc) ;;
        *.zst) dec=(zstd -dcq) ;;
        *)     dec=(cat) ;;
    esac
    base=${f%.gz}; base=${base%.bz2}; base=${base%.xz}; base=${base%.zst}
    case $base in
        *.fa|*.fasta|*.fna) fasta=1 ;;
        *)                  fasta=0 ;;
    esac
    if [ "$fasta" -eq 0 ] && [[ $f == *.gz ]]; then
        ln -s "$f" "$dest"; printf '%s\tok\t\n' "$f" >> "$out"; continue
    fi
    # FASTA -> FASTQ: one record per header, sequence lines joined, a
    # placeholder quality of the same length (D79).
    conv=(cat)
    [ "$fasta" -eq 1 ] && conv=(awk 'function flush() { if (h != "") { q = s; gsub(/./, "I", q); print "@" h "\n" s "\n+\n" q } }
        /^>/ { flush(); h = substr($0, 2); sub(/\r$/, "", h); s = ""; next } { gsub(/[ \t\r]/, ""); s = s $0 } END { flush() }')
    if err=$( { "${dec[@]}" "$f" | "${conv[@]}" | bgzip -c -@ "$threads" > "$dest"; } 2>&1 ); then
        printf '%s\tok\t\n' "$f" >> "$out"
    else
        first=$(printf '%s\n' "$err" | grep -m1 -v '^$' | tr '\t' ' ')
        printf '%s\tunreadable\t%s\n' "$f" "could not be decompressed (${dec[0]}): ${first:-the file is truncated or corrupt}" >> "$out"
        rm -f -- "$dest"
    fi
done
