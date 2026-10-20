#!/usr/bin/env bash
#*******************************************************************************
# fetch_sources.sh: download the start of every layer D source run (D43)
#
# Streams only the FIRST reads of each run from ENA: 100,000 pairs per
# Illumina run, 12,000 reads per nanopore run; the rest is never downloaded.
# Already-fetched files are reused. Checksums of the uncompressed reads go to
# <out_dir>/fetch_md5.tsv, so a run can be compared with any other.
#
# Usage:  bash validation/layer_D/fetch_sources.sh <sources.tsv> <out_dir>
#*******************************************************************************
set -u
sources=$1
out=$2
mkdir -p "$out"
ILM_PAIRS=100000
ONT_READS=12000
md5cmd() { if command -v md5sum >/dev/null; then md5sum | cut -d' ' -f1; else md5 -q; fi; }
printf 'source\tplatform\tfile\treads\tmd5_uncompressed\n' > "$out/fetch_md5.tsv"

# One file: first $2 reads of ENA URL $1 into $3 (gzipped). curl stops with
# "error 23" once head has enough; that is expected.
fetch() {
    if [ ! -s "$3" ]; then
        curl -sS --retry 3 --max-time 1800 "https://$1" 2>/dev/null | gunzip -c 2>/dev/null \
            | head -n $(( $2 * 4 )) | gzip -n > "$3.tmp" && mv "$3.tmp" "$3"
    fi
    local n; n=$(( $(gunzip -c "$3" | wc -l) / 4 ))
    printf '%s\t%s\t%s\t%s\t%s\n' "$4" "$5" "$(basename "$3")" "$n" "$(gunzip -c "$3" | md5cmd)" >> "$out/fetch_md5.tsv"
    [ "$n" -gt 0 ] || { echo "ERROR: no reads fetched from $1" >&2; exit 1; }
    echo "  $(basename "$3"): $n reads"
}
urls() { curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$1&result=read_run&fields=fastq_ftp&format=tsv" | tail -1 | cut -f2; }

# Process substitution (not a pipe), so "exit 1" inside the loop stops the script.
while IFS=$'\t' read -r src _org _iso _study ilm ont _role; do
    echo "$src"
    IFS=';' read -r -a f <<< "$(urls "$ilm")"
    r1=$(printf '%s\n' "${f[@]}" | grep '_1\.fastq\.gz$'); r2=$(printf '%s\n' "${f[@]}" | grep '_2\.fastq\.gz$')
    fetch "$r1" $ILM_PAIRS "$out/${src}_ILM_R1.fastq.gz" "$src" illumina
    fetch "$r2" $ILM_PAIRS "$out/${src}_ILM_R2.fastq.gz" "$src" illumina
    fetch "$(urls "$ont" | cut -d';' -f1)" $ONT_READS "$out/${src}_ONT.fastq.gz" "$src" ont
done < <(tail -n +2 "$sources")
