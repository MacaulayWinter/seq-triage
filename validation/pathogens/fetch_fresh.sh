#!/usr/bin/env bash
# fetch_fresh.sh  (D96, P7): the first 100,000 read pairs of each run with the
# given role in fresh/panel.tsv, from ENA's FTP over HTTPS, and a samplesheet
# (declared isolate, Illumina). Usage: fetch_fresh.sh <panel.tsv> <role> <out_dir>
set -euo pipefail
panel=$1; role=$2; out=$3
mkdir -p "$out"
sheet="$out/samplesheet_${role}.csv"
echo "sample,fastq_1,fastq_2,platform,sample_type" > "$sheet"
awk -F'\t' -v r="$role" 'NR==1{for(i=1;i<=NF;i++)h[$i]=i; next} $h["role"]==r{print $h["run_accession"]"\t"$h["fastq_ftp"]}' "$panel" |
while IFS=$'\t' read -r run ftp; do
  IFS=';' read -ra urls <<< "$ftp"
  [ "${#urls[@]}" -ge 2 ] || { echo "$run: not two files, skipped" >&2; continue; }
  n=0
  for u in "${urls[@]: -2}"; do
    n=$((n+1))
    f="$out/${run}_R${n}.fastq.gz"
    [ -s "$f" ] || { curl -sf --retry 2 "https://${u}" | gzip -dc 2>/dev/null | head -n 400000 | gzip > "$f" || true; }
  done
  echo "${run},$(pwd)/$out/${run}_R1.fastq.gz,$(pwd)/$out/${run}_R2.fastq.gz,illumina,isolate" >> "$sheet"
done
echo "fetched: $(($(wc -l < "$sheet") - 1)) runs"
