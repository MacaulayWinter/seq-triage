#!/usr/bin/env bash
#*******************************************************************************
# kraken2_batch.sh: classify many samples with one database load (D50)
#
# Kraken2 loads its database afresh each time it runs, and a big database
# (full Standard: ~94 GB) then dominates the run time and the I/O. Here the
# database is read ONCE into the node's file cache, and every sample is
# classified with --memory-mapping, so Kraken2 uses that cached copy instead
# of loading its own. Same command and output names as the nf-core
# kraken2/kraken2 module, so the two modes are interchangeable.
#
# Usage:  kraken2_batch.sh <db_dir> <threads> <manifest.tsv>
#         manifest: sample <TAB> single_end (true|false) <TAB> file[,file2]
#*******************************************************************************
set -euo pipefail
db=$1; threads=$2; manifest=$3

echo "Reading the database into the file cache once: $(du -sh "$db" | cut -f1)"
cat "$db"/*.k2d > /dev/null

while IFS=$'\t' read -r sample single_end files; do
    [ -n "$sample" ] || continue
    paired=""; [ "$single_end" = "true" ] || paired="--paired"
    echo "Classifying $sample"
    # shellcheck disable=SC2086  (the file list must split on commas -> words)
    kraken2 \
        --db "$db" \
        --memory-mapping \
        --threads "$threads" \
        --report "${sample}.kraken2.report.txt" \
        --gzip-compressed \
        --output /dev/null \
        $paired \
        ${files//,/ }
done < "$manifest"
