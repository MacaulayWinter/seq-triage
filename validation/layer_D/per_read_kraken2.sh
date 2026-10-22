#!/usr/bin/env bash
#SBATCH --job-name=k2_perread
#SBATCH --time=4:00:00
#SBATCH --mem=110G
#SBATCH --cpus-per-task=8
#SBATCH --output=k2_perread_%j.log
#*******************************************************************************
# per_read_kraken2.sh: Kraken2's call for every read of one FASTQ, against one
#                      or more databases (CLAUDE.md D56 follow-up)
#
# WHY           Layer D round 2: the strain mixture on nanopore showed 6.4%
#               K. pneumoniae with the full Standard database and 0.0% with
#               Standard-8. Competitive mapping on the Mac showed the source
#               reads (AH21H052) are E. coli (only 5 of 3,600 align well to
#               K. pneumoniae), so the question is WHICH reads the full
#               database calls K. pneumoniae. This writes one row per read:
#               read ID and taxon ID. Public reads, so the table holds no
#               private data; it is small (about 0.3 MB per database).
# WHAT IT DOES  Builds a conda environment with the pipeline's Kraken2 version
#               (2.1.6) once, reads each database into the file cache, then
#               classifies the reads with the same options as the pipeline
#               (--memory-mapping; no confidence threshold).
# WRITES        <out_prefix>.<db name>.read_taxid.tsv.gz and .report.txt for
#               each database.
#
# Usage:  sbatch -p <partition> -A <account> per_read_kraken2.sh <reads.fastq.gz> <out_prefix> <db_dir> [<db_dir> ...]
#*******************************************************************************
set -euo pipefail
reads=$1; out=$2; shift 2
# Check every database folder first, so a wrong path stops the job before
# any database is loaded (the first submission read a settings comment into
# the path and failed with a bare cat error).
for db in "$@"; do
    if [ ! -f "$db/hash.k2d" ]; then
        echo "Not a Kraken2 database folder (no hash.k2d): [$db]" >&2
        exit 1
    fi
done
env=${K2_ENV:-$HOME/.conda/envs/seqtriage_k2}
if [ ! -x "$env/bin/kraken2" ]; then
    conda create -y -q -p "$env" -c conda-forge -c bioconda kraken2=2.1.6
fi
"$env/bin/kraken2" --version | head -1
for db in "$@"; do
    name=$(basename "$db")
    echo "== $name: loading into the file cache"
    cat "$db"/*.k2d > /dev/null
    "$env/bin/kraken2" --db "$db" --memory-mapping --threads "${SLURM_CPUS_PER_TASK:-8}" \
        --output "$out.$name.perread.tsv" --report "$out.$name.report.txt" "$reads"
    # Keep read ID and taxon ID only (columns 2 and 3 of Kraken2's output).
    cut -f2,3 "$out.$name.perread.tsv" | gzip -n > "$out.$name.read_taxid.tsv.gz"
    rm "$out.$name.perread.tsv"
done
echo "Done: $(ls "$out".*.read_taxid.tsv.gz | tr '\n' ' ')"
