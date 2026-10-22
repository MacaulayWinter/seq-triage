#!/usr/bin/env bash
#SBATCH --job-name=seqtriage-layerE
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=layer_E_%j.log
#*******************************************************************************
# run_layer_E.sh: validation layer E (sample type), start to finish (D65)
#
# On an HPC (the user runs this; Claude never acts on the HPC):
#   sbatch -p <partition> -A <account> validation/layer_E/run_layer_E.sh <settings.config> <work_dir>
# <settings.config> is the layer D one: it names the Kraken2 database, the
# partition and the account. This job only conducts: Nextflow submits each
# pipeline step as its own SLURM job.
#
# Steps: 1 fetch the start of each panel run from ENA, declaring ENA's
#        platform (sample_type is NOT declared: the reads must decide);
#        2 run seq-triage; 3 record the database; 4 pack ONLY the result
#        tables into <work_dir>/layer_E_results.tgz (no reads).
# Run from the repository root.
#*******************************************************************************
set -euo pipefail

settings=$(realpath "$1")
work=$(mkdir -p "$2" && realpath "$2")
PROFILE=${PROFILE:-conda,slurm}
repo=$(pwd)
[ -f "$repo/main.nf" ] || { echo "Run this from the seq-triage repository root." >&2; exit 1; }
echo "seq-triage layer E | profile: $PROFILE | work: $work | commit: $(git -C "$repo" rev-parse --short HEAD)"

if [[ "$PROFILE" == *conda* ]]; then
    eval "$(conda shell.bash hook)"
    conda env list | grep -q '^seq-triage-dev ' || conda env create -f "$repo/environment.yml"
    conda activate seq-triage-dev
    export NXF_CONDA_CACHEDIR=${NXF_CONDA_CACHEDIR:-$work/conda-envs}
fi

echo "== 1. Fetch panel reads"
DECLARE=ena bash "$repo/validation/fetch_panel.sh" "$repo/validation/outcomes/layer_E/panel.tsv" "$work/panel"

echo "== 2. Run seq-triage"
# Own launch folder, so -resume cannot pick up another run (layer D, D50).
( cd "$work" && nextflow run "$repo/main.nf" -profile "$PROFILE" -c "$settings" \
    --input "$work/panel/samplesheet.csv" --outdir "$work/results" -w "$work/nf-work" -resume )

echo "== 3. Record the database"
db=$(sed -n "s/^[[:space:]]*params\.kraken2_db[[:space:]]*=[[:space:]]*['\"]\([^'\"]*\)['\"].*/\1/p" "$settings")
[ -n "$db" ] || { echo "ERROR: no params.kraken2_db line found in $settings" >&2; exit 1; }
{
    printf 'kraken2_db_name\t%s\n' "$(basename "$db")"
    printf 'hash_k2d_bytes\t%s\n' "$(stat -c %s "$db/hash.k2d" 2>/dev/null || stat -f %z "$db/hash.k2d")"
    printf 'files_dated\t%s\n' "$(date -r "$db/hash.k2d" +%Y-%m-%d 2>/dev/null || true)"
} > "$work/results/database_info.tsv"

# qc_values.tsv carries every Kraken2 row, so a later rule can be re-scored
# from the tables alone (the first layer E run lacked it, D67).
echo "== 4. Pack result tables (no reads)"
tar -czf "$work/layer_E_results.tgz" -C "$work" \
    results/triage_summary.tsv results/triage_metrics.tsv results/qc_values.tsv results/database_info.tsv \
    results/pipeline_info/software_versions.tsv panel/fetch_log.tsv
ls -lh "$work/layer_E_results.tgz"
echo "Done. Copy layer_E_results.tgz back to your Mac for scoring."
