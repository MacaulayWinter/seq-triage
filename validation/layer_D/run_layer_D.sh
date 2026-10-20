#!/usr/bin/env bash
#SBATCH --job-name=seqtriage-layerD
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=layer_D_%j.log
#*******************************************************************************
# run_layer_D.sh: validation layer D, start to finish (CLAUDE.md D43)
#
# On an HPC (the user runs this; Claude never acts on the HPC):
#   sbatch -p <partition> -A <account> validation/layer_D/run_layer_D.sh <settings.config> <work_dir>
# This job is only the conductor: Nextflow submits each pipeline step as its
# own SLURM job, using the account and partition in <settings.config>.
#
# Locally, as a rehearsal (Docker, tiny test database):
#   PROFILE=docker,emulate_amd64 bash validation/layer_D/run_layer_D.sh <settings.config> <work_dir>
#
# Steps: 1 fetch the start of each source run from ENA; 2 build the mixtures;
#        3 run seq-triage on them; 4 record which database was used;
#        5 pack ONLY the result tables into <work_dir>/layer_D_results.tgz
#          (no reads), to copy back for scoring.
# Run from the repository root.
#*******************************************************************************
set -euo pipefail

settings=$(realpath "$1")
work=$(mkdir -p "$2" && realpath "$2")
PROFILE=${PROFILE:-conda,slurm}
repo=$(pwd)
[ -f "$repo/main.nf" ] || { echo "Run this from the seq-triage repository root." >&2; exit 1; }
echo "seq-triage layer D | profile: $PROFILE | work: $work | commit: $(git -C "$repo" rev-parse --short HEAD)"

# Conda environments: the launch environment (Nextflow, Java) and a tiny
# R environment for the mixture builder. Created on first use.
if [[ "$PROFILE" == *conda* ]]; then
    eval "$(conda shell.bash hook)"
    conda env list | grep -q '^seq-triage-dev ' || conda env create -f "$repo/environment.yml"
    conda env list | grep -q '^seq-triage-validation ' || conda env create -f "$repo/validation/layer_D/environment.yml"
    conda activate seq-triage-validation
    RSCRIPT=$(command -v Rscript)
    conda activate seq-triage-dev
    export NXF_CONDA_CACHEDIR=${NXF_CONDA_CACHEDIR:-$work/conda-envs}
else
    RSCRIPT=${RSCRIPT:-Rscript}
fi

echo "== 1. Fetch source reads"
bash "$repo/validation/layer_D/fetch_sources.sh" "$repo/validation/layer_D/sources.tsv" "$work/sources"

echo "== 2. Build mixtures"
"$RSCRIPT" "$repo/validation/layer_D/build_mixtures.R" "$work/sources" "$repo/validation/layer_D/mixtures.tsv" "$work/mixtures"

echo "== 3. Run seq-triage"
nextflow run "$repo/main.nf" -profile "$PROFILE" -c "$settings" \
    --input "$work/mixtures/samplesheet.csv" --outdir "$work/results" -w "$work/nf-work" -resume

echo "== 4. Record the database"
# Read the path straight from the settings file (its format is ours). Not
# via `nextflow config`: that command does not accept -c (found in the
# local rehearsal, which then stopped here silently).
db=$(sed -n "s/^[[:space:]]*params\.kraken2_db[[:space:]]*=[[:space:]]*['\"]\([^'\"]*\)['\"].*/\1/p" "$settings")
[ -n "$db" ] || { echo "ERROR: no params.kraken2_db line found in $settings" >&2; exit 1; }
{
    printf 'kraken2_db_name\t%s\n' "$(basename "$db")"
    if [ -d "$db" ]; then
        printf 'hash_k2d_bytes\t%s\n' "$(stat -c %s "$db/hash.k2d" 2>/dev/null || stat -f %z "$db/hash.k2d")"
        printf 'libraries\t%s\n' "$( [ -f "$db/library_report.tsv" ] && cut -f1 "$db/library_report.tsv" | grep -v '^#' | sort | uniq -c | awk '{printf "%s:%s ", $2, $1}')"
        printf 'files_dated\t%s\n' "$(date -r "$db/hash.k2d" +%Y-%m-%d 2>/dev/null || true)"
    fi
} > "$work/results/database_info.tsv"
cat "$work/results/database_info.tsv"

echo "== 5. Pack result tables (no reads)"
tar -czf "$work/layer_D_results.tgz" -C "$work" \
    results/triage_summary.tsv results/triage_metrics.tsv results/qc_values.tsv \
    results/database_info.tsv results/pipeline_info/software_versions.tsv \
    mixtures/truth.tsv sources/fetch_md5.tsv
ls -lh "$work/layer_D_results.tgz"
echo "Done. Copy layer_D_results.tgz back to your Mac for scoring."
