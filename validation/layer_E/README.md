# Validation layer E: sample type (isolate, metagenome, amplicon)

Does seq-triage call the library type right? Truth is the library record in
ENA: `WGS` + `GENOMIC` = isolate, `WGS` + `METAGENOMIC` = metagenome,
`AMPLICON` = amplicon. Design and pass mark (fixed before any run): CLAUDE.md
D65 and [`docs/validation.md`](../../docs/validation.md).

| File | What it is |
|---|---|
| `strata.psv` | The 7 strata (57 runs): ENA query, platform, true type |
| `../outcomes/layer_E/panel.tsv` | The panel, drawn and committed before any reads were fetched |
| `run_layer_E.sh` | Runs everything on the HPC, start to finish |
| `score_layer_E.R` | Scores the results (on your Mac) |

Isolate vs metagenome needs the contamination screen, so this runs on the HPC
with a Kraken2 database, using the **layer D settings file**.

**1. Copy the code to the HPC** (on your Mac):

```bash
rsync -av --exclude work --exclude .nextflow --exclude .nf-test --exclude 'results' --exclude '.nextflow.log*' ~/Documents/Bioinformatics_pipelines/seq-triage isca:~/
```

**2. Submit** (on the HPC, from `~/seq-triage`):

```bash
sbatch -p YOUR_PARTITION -A YOUR_ACCOUNT validation/layer_E/run_layer_E.sh validation/layer_D/settings.config ~/seqtriage_layer_E
```

**3. Bring the results back and score** (on your Mac):

```bash
mkdir -p ~/seqtriage-scratch/validation/layer_E && scp isca:~/seqtriage_layer_E/layer_E_results.tgz ~/seqtriage-scratch/validation/layer_E/ && tar -xzf ~/seqtriage-scratch/validation/layer_E/layer_E_results.tgz -C ~/seqtriage-scratch/validation/layer_E
```

```bash
cd ~/Documents/Bioinformatics_pipelines/seq-triage && Rscript validation/layer_E/score_layer_E.R validation/outcomes/layer_E/panel.tsv ~/seqtriage-scratch/validation/layer_E layer_E_scores.tsv
```

The results file holds tables only: verdicts, metrics, detected types,
software versions, the database's name and size, and checksums of the reads.

## Round 2 (D68)

Same steps, with the round 2 panel. On the HPC:

```bash
cd ~/seq-triage && PANEL=validation/outcomes/layer_E_round2/panel.tsv sbatch -p YOUR_PARTITION -A YOUR_ACCOUNT --export=ALL validation/layer_E/run_layer_E.sh validation/layer_D/settings.config ~/seqtriage_layer_E_r2
```

Score on your Mac with `validation/outcomes/layer_E_round2/panel.tsv` and
`~/seqtriage_layer_E_r2/layer_E_results.tgz`.

## Round 4 (D87)

Tests D86 and the library-structure rules. Same steps, with the round 4
panel and SortMeRNA's rRNA database. On the HPC:

```bash
cd ~/seq-triage && PANEL=validation/outcomes/layer_E_round4/panel.tsv EXTRA_ARGS="--rrna_db https://github.com/sortmerna/sortmerna/releases/download/v7.0.0/smr_v4.3_fast_db.fasta.gz" sbatch -p YOUR_PARTITION -A YOUR_ACCOUNT --export=ALL validation/layer_E/run_layer_E.sh validation/layer_D/settings.config ~/seqtriage_layer_E_r4
```

Score on your Mac:

```bash
Rscript validation/layer_E/score_layer_E_round4.R validation/outcomes/layer_E_round4/panel.tsv ~/seqtriage-scratch/validation/layer_E_r4_results layer_E_r4
```

## Round 7 (D94)

G-dominated tiny runs only, chosen by a fixed rule; run on the Mac (the runs
are tiny). From the repository root:

```bash
python3 validation/layer_E/round7/screen_pool.py <pool_panel.tsv> validation/outcomes/layer_E_round7/panel.tsv
bash validation/fetch_panel.sh validation/outcomes/layer_E_round7/panel.tsv ~/seqtriage-scratch/validation/layer_E_r7/panel
nextflow run . -profile docker --input ~/seqtriage-scratch/validation/layer_E_r7/panel/samplesheet.csv --outdir ~/seqtriage-scratch/validation/layer_E_r7/results
Rscript validation/layer_E/score_layer_E_round7.R validation/outcomes/layer_E_round7/panel.tsv ~/seqtriage-scratch/validation/layer_E_r7 layer_E_r7
```

## Round 6 (D93)

After round 5's mark C failed. Same steps; the bisulfite share is the fourth
argument of the scorer. On the HPC:

```bash
cd ~/seq-triage && PANEL=validation/outcomes/layer_E_round6/panel.tsv sbatch -p YOUR_PARTITION -A YOUR_ACCOUNT --export=ALL validation/layer_E/run_layer_E.sh validation/layer_D/settings.config ~/seqtriage_layer_E_r6
```

Bring back and score on your Mac:

```bash
mkdir -p ~/seqtriage-scratch/validation/layer_E_r6_results && scp isca:~/seqtriage_layer_E_r6/layer_E_results.tgz ~/seqtriage-scratch/validation/layer_E_r6_results/ && tar -xzf ~/seqtriage-scratch/validation/layer_E_r6_results/layer_E_results.tgz -C ~/seqtriage-scratch/validation/layer_E_r6_results
```

```bash
Rscript validation/layer_E/score_layer_E_round4.R validation/outcomes/layer_E_round6/panel.tsv ~/seqtriage-scratch/validation/layer_E_r6_results layer_E_r6 0.75
```

## Round 5 (D90)

Confirms the bisulfite fix on unseen runs, after round 4's mark B failed. On
the HPC:

```bash
cd ~/seq-triage && PANEL=validation/outcomes/layer_E_round5/panel.tsv sbatch -p YOUR_PARTITION -A YOUR_ACCOUNT --export=ALL validation/layer_E/run_layer_E.sh validation/layer_D/settings.config ~/seqtriage_layer_E_r5
```

Bring back and score on your Mac:

```bash
mkdir -p ~/seqtriage-scratch/validation/layer_E_r5_results && scp isca:~/seqtriage_layer_E_r5/layer_E_results.tgz ~/seqtriage-scratch/validation/layer_E_r5_results/ && tar -xzf ~/seqtriage-scratch/validation/layer_E_r5_results/layer_E_results.tgz -C ~/seqtriage-scratch/validation/layer_E_r5_results
```

```bash
Rscript validation/layer_E/score_layer_E_round4.R validation/outcomes/layer_E_round5/panel.tsv ~/seqtriage-scratch/validation/layer_E_r5_results layer_E_r5
```

## Round 3 (D74)

Tests D73. Same steps, with the round 3 panel. On the HPC:

```bash
cd ~/seq-triage && PANEL=validation/outcomes/layer_E_round3/panel.tsv sbatch -p YOUR_PARTITION -A YOUR_ACCOUNT --export=ALL validation/layer_E/run_layer_E.sh validation/layer_D/settings.config ~/seqtriage_layer_E_r3
```

Score on your Mac with `validation/outcomes/layer_E_round3/panel.tsv` and
`~/seqtriage_layer_E_r3/layer_E_results.tgz`.
