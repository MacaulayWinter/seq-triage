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
