# Validation layer D: contamination screen accuracy

Does the contamination screen report the right amount of host DNA and of a
second organism? Each test sample is a **mixture made on purpose**, so the
true answer is known. Design and pass marks (fixed before any run):
[`mixtures.tsv`](mixtures.tsv), [`score_layer_D.R`](score_layer_D.R) and
CLAUDE.md D43.

| File | What it is |
|---|---|
| `sources.tsv` | The 6 public source runs (5 bacterial isolates from one study, plus human HG002) |
| `mixtures.tsv` | The 13 mixtures, each built on Illumina and nanopore, with expected verdicts |
| `fetch_sources.sh` | Streams the first reads of each source from ENA |
| `build_mixtures.R` | Builds the mixtures and `truth.tsv` |
| `run_layer_D.sh` | Runs everything, start to finish |
| `score_layer_D.R` | Scores the results against the known answers |
| `settings.example.config` | Template for your own settings (account, partition, database) |

## Running it on an HPC (SLURM, conda)

Everything below runs on the HPC, in your own account. Only a small table of
results comes back.

**1. Get the code** (once):

```bash
git clone https://github.com/MacaulayWinter/seq-triage.git
cd seq-triage
```

The repository is private, so git asks you to sign in. Use a GitHub
personal access token, or set up an SSH key on the HPC.

**2. Make your settings file** (once):

```bash
cp validation/layer_D/settings.example.config validation/layer_D/settings.config
nano validation/layer_D/settings.config
```

Fill in your Kraken2 database folder, your partition and your account.
`settings.config` is ignored by git, so it can never be committed.

**3. Submit the run:**

```bash
sbatch -p YOUR_PARTITION -A YOUR_ACCOUNT validation/layer_D/run_layer_D.sh validation/layer_D/settings.config ~/seqtriage_layer_D
```

The first run creates the conda environments, which takes about 20-40
minutes. Then the pipeline runs on 26 samples. Watch progress with
`squeue -u $USER` and `tail -f layer_D_*.log`.

**4. Bring the results back to your Mac.** Run this **on your Mac**, with the
VPN connected:

```bash
scp isca:~/seqtriage_layer_D/layer_D_results.tgz ~/seqtriage-scratch/validation/
```

The file holds result tables only: verdicts, metrics, the run's
software versions, the database's name and size, and checksums of the input
reads. No reads.

**5. Score it**, on your Mac:

```bash
mkdir -p ~/seqtriage-scratch/validation/layer_D_isca && tar -xzf ~/seqtriage-scratch/validation/layer_D_results.tgz -C ~/seqtriage-scratch/validation/layer_D_isca
Rscript validation/layer_D/score_layer_D.R ~/seqtriage-scratch/validation/layer_D_isca layer_D_scores.tsv
```

## If something goes wrong

The run log `layer_D_<jobid>.log` and Nextflow's `.nextflow.log` say where
it stopped. Re-submitting the same command resumes from the last finished
step (`-resume`), and already-downloaded source reads are reused.

## Round 2: new isolates, two databases

Round 2 (CLAUDE.md D49) repeats the same 13 mixtures and pass marks with
**different isolates** ([`round2/`](round2/)). It is the real test of the
1% trace-species rule (D47), which was chosen after seeing round 1. Running
it against both the 8 GB-capped database and the full one shows how much of
the background and the *Enterobacter* under-detection comes from the cap.

**1. Copy the updated code to the HPC** (on your Mac):

```bash
rsync -av --exclude work --exclude .nextflow --exclude .nf-test --exclude 'results' --exclude '.nextflow.log*' ~/Documents/Bioinformatics_pipelines/seq-triage isca:~/
```

`settings.config` lives only on the HPC, and rsync leaves it alone.

**2. Make a second settings file for the full database** (on the HPC). Kraken2
loads the whole database into memory, so it needs about 100 GB:

```bash
cd ~/seq-triage && sed -e "s#^params.kraken2_db .*#params.kraken2_db = '$HOME/kraken2_db_full'   // Standard (full), Oct 2025#" -e 's/memory = 16.GB/memory = 100.GB/' validation/layer_D/settings.config > validation/layer_D/settings_full.config && cat validation/layer_D/settings_full.config
```

**3. Submit both runs:**

```bash
DESIGN=validation/layer_D/round2 sbatch -p YOUR_PARTITION -A YOUR_ACCOUNT --export=ALL validation/layer_D/run_layer_D.sh validation/layer_D/settings.config ~/seqtriage_layer_D_r2_std8
DESIGN=validation/layer_D/round2 sbatch -p YOUR_PARTITION -A YOUR_ACCOUNT --export=ALL validation/layer_D/run_layer_D.sh validation/layer_D/settings_full.config ~/seqtriage_layer_D_r2_full
```

**4. Bring both results back** (on your Mac), then score them as in step 5
above:

```bash
mkdir -p ~/seqtriage-scratch/validation/r2 && scp isca:~/seqtriage_layer_D_r2_std8/layer_D_results.tgz ~/seqtriage-scratch/validation/r2/std8.tgz && scp isca:~/seqtriage_layer_D_r2_full/layer_D_results.tgz ~/seqtriage-scratch/validation/r2/full.tgz
```
