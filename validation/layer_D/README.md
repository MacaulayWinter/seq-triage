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
