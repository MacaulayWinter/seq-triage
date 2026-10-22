# seq-triage

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A526.10.0-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![CI](https://github.com/MacaulayWinter/seq-triage/actions/workflows/ci.yml/badge.svg)](https://github.com/MacaulayWinter/seq-triage/actions/workflows/ci.yml)

A Nextflow DSL2 pipeline that inspects incoming sequencing data before any
analysis starts. It validates the sample metadata, works out from the reads
themselves which platform produced them, routes each sample to the right QC,
and writes a plain-language summary that says which samples are ready and why.

> **Status: v0.1.0 released (2026-10-22); validated (round 3 met every mark).** Every step runs
> end to end and is tested in CI on every change. Round 1 (40 public runs,
> 29 damaged samples) passed all of its pre-set pass marks. Round 2 (115
> runs, 31 damaged samples) passed layers B and C and failed layer A on 2 of
> 115 runs. One was an ENA labelling error; the other led to a stricter rule
> (D58), chosen after the fact. Round 3 (54 runs from unseen studies,
> pre-registered) then tested it and met every pass mark
> ([results](docs/validation.md)). BAM input (v0.2) is built; see
> [Roadmap](#roadmap).

## Why

Sequencing data often arrives with incomplete or wrong metadata: a run marked
Illumina that is really nanopore, a paired-end sample listed with one file, a
sample ID that doesn't match its files. Analysing such data with the wrong
assumptions produces results that look plausible and are wrong.

seq-triage checks the metadata against the reads before anything else runs.
Its rule is conservative: **a sample is only marked PASS when the reads
independently confirm its metadata.** When the evidence is weak, the report
says "uncertain" rather than guessing.

## How it works

```
samplesheet ─► validate (nf-schema) ─┬─► read statistics (seqkit) ─► detect platform
                                     └─► screen metadata ───────────► flag inconsistencies
                                                    ▼
                     route each sample by confidence and agreement
                 Illumina ─► fastp + FastQC        ONT ─► NanoPlot
                                                    ▼
                         MultiQC  +  triage summary report
```

| Detection | Declared metadata | QC routed by | Best possible verdict |
|---|---|---|---|
| high confidence | agrees | detected | PASS |
| high confidence | missing | detected | WARN (inferred, not declared) |
| high confidence | contradicts | detected | WARN (metadata contradicts reads) |
| low or medium, same platform or uncertain | present | declared | WARN (not confirmed by reads) |
| low or medium, other platform | contradicts | not routed | UNCERTAIN (review: neither is certain enough to pick the QC) |
| uncertain | missing | not routed | UNCERTAIN |
| PacBio HiFi (read names ending `/ccs`, or a BAM written by `ccs`) | as for any platform | detected | PASS if declared `pacbio`, else WARN; NanoPlot QC with HiFi limits |
| PacBio other than HiFi, MGI/DNBSEQ or Element (named) | any | not routed | UNCERTAIN (platform not handled) |
| none at high confidence | `pacbio` | not routed | UNCERTAIN (only HiFi is QC'd) |
| empty file | any | not routed | FAIL |

**How the platform is detected.** Read lengths make the call: Illumina reads
are short and nearly all the same length; nanopore reads are long and vary
widely. Independent evidence then sets the confidence:

- **read names**, which each instrument writes in its own fixed format;
- **library adapters**, synthetic sequences each platform's kits attach to
  the DNA, often still present in some reads;
- **quality values**: NovaSeq and NextSeq use only a handful of distinct scores.

PacBio, MGI/DNBSEQ and Element are recognised by their read names (and
PacBio and MGI by their adapters) and reported by name, without QC.
Illumina adapters never raise confidence, because Element and MGI
instruments can sequence converted Illumina libraries. The rules and their
calibration are in CLAUDE.md (D33, D40); how well they work on data they
were not tuned on is in [docs/validation.md](docs/validation.md).

QC metrics can lower a verdict but never raise it. If quality cannot be
assessed (every base has the same score), a sample cannot PASS. Every metric in the report
comes with an explanation of what it measures and how to read it, from
[`assets/metric_guide.tsv`](assets/metric_guide.tsv). The thresholds there
are provisional defaults and are labelled as such.

## Requirements

- [conda](https://conda-forge.org/download/) (Miniforge recommended)
- [Docker](https://www.docker.com/) (Docker Desktop on macOS)

Everything else is installed by the environment file below.

## Installation

```bash
git clone https://github.com/MacaulayWinter/seq-triage.git
cd seq-triage
conda env create -f environment.yml
conda activate seq-triage-dev
```

This installs Java 17, Nextflow 26.10.2 and nf-core tools 4.1.0, with every
version pinned. For a byte-identical environment (every dependency, with
checksums), use the lockfile for your system instead:

```bash
conda create -n seq-triage-dev --file envs/seq-triage-dev.linux-64.lock.txt        # Linux and WSL2 (Intel/AMD)
conda create -n seq-triage-dev --file envs/seq-triage-dev.linux-aarch64.lock.txt   # Linux on ARM
conda create -n seq-triage-dev --file envs/seq-triage-dev.osx-arm64.lock.txt       # Apple Silicon Mac
``` The analysis tools themselves run in version-pinned
containers and need no installation. See
[docs/reproducibility.md](docs/reproducibility.md) for how each layer is
pinned.

## Usage

```bash
nextflow run main.nf -profile docker --input samplesheet.csv --outdir results
```

`--help` lists every parameter with its description and default;
`--help <name>` explains one in full.

**On Apple Silicon Macs**, add `arm64`:

```bash
nextflow run main.nf -profile docker,arm64 --input samplesheet.csv --outdir results
```

seqkit, fastp, FastQC, MultiQC and the R steps then run natively. NanoPlot,
Kraken2, Bracken and untar have no Apple Silicon images yet, so they still
run under Rosetta emulation. Enable it first: Docker Desktop → Settings →
General → *Use Rosetta for x86_64/amd64 emulation on Apple Silicon*. Put
`arm64` after `docker`.

On the test profile this is 32% faster than emulating everything (99 s
against 145 s), and the results are **identical**: verdicts, metrics and
every raw QC number. `-profile docker,emulate_amd64` (everything emulated)
still works.

If the first run stops with Docker's `failed to extract layer … not found`,
two tasks tried to download the same new image at once. Pull it once by
hand, then rerun with `-resume`:

```bash
docker pull --platform linux/arm64 docker.io/rocker/r-ver:4.5.2@sha256:fd4ccdd3a4a6f7ef805e2daeee2a0fe3bf126bc231f36351223baecf5a595a4c
```

**On Windows**, use WSL2 (Windows Subsystem for Linux). Nextflow and most
bioinformatics tools do not run on Windows itself, but WSL2 is a full Linux
built into Windows 10 (version 2004 or later) and 11, and inside it
seq-triage runs exactly as on any Linux machine:

1. In PowerShell as administrator, install WSL2 with Ubuntu, then restart:
   `wsl --install -d Ubuntu-24.04`
2. Install [Docker Desktop](https://docs.docker.com/desktop/setup/install/windows-install/)
   with the WSL 2 backend, and switch on *Settings → Resources → WSL
   integration* for Ubuntu.
3. Open Ubuntu, install [Miniforge](https://conda-forge.org/download/), and
   continue with [Installation](#installation) and the Linux commands above,
   **inside Ubuntu**.

Keep the repository and your data in the Linux file system (`~/...`), not
under `/mnt/c/...`: Windows drives seen from WSL2 are much slower and do not
handle Linux file permissions. These steps follow the official Microsoft,
Docker and Nextflow instructions; they have not been tried on a Windows
machine by this project, so please report any step that needs changing.

**On an HPC or any machine without Docker** (no Apptainer/Singularity
either), use conda for the tools and SLURM for the jobs:

```bash
nextflow run main.nf -profile conda,slurm -c my_cluster.config --input samplesheet.csv --outdir results
```

where `my_cluster.config` holds your own cluster settings, which never
belong in this repository:

```groovy
process.queue          = 'my_partition'
process.clusterOptions = '-A my_account'
```

Each tool's environment is built once from its pinned recipe; set
`NXF_CONDA_CACHEDIR` to a shared folder to reuse them between runs. CI runs
the full test with conda as well as Docker, and expects identical verdicts.

### Contamination screen (optional)

Give a Kraken2 database to also check each sample for human DNA and for a
second organism:

```bash
nextflow run main.nf -profile docker --input samplesheet.csv --outdir results \
    --kraken2_db /path/to/kraken2_db          # a directory, or a .tar.gz of one
```

The report then adds the human fraction, the share of reads from species
other than the most common one, the unclassified fraction and the top
species. Illumina samples get Bracken's species-level re-estimate. Results
depend on the database, and the screen **cannot detect contamination by
another strain of the same species**. Without a database the report says
"not screened". On a laptop, add `--kraken2_memory_mapping` to read the
database from disk instead of loading it into memory. With a big database
on an HPC (the full Standard database is ~94 GB), add `--kraken2_batch`:
every sample is classified in one task, so the database is read once
instead of once per sample. Results are identical either way.

### Samplesheet

A CSV with one row per sample:

```csv
sample,fastq_1,fastq_2,bam,platform,sample_type,molecule
ISO_01,ISO_01_R1.fastq.gz,ISO_01_R2.fastq.gz,,illumina,isolate,dna
MTX_01,MTX_01_R1.fastq.gz,MTX_01_R2.fastq.gz,,illumina,metagenome,rna
ONT_01,ONT_01.fastq.gz,,,ont,,
ONT_02,,,ONT_02.bam,ont,,
UNKNOWN_01,UNKNOWN_01.fastq.gz,,,,,
```

| Column | Required | Notes |
|---|---|---|
| `sample` | yes | Unique ID. |
| `fastq_1` | yes, unless `bam` | Gzipped FASTQ. |
| `fastq_2` | no | Second mate for paired-end data. Leave blank for single-end or ONT. |
| `bam` | instead of FASTQ | A BAM file, aligned or unaligned: dorado and PacBio deliver reads as unaligned BAM. Give either FASTQ or one BAM per row. The BAM is checked from header to end, then converted to FASTQ (`samtools collate \| samtools fastq`); paired or single is decided from the reads' own READ1/READ2 flags, and secondary and supplementary alignments are left out, so each read counts once. The header is read too: software that only runs on one platform's data (`@PG` ccs or lima for PacBio, dorado or guppy for nanopore) counts as evidence, like read names; `@RG PL`, which people type, stands in for a blank `platform`, and a samplesheet that disagrees with it is flagged for review. CRAM is not accepted yet. |
| `platform` | no | `illumina`, `ont` or `pacbio`. Leave blank to have it detected. PacBio is QC'd only when the reads are HiFi. |
| `sample_type` | no | `isolate`, `metagenome`, `amplicon` or `unknown`. Amplicon samples are not judged on read N50 or duplication, which amplicon libraries fail by design. |
| `molecule` | no | `dna` or `rna`. Use `rna` for RNA-seq, metatranscriptomes and direct RNA, even though Illumina RNA reads are cDNA. RNA samples are not judged on duplication, and FastQC's base-content and overrepresented-sequence modules are marked as known false alarms. RNA mostly cannot be detected from the reads, so it must be declared; direct-RNA reads written with U are flagged if `molecule` says `dna`. |

Unknown parameters and unknown samplesheet columns are errors, not warnings.
This stops a typo such as `fastq2` from silently turning paired data into
single-end.

### Outputs

| File (in `--outdir`) | What it is |
|---|---|
| `triage_summary.md` | **Start here.** One line per sample with its verdict, then each sample's reasons, its QC numbers with their limits, and a *How to read these numbers* section explaining every metric shown (N50, Q30 and so on) |
| `triage_summary.tsv` | The same verdicts and reasons as a table |
| `triage_metrics.tsv` | One row per sample and metric: value, limits, result and note |
| `qc_values.tsv` | Every raw number the report judged, traced to the tool and file it came from |
| `multiqc/` | MultiQC report over all QC tools |
| `fastp/`, `fastqc/`, `nanoplot/` | Each tool's own reports (no reads are ever copied here) |
| `integrity/` | Whether each file was readable and well-formed FASTQ. A damaged file gets a FAIL row saying exactly what is wrong, and the other samples still finish |
| `metadata/`, `detection/` | Metadata screen flags, and the platform evidence for each sample |
| `pipeline_info/` | Software versions actually used, run timeline, resource report, workflow diagram (`pipeline_dag.html`) |

Verdicts:

- **PASS:** the reads confirm the metadata and every QC number is within its limits.
- **WARN:** usable, but something needs a look; the report says what.
- **FAIL:** too little or too poor data to use.
- **UNCERTAIN:** the platform could not be established, so no QC was run.

## Testing

The `test` profile runs the whole pipeline on ten small public samples
(~130 MB), pinned by commit. Eight come from
[nf-core/test-datasets](https://github.com/nf-core/test-datasets): a
simulated Illumina pair, real ONT amplicon reads, simulated ONT reads with
placeholder quality, Illumina reads mislabelled as ONT, a sample with no
metadata, single-end Illumina, a tiny real Illumina SARS-CoV-2 amplicon pair (so the
contamination screen's species step runs), and PacBio HiFi reads, which seq-triage must
report as *uncertain* rather than force into a platform. Two come from
[seq-triage-test-data](https://github.com/MacaulayWinter/seq-triage-test-data):
real NovaSeq and R10.4.1 nanopore reads of one *E. coli* isolate, subsampled
from ENA and large enough to pass every threshold. Each sample has a known
correct answer.

```bash
nextflow run main.nf -profile test,docker --outdir results              # Linux
nextflow run main.nf -profile test,docker,emulate_amd64 --outdir results # Apple Silicon
Rscript tests/e2e/check_test_profile.R results
```

The check confirms that every sample got a verdict, that each verdict matches
`tests/e2e/expected_verdicts.tsv`, that no reads were copied into `results/`,
that every tool's version was recorded, that every sample sent to QC has
judged metrics, and that MultiQC read every tool's output.

Unit tests, each on synthetic inputs built inside the test (no data needed):

```bash
Rscript tests/screen_metadata/test_screen_metadata.R   # metadata screen
Rscript tests/detect_platform/test_detect_platform.R   # platform detection rules
Rscript tests/triage_report/test_triage_report.R       # QC limits and verdicts
```

**nf-test**, Nextflow's own testing framework (as used by nf-core), tests
the pipeline's processes and the whole workflow:

```bash
nf-test test --profile docker          # Apple Silicon: --profile docker,arm64
```

| Test | Checks |
|---|---|
| `CHECK_FASTQ` | a good pair passes; malformed FASTQ inside valid gzip, and a truncated gzip, are caught with the right message |
| `CHECK_BAM` | an intact unaligned BAM passes; a truncated BAM, and a BAM damaged part-way through (which samtools quickcheck alone passes), are caught |
| `DETECT_PLATFORM` | real Illumina and nanopore reads, through the real seqkit step, are called correctly at high confidence |
| `SCREEN_METADATA` | the synthetic samplesheet gives exactly the saved output |
| whole pipeline | `-profile test` succeeds, and every sample's verdict and every tool's version match the saved snapshot |

A change that alters any verdict or version fails the pipeline test until
the snapshot is deliberately updated (`nf-test test --update-snapshot`), so
such changes are always visible in review.

GitHub Actions runs all of these on every pull request, and changes reach
`main` only after they pass. `tests/survey/check_survey.R` checks detection
against 12 real public samples; it runs locally, because the data stays
outside the repository.

## Validation

Tests check the code does what it was designed to do. Validation checks the
design is right, on data it was never tuned on, against answers known
independently of the pipeline. The pass marks were fixed before any result
was seen. Round 1 (full write-up in [docs/validation.md](docs/validation.md)):

| Layer | Question | Result |
|---|---|---|
| A | Is the platform call right? (40 public ENA runs, 18 kinds of data, chosen by seeded rule) | **0 wrong calls at high confidence**; 29 of 40 correct at high confidence; the rest medium, low or *uncertain* |
| B | Are the reported numbers right? (recomputed independently in base R) | **542 / 542** agree |
| C | Do damaged files get the verdict they should? (29 cases, answers written in advance) | **29 / 29** |

Validation found and fixed a silent failure in a dependency (seqkit reports
success on a truncated gzip), pinned down how NanoPlot averages quality, and
found one run whose ENA instrument label is wrong; that correction is
recorded with its evidence.

**Round 2** (2026-10-21, 115 runs from a whole-ENA draw, 31 damaged
samples, pre-registered): **B 1,158 / 1,158 and C 30 / 30 pass. A fails as
registered**, with 2 wrong calls at high confidence out of 115. One is very
likely another ENA labelling error; the other exposes a real limit: when the
archive has stripped the read names, binned quality alone can raise an MGI
run to "Illumina, high". Both are set out with their evidence in
[docs/validation.md](docs/validation.md). Decided afterwards (D58): the
first is recorded as an ENA error; binned quality alone now gives at most
medium confidence. Re-scored that way, no call is wrong at high confidence,
at the cost of 7 correct calls dropping to medium. A rule chosen after
seeing the data cannot be tested by the same data, so round 2 stays a fail.

**Round 3** (2026-10-22, 54 runs from studies neither earlier round used,
pre-registered to test that rule): **every mark met**, with 0 wrong calls at
high confidence and 0 unhandled runs sent to QC; 411 / 411 recalculations
agree. The new rule stopped one wrong call it had not been designed on: an
archived MGI run that the old rule would have called Illumina with high
confidence.

## Roadmap

**v0.1**

- [x] Reproducible launch environment (`environment.yml` + lockfile)
- [x] Configuration, resource labels and retry policy
- [x] Parameter validation
- [x] Samplesheet validation and metadata screening
- [x] Platform and layout detection with confidence levels
- [x] Routed QC (fastp + FastQC, NanoPlot) and MultiQC
- [x] Triage summary report: QC metrics judged and explained, `molecule` column for RNA
- [x] `-profile test` on small public data, run in GitHub Actions
- [x] Validation round 1: 40 public runs + 29 damaged samples, all pass marks met
- [x] PacBio, MGI/DNBSEQ and Element recognised and named (no QC yet)
- [x] Validation round 2 (115 runs, whole-ENA draw) and round 3 (54 runs, testing the rule round 2 led to)

**v0.2** (planned)

- BAM/uBAM input (current nanopore and PacBio data is delivered as BAM): **built** (D55); CRAM to follow
- QC for PacBio HiFi: **built** (D64); PacBio CLR stays recognised, not QC'd
- rRNA fraction for metatranscriptomes
- Isolate vs metagenome vs amplicon classification
- Host and contaminant screen, validated with human reads spiked into
  bacterial reads at known fractions
- Recommendation of the right downstream nf-core pipeline, with a
  ready-to-run parameters file

## Data policy

No real research data is stored in this repository. Test data comes from
public sources and is referenced by URL. `.gitignore` blocks sequencing file
formats anywhere in the tree as a second safeguard.

### Repository guard

A pre-commit hook refuses any commit that holds sequencing data (by file
type or by content), a file over 1 MB, an HPC account or cluster path, or a
non-noreply author email. Turn it on once per clone:

```bash
git config core.hooksPath .githooks
```

Add your own private strings (account names, usernames), one per line, to
`.githooks/private_patterns.txt`; that file is git-ignored. CI runs the same
check on every tracked file, in case the hook was skipped.

To merge a pull request only once CI passes (GitHub's branch protection is
not available for private repositories on the free plan):

```bash
bash .githooks/merge_when_green.sh <PR number>
```
## Citations

If you use seq-triage, please cite the tools it is built on, listed in
[CITATIONS.md](CITATIONS.md).

## Author

Macaulay Winter ([@MacaulayWinter](https://github.com/MacaulayWinter))

## Licence

[MIT](LICENSE)
