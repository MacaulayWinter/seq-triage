# seq-triage

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A526.10.0-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![CI](https://github.com/MacaulayWinter/seq-triage/actions/workflows/ci.yml/badge.svg)](https://github.com/MacaulayWinter/seq-triage/actions/workflows/ci.yml)

A Nextflow DSL2 pipeline that inspects incoming sequencing data before any
analysis starts. From the reads themselves it works out what they are
(platform, layout, library structure, isolate / metagenome / amplicon,
contamination, rRNA share, depth), checks that against the sample
metadata, runs the right QC, and writes a plain-language report: which
samples are ready and why, and **what can be done with each one next**,
with the nf-core pipeline and a ready-to-fill samplesheet.

> **Status:** v0.2.0. Every change is tested end to end in CI. Validated in
> pre-registered rounds on public data it was not tuned on, failures
> included: see [Validation](#validation).

**Quick start** (needs [conda](https://conda-forge.org/download/) and
[Docker](https://www.docker.com/)):

```bash
git clone https://github.com/MacaulayWinter/seq-triage.git && cd seq-triage
conda env create -f environment.yml && conda activate seq-triage-dev
nextflow run main.nf -profile test,docker --outdir results   # Apple Silicon: test,docker,arm64
```

Then open `results/triage_summary.md`. For your own data, see [Usage](#usage).

**What it reads:** FASTQ (plain, `.gz`, `.bz2`, `.xz`, `.zst`), FASTA reads,
BAM, SAM and CRAM; single, paired or interleaved; one file or several per
sample (lanes, a nanopore barcode folder). **Platforms QC'd:** Illumina,
MGI/DNBSEQ, Element, Oxford Nanopore and PacBio HiFi; other platforms are
named and left UNCERTAIN rather than guessed. **Libraries recognised from
the reads:** isolate, metagenome, amplicon, 10x single-cell, bisulfite /
EM-seq; human and other targets via `target_taxid`.

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
samplesheet ─► validate (nf-schema; rows of one sample must agree) ─► screen metadata
     │
     ▼  per file: convert to .fastq.gz (BAM/SAM/CRAM, FASTA, .bz2/.xz/.zst) ─► integrity check
     ▼  per lane: read statistics (seqkit) ─► detect platform, layout, library structure
     ▼  combine lanes ─► route by confidence and agreement with the metadata
     ▼  merge lanes / split interleaved ─► QC
            short reads (Illumina, MGI/DNBSEQ, Element): fastp + FastQC
            long reads (nanopore, PacBio HiFi): NanoPlot
     ▼  optional screens: Kraken2 + Bracken (contamination, sample type), SortMeRNA (rRNA)
     ▼  triage report: verdicts, explained metrics, sample type, depth   + MultiQC
     ▼  next steps: nf-core pipeline + samplesheet per sample, and the analysis menu
```

| Detection | Declared metadata | QC routed by | Best possible verdict |
|---|---|---|---|
| high confidence | agrees | detected | PASS |
| high confidence | missing | detected | WARN (inferred, not declared) |
| high confidence | contradicts | detected | WARN (metadata contradicts reads) |
| low or medium, same platform or uncertain | present | declared | WARN (not confirmed by reads) |
| low or medium, other read class (short vs long) | contradicts | not routed | UNCERTAIN (review: neither is certain enough to pick the QC) |
| uncertain | missing | not routed | UNCERTAIN |
| PacBio HiFi (read names ending `/ccs`, or a BAM written by `ccs`) | as for any platform | detected | PASS if declared `pacbio`, else WARN; NanoPlot QC with HiFi limits |
| MGI/DNBSEQ or Element (named by read names or adapters) | as for any platform | detected | short-read QC (fastp + FastQC), their own guide rows |
| PacBio other than HiFi (named) | any | not routed | UNCERTAIN (platform not handled) |
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
PacBio and MGI by their adapters). MGI/DNBSEQ and Element get the same
short-read QC as Illumina; PacBio is QC'd only as HiFi. Illumina adapters
never raise confidence, because Element and MGI
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
```

The analysis tools themselves run in version-pinned
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

seqkit, fastp, FastQC, MultiQC, samtools and the R steps then run natively.
NanoPlot, Kraken2, Bracken, untar, SortMeRNA and the file-merging step have
no Apple Silicon images yet, so they still run under Rosetta emulation. Enable it first: Docker Desktop → Settings →
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

### rRNA screen (optional)

For RNA samples (`molecule = rna`), give an rRNA reference to measure how
much of the library is ribosomal RNA, which tells you whether rRNA
depletion worked:

```bash
nextflow run main.nf -profile docker --input samplesheet.csv --outdir results \
    --rrna_db https://github.com/sortmerna/sortmerna/releases/download/v7.0.0/smr_v4.3_fast_db.fasta.gz
```

SortMeRNA (4.3.7) classifies each RNA sample against it; the report adds
`rrna_fraction` (warn above 25%, fail above 75%, provisional; not judged
for amplicons, which are rRNA by design). `--rrna_db` takes one FASTA or a
`.txt` manifest of several (paths or URLs); the index is built once per
run. The `smr_v4.3_*` databases match the SortMeRNA version the pipeline
runs. Without `--rrna_db`, RNA samples say "not screened".

### Samplesheet

A CSV with one row per sample, or several rows with the same `sample` for
a sample whose reads are in several files:

```csv
sample,fastq_1,fastq_2,bam,fastq_dir,platform,sample_type,molecule
ISO_01,ISO_01_R1.fastq.gz,ISO_01_R2.fastq.gz,,,illumina,isolate,dna
LANES_01,LANES_01_L001_R1_001.fastq.gz,LANES_01_L001_R2_001.fastq.gz,,,illumina,,
LANES_01,LANES_01_L002_R1_001.fastq.gz,LANES_01_L002_R2_001.fastq.gz,,,illumina,,
MTX_01,MTX_01_R1.fastq.gz,MTX_01_R2.fastq.gz,,,illumina,metagenome,rna
ONT_01,ONT_01.fastq.gz,,,,ont,,
ONT_02,,,ONT_02.bam,,ont,,
ONT_03,,,,fastq_pass/barcode03,ont,,
UNKNOWN_01,UNKNOWN_01.fastq.gz,,,,,,
```

**Several files for one sample** (lanes, a top-up run, a nanopore barcode
folder): repeat the `sample` on each row, or give a folder in `fastq_dir`.
The rows must agree on `platform`, `sample_type` and `molecule`, and be all
paired or all single; a file may appear only once; a BAM must be a sample's
only row. Otherwise the run stops before reading anything, naming the
sample. Every file is checked on its own, and **one damaged file fails the
whole sample**, naming the file, rather than silently dropping a lane. The
platform is detected per lane and the lanes must agree (the sample takes
the weakest lane's confidence; lanes that disagree make it uncertain).
Then the files are joined (nf-core `cat/fastq`) and QC, the contamination
screen and the sample type judge the whole sample.

| Column | Required | Notes |
|---|---|---|
| `sample` | yes | Sample ID. Repeat it on several rows for a sample in several files. |
| `fastq_1` | yes, unless `bam` or `fastq_dir` | FASTQ: `.fastq`/`.fq`, plain or compressed as `.gz`, `.bz2`, `.xz` or `.zst`. Anything that is not `.fastq.gz` is rewritten as `.fastq.gz` in the work folder first; a file that will not decompress becomes a FAIL row. FASTA reads (`.fa`, `.fasta`, `.fna`) are accepted too: they have no quality values, so a placeholder is written, quality is not judged and the sample cannot PASS. Raw signal (POD5/FAST5), BCL run folders and `.sra` files are not reads yet: the run stops and says how to convert them. |
| `fastq_2` | no | Second mate for paired-end data. Leave blank for single-end or ONT. |
| `bam` | instead of FASTQ | A BAM, SAM or CRAM file, aligned or unaligned (CRAM needs the reference it was written against: `--cram_fasta`, or samtools fetches it from EBI's reference registry; if neither works, the sample is a FAIL row saying so): dorado and PacBio deliver reads as unaligned BAM. Give either FASTQ or one BAM per row. The BAM is checked from header to end, then converted to FASTQ (`samtools collate \| samtools fastq`); paired or single is decided from the reads' own READ1/READ2 flags, and secondary and supplementary alignments are left out, so each read counts once. The header is read too: software that only runs on one platform's data (`@PG` ccs or lima for PacBio, dorado or guppy for nanopore) counts as evidence, like read names; `@RG PL`, which people type, stands in for a blank `platform`, and a samplesheet that disagrees with it is flagged for review. |
| `fastq_dir` | instead of FASTQ files | A folder of one sample's FASTQ files (every `.fastq`/`.fq` in it, plain or compressed), each a file of single-end or long reads: for example a nanopore `fastq_pass/barcode03`. A local or cluster path (a web address cannot be listed as a folder). |
| `platform` | no | `illumina`, `dnbseq` (also `mgi`, `bgiseq`), `element` (also `aviti`), `ont` or `pacbio`. Leave blank to have it detected. MGI/DNBSEQ and Element get the same short-read QC as Illumina; PacBio is QC'd only when the reads are HiFi. |
| `sample_type` | no | `isolate`, `metagenome`, `amplicon` or `unknown`. Amplicon samples are not judged on read N50 or duplication, which amplicon libraries fail by design. |
| `target_taxid` | no | The NCBI taxonomy ID of the organism the sample is meant to be: `9606` human, `4932` *S. cerevisiae*, `10090` mouse, `562` *E. coli*. A number, not a name, because names change between database versions. With the contamination screen on, it sets what counts as off-target (anything that is not the target), reports `target_fraction`, and flags a target that is not the most common species (`TARGET_NOT_DOMINANT`, at most WARN) or that no read was classified to (`TARGET_NOT_SEEN`: contamination cannot be judged, at most WARN; Standard-8 holds no eukaryote but human, so use e.g. PlusPF for fungi). With `9606`, human reads are the sample, not contamination. With no target declared, a sample whose classified reads are at least half human is UNCERTAIN (`HOST_DOMINANT`): a human sample, or one swamped by human DNA, which only you can say. |
| `genome_size` | no | Expected genome size in bases (`5200000`, `5.2m`, `3g`), for the depth-of-coverage estimate. Without it, the NCBI expected size of the declared `target_taxid`, or of the most common species, is used (`assets/species_genome_size.tsv`, NCBI `species_genome_size.txt`). Depth is information only, reported for one organism's DNA; it feeds the next-pipeline advice (under 30x: an assembly will be fragmented). |
| `group` | no | Samples of one isolate or biological sample sequenced separately share a group: short and long reads of one isolate (the advice becomes one nf-core/bacass **hybrid** assembly, which closes plasmids), or the RNA and DNA of one metagenome (activity vs abundance becomes ready). Each sample is still QC'd on its own. |
| `molecule` | no | `dna` or `rna`. Use `rna` for RNA-seq, metatranscriptomes and direct RNA, even though Illumina RNA reads are cDNA. RNA samples are not judged on duplication, and FastQC's base-content and overrepresented-sequence modules are marked as known false alarms. RNA mostly cannot be detected from the reads, so it must be declared; direct-RNA reads written with U are flagged if `molecule` says `dna`. |
| `library` | no | `bisulfite` for a methylation library (`wgbs`, `rrbs`, `pbat`, `em-seq` and `methylation` are read as it). The reads show a bisulfite library only when its conversion is visible in them; in validation that was 19 of 23 ENA bisulfite runs (D93), so declare it when you know. A declared library is treated like a detected one: contamination is not judged and the next pipeline is nf-core/methylseq. If the reads show no conversion, the report says so. |

**Library structure is read from the reads** (D78), whatever was declared:

| Seen in the reads | Flag | What happens |
|---|---|---|
| Every read 30 bp or shorter | `INDEX_READS` | Index or barcode reads, not the sample's DNA: no QC, UNCERTAIN, give R1/R2 instead. |
| R1 30 bp or shorter, R2 longer | `BARCODE_READ` | A 10x Genomics-style single-cell library (cell barcode + UMI in R1): QC on R2; next pipeline nf-core/scrnaseq. |
| fastq_2 30 bp or shorter, fastq_1 longer | `INDEX_READ_AS_MATE` | An index (I1) read given as R2: QC on fastq_1 only, at most WARN. |
| One strand with little C (or G), and the converted base (T, or A) the commonest: under 10% and under half of the other base | `BISULFITE` | A bisulfite or EM-seq methylation library whose conversion shows in the reads: contamination not judged; next pipeline nf-core/methylseq. Some bisulfite runs show no conversion, so declare `library` if you know it. |
| One file whose consecutive reads share a name | `INTERLEAVED` | Both mates in one file: split into R1 and R2 before QC (never published). |

Plate-based single-cell (Smart-seq) has no barcode read, so it looks like
bulk RNA-seq; RRBS and EM-seq cannot be told from WGBS by the reads.

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
| `recommendations/analyses.md` | **What can be done with each sample** (D82): every analysis the data supports (assembly, typing, AMR, plasmids, comparative genomics, profiling, binning, resistome, expression, single-cell, methylation, variant calling), each marked ready, possible (with the caveat) or not now (with what is missing: depth, reads, replicates, more samples, a declaration), and the nf-core pipeline where one exists. Also as `analyses.tsv`. |
| `recommendations/` | The next nf-core pipeline for each ready sample (`recommendations.tsv`: a pipeline and pinned version, or the reason there is none), and one folder per rule (e.g. `bacass_short/`) with a `samplesheet.csv` in that pipeline's own format, a `params.yaml` and a `README.md` with the command. For the batch as a whole, `pathogensurveillance/` holds one nf-core/pathogensurveillance samplesheet of every isolate ready for assembly (species ID, reference choice and phylogeny across the batch). The pipelines take gzipped FASTQ, so samples given in another format get the reason instead of a samplesheet. Values only you know (primers, databases, reference) are never guessed: the pipeline will not run until you add them |
| `pipeline_info/` | Software versions actually used, run timeline, resource report, workflow diagram (`pipeline_dag.html`) |

Verdicts:

- **PASS:** the reads confirm the metadata and every QC number is within its limits.
- **WARN:** usable, but something needs a look; the report says what.
- **FAIL:** too little or too poor data to use.
- **UNCERTAIN:** something only you can settle: the platform could not be
  established, the reads contradict the metadata, only index reads were
  given, or a sample is mostly human with no `target_taxid`. The reason
  says which; no verdict is guessed.

## Testing

The `test` profile runs the whole pipeline on 21 small public samples,
pinned by commit, each with a known correct answer. From
[nf-core/test-datasets](https://github.com/nf-core/test-datasets): simulated
and real Illumina and nanopore reads, Illumina reads mislabelled as ONT, a
sample with no metadata, a malformed FASTQ, aligned BAM, unaligned PacBio
HiFi BAM, plain FASTQ, CRAM, two lanes of one sample, human DNA and RNA-seq
reads, and real 10x Genomics single-cell reads. From
[seq-triage-test-data](https://github.com/MacaulayWinter/seq-triage-test-data):
real NovaSeq and R10.4.1 nanopore reads of one *E. coli* isolate (linked as
one `group`, so the hybrid-assembly advice runs), and real MGI/DNBSEQ and
Element isolates. The test Kraken2 database and rRNA reference are tiny
stand-ins: they test the plumbing, not accuracy (that is what
[Validation](#validation) is for).

```bash
nextflow run main.nf -profile test,docker --outdir results         # Linux
nextflow run main.nf -profile test,docker,arm64 --outdir results   # Apple Silicon
Rscript tests/e2e/check_test_profile.R results
```

The check confirms that every sample got a verdict, that each verdict matches
`tests/e2e/expected_verdicts.tsv`, that no reads were copied into `results/`,
that every tool's version was recorded, that every sample sent to QC has
judged metrics, that MultiQC read every tool's output, that every sample has
a next-pipeline recommendation or a reason, and that every sample has an
analysis menu.

Unit tests, each on synthetic inputs built inside the test (no data needed):

```bash
Rscript tests/screen_metadata/test_screen_metadata.R     # metadata screen
Rscript tests/detect_platform/test_detect_platform.R     # platform and library detection
Rscript tests/combine_detection/test_combine_detection.R # combining the lanes of a sample
Rscript tests/triage_report/test_triage_report.R         # QC limits, verdicts, type, depth
Rscript tests/recommend/test_recommend.R                 # next pipelines and the analysis menu
bash tests/normalise/test_normalise.sh                   # format conversion (in the samtools image)
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

Tests check that the code does what it was designed to do. Validation checks
that the design is right, on public data it was never tuned on, against
answers known independently of the pipeline. Every round's pass marks were
fixed (and committed) before its data was drawn, and a failed round stays
recorded as failed: a rule chosen after seeing the data is only trusted once
a later round tests it. Full write-up: [docs/validation.md](docs/validation.md).

| Round | What it tests | Data | Result |
|---|---|---|---|
| Round 1 | Platform calls (A), every reported number (B), damaged files (C) | 40 ENA runs, 29 damaged samples | A 0 wrong at high, B 526/526; **C 28/29** (placeholder quality could PASS; fixed by D40, re-run 29/29) |
| Round 2 | as round 1, whole-ENA draw | 115 runs, 31 damaged samples | B 1158/1158, C 30/30; **A failed: 2 wrong at high confidence** (one an ENA labelling error; one led to rule D58) |
| Round 3 | D58 on unseen studies | 54 runs | **all pass marks met**; D58 stopped a wrong call it was not designed on |
| Layer D, rounds 1-2 | Contamination estimates on mixtures with known answers | 13 mixtures x 2 platforms, 2 databases | host DNA right in every case; **all-samples pass mark not met** (best: 22/24 with the full database) |
| Layer E | Sample type: isolate / metagenome / amplicon | 57 runs | **failed: 4 wrong at high confidence** (fragmented viral amplicons; a mostly-unclassified soil metagenome); led to D67 |
| Layer E round 2 | D67 on unseen studies | 49 runs | **failed: 7 wrong at high**, all true isolates called metagenome or amplicon (D67's own target held: nothing wrongly called isolate); led to D73 |
| Layer E round 3 | D73 on unseen studies | 50 runs | **failed: 4 wrong at high** (3 archaea, 1 phage); D73 prevented 6 others, and cost 6 correct metagenomes their high confidence |
| Layer E round 4 | D86, plus library structure (D78), rRNA (D81) and RNA starts (D85), on unseen studies | 74 runs | pre-registered |

Along the way validation found: a silent failure in a dependency (seqkit
reports success on a truncated gzip), two runs whose ENA instrument labels
are wrong, a full-size Kraken2 database that files one E. coli lineage
under *K. pneumoniae*, and a taxonomy change (NCBI's 2024 virus
reorganisation) that broke a rule on real data while synthetic tests passed.

## Known limitations

- **Platforms:** PacBio CLR, Ion Torrent, Ultima and older platforms are not
  QC'd. Without read names, an MGI or Element run can only be called a
  short-read platform at medium confidence.
- **Inputs:** raw signal (FAST5/POD5), BCL run folders and SRA files are
  named and refused with the conversion to do first. Index (I1/I2) files
  are recognised but not used. The next-pipeline samplesheet cannot be
  written for a `fastq_dir` folder, for nf-core/ampliseq with several files
  per sample, or for interleaved or FASTA reads (the advice says what to do
  first).
- **Sample kinds:** human and other eukaryotic samples need `target_taxid`
  (D76); a eukaryote missing from the Kraken2 database cannot be checked for
  contamination. Plate-based single-cell (Smart-seq), Hi-C and ChIP/ATAC
  libraries are not recognised; RRBS and EM-seq cannot be told from WGBS.
  Long-read amplicons cannot be recognised from the reads.
- **Sample type:** an organism the Kraken2 database knows poorly (archaea,
  novel bacteria, some phages) can look like a metagenome: its reads scatter
  across many species and genera. Since D86 a metagenome call is high only
  when the reads themselves (k-mer recurrence and GC spread, no database)
  also look like many genomes; on held-out data this removed 3 of 4 such
  errors, but it is not yet confirmed by a pre-registered round. Declare
  `sample_type` when you know it.
- **Contamination:** species level only (another strain of the same species
  is invisible); same-genus contamination is under-estimated on short reads;
  estimates depend on the Kraken2 database (see layer D).
- **Thresholds** are provisional defaults, each labelled as such in
  `assets/metric_guide.tsv`.

## How this was built

seq-triage was built by Macaulay Winter with
[Claude Code](https://claude.com/claude-code) as a coding assistant; commits
it wrote carry a `Co-Authored-By: Claude` line. The design decisions and
their reasons are recorded as they were made in [CLAUDE.md](CLAUDE.md)
(D1-D87), including every rule that validation overturned and every CI
failure and how it was fixed.

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
- [x] PacBio, MGI/DNBSEQ and Element recognised and named
- [x] Validation round 2 (115 runs, whole-ENA draw) and round 3 (54 runs, testing the rule round 2 led to)

**v0.2.0**

- [x] BAM, SAM and CRAM input; plain, `.bz2`, `.xz`, `.zst` and FASTA reads (D55, D69, D79)
- [x] Several files per sample: lanes and nanopore barcode folders (D75)
- [x] QC for MGI/DNBSEQ, Element and PacBio HiFi (D64, D69)
- [x] Library structure from the reads: index, 10x single-cell, bisulfite, interleaved (D78)
- [x] Isolate / metagenome / amplicon from the reads, validated in layer E rounds 1-3 (D65-D74, D84)
- [x] Contamination screen, validated as layer D (D41); human and other targets (D76)
- [x] rRNA share of RNA samples (D81); genome size and depth of coverage (D80)
- [x] Next nf-core pipeline with a ready-to-fill samplesheet (D66), and the analysis menu (D82)
- [x] Linked samples: hybrid assembly, RNA/DNA pairs (D83)

**Next**

- A per-read GC-spread signal for isolates the database knows poorly (layer E rounds 2-3)
- QC for declared Ion Torrent samples
- A validation round for the library-structure, depth and rRNA rules (D78-D81)
- Strain-level contamination

## Data policy

No real research data is stored in this repository. Test data comes from
public sources and is referenced by URL. `.gitignore` blocks sequencing file
formats anywhere in the tree as a second safeguard.

**Human reads are personal data.** seq-triage never publishes reads: its
results hold counts, QC numbers and species names only. The work folder
does hold copies of the reads while it runs, so keep it, like the input,
where your data-protection rules allow, and delete it after the run.

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

`main` is protected: GitHub refuses a merge until every CI check has passed
(a check skipped because a pull request changes only documentation counts as
passed). To wait for the checks and merge in one step:

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
