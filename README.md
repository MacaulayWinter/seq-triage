# seq-triage

A Nextflow DSL2 pipeline that inspects incoming sequencing data before any
analysis starts. It validates the sample metadata, works out from the reads
themselves which platform produced them, routes each sample to the right QC,
and writes a plain-language summary that says which samples are ready and why.

> **Status: under development (v0.1).** The configuration, parameter
> validation and reproducible environment are in place. The analysis steps
> are being built now. See [Roadmap](#roadmap) for what works today.

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
| low or medium | present | declared | WARN (not confirmed by reads) |
| uncertain | missing | not routed | UNCERTAIN |
| empty file | any | not routed | FAIL |

QC metrics can lower a verdict but never raise it. Every metric in the report
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
version pinned. The analysis tools themselves run in version-pinned
containers and need no installation. See
[docs/reproducibility.md](docs/reproducibility.md) for how each layer is
pinned.

## Usage

```bash
nextflow run main.nf -profile docker --input samplesheet.csv --outdir results
```

**On Apple Silicon Macs**, add `emulate_amd64`. Most bioinformatics
containers are built for Intel/AMD (amd64) processors only, and this runs
them under Rosetta emulation:

```bash
nextflow run main.nf -profile docker,emulate_amd64 --input samplesheet.csv --outdir results
```

Enable Rosetta first: Docker Desktop → Settings → General → *Use Rosetta for
x86_64/amd64 emulation on Apple Silicon*. Profile order matters: put
`emulate_amd64` after `docker`.

### Samplesheet

A CSV with one row per sample:

```csv
sample,fastq_1,fastq_2,platform,sample_type
ISO_01,ISO_01_R1.fastq.gz,ISO_01_R2.fastq.gz,illumina,isolate
ONT_01,ONT_01.fastq.gz,,ont,
UNKNOWN_01,UNKNOWN_01.fastq.gz,,,
```

| Column | Required | Notes |
|---|---|---|
| `sample` | yes | Unique ID. |
| `fastq_1` | yes | Gzipped FASTQ. |
| `fastq_2` | no | Second mate for paired-end data. Leave blank for single-end or ONT. |
| `platform` | no | `illumina` or `ont`. Leave blank to have it detected. |
| `sample_type` | no | `isolate`, `metagenome`, `amplicon` or `unknown`. Recorded only in v0.1. |

Unknown parameters and unknown samplesheet columns are errors, not warnings.
This stops a typo such as `fastq2` from silently turning paired data into
single-end.

## Roadmap

**v0.1** (in progress)

- [x] Reproducible launch environment (`environment.yml` + lockfile)
- [x] Configuration, resource labels and retry policy
- [x] Parameter validation
- [ ] Samplesheet validation and metadata screening
- [ ] Platform and layout detection with confidence levels
- [ ] Routed QC (fastp + FastQC, NanoPlot) and MultiQC
- [ ] Triage summary report with metric interpretation
- [ ] `-profile test` on small public data, run in GitHub Actions

**v0.2** (planned)

- Isolate vs metagenome vs amplicon classification
- Host and contaminant screen
- Recommendation of the right downstream nf-core pipeline, with a
  ready-to-run parameters file

## Data policy

No real research data is stored in this repository. Test data comes from
public sources and is referenced by URL. `.gitignore` blocks sequencing file
formats anywhere in the tree as a second safeguard.

## Citations

If you use seq-triage, please cite the tools it is built on:

- **Nextflow**: Di Tommaso P, et al. (2017) Nextflow enables reproducible
  computational workflows. *Nature Biotechnology*.
  [doi:10.1038/nbt.3820](https://doi.org/10.1038/nbt.3820)
- **nf-core**: Ewels PA, et al. (2020) The nf-core framework for
  community-curated bioinformatics pipelines. *Nature Biotechnology*.
  [doi:10.1038/s41587-020-0439-x](https://doi.org/10.1038/s41587-020-0439-x)
- **SeqKit**: Shen W, et al. (2016) SeqKit: a cross-platform and ultrafast
  toolkit for FASTA/Q file manipulation. *PLOS ONE*.
  [doi:10.1371/journal.pone.0163962](https://doi.org/10.1371/journal.pone.0163962)
- **fastp**: Chen S, et al. (2018) fastp: an ultra-fast all-in-one FASTQ
  preprocessor. *Bioinformatics*.
  [doi:10.1093/bioinformatics/bty560](https://doi.org/10.1093/bioinformatics/bty560)
- **FastQC**: Andrews S. (2010) FastQC: a quality control tool for high
  throughput sequence data.
  [bioinformatics.babraham.ac.uk](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- **NanoPlot**: De Coster W, Rademakers R. (2023) NanoPack2:
  population-scale evaluation of long-read sequencing data. *Bioinformatics*.
  [doi:10.1093/bioinformatics/btad311](https://doi.org/10.1093/bioinformatics/btad311)
- **MultiQC**: Ewels P, et al. (2016) MultiQC: summarize analysis results for
  multiple tools and samples in a single report. *Bioinformatics*.
  [doi:10.1093/bioinformatics/btw354](https://doi.org/10.1093/bioinformatics/btw354)

## Author

Macaulay Winter ([@MacaulayWinter](https://github.com/MacaulayWinter))

## Licence

[MIT](LICENSE)
