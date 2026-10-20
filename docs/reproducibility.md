# Reproducibility

seq-triage pins software at four layers. Each layer answers a different
question, and the gaps between them are where pipelines usually stop being
reproducible.

| Layer | What it pins | Where it is recorded | Who uses it |
|---|---|---|---|
| 1. Launch environment | Java, Nextflow, nf-core tools, nf-test | `environment.yml` (exact versions, any OS) and `envs/*.lock.txt` (exact builds plus checksums, one platform) | Anyone running or developing the pipeline |
| 2. Nextflow itself | The workflow engine version | `manifest.nextflowVersion` in `nextflow.config`. CI sets `NXF_VER` | Nextflow refuses to run on an incompatible version |
| 3. Analysis tools | seqkit, fastp, FastQC, NanoPlot, MultiQC, R | One container per process, pinned by exact tag. nf-core module versions are frozen by git commit in `modules.json` | Every pipeline run |
| 4. What actually ran | The versions reported by each tool at run time | `results/pipeline_info/software_versions.tsv` | The record you cite in a methods section |

## Layer 1: the launch environment

```bash
# Portable: same versions on Linux, Intel Mac and Apple Silicon
conda env create -f environment.yml
conda activate seq-triage-dev

# Byte-identical rebuild on Apple Silicon (every dependency, with md5 checksums)
conda create -n seq-triage-dev --file envs/seq-triage-dev.osx-arm64.lock.txt
```

`environment.yml` pins the four packages we depend on directly. Their own
dependencies (Python libraries and so on) are chosen by the solver, so a
rebuild months later can differ in those. The lockfile closes that gap: it
lists every one of the 166 packages by exact URL and checksum, so conda
installs exactly those files and does no solving at all. Lockfiles are
platform-specific, which is why both files exist.

`environment.yml` was confirmed to solve on linux-64, linux-aarch64, osx-64
and osx-arm64 with identical versions of the pinned packages (2026-10-20).

**Why the analysis tools are not in this environment:** they run inside
containers (layer 3). Putting them here as well would leave two sources of
truth for their versions, and they could drift apart.

## Layer 3: containers

Each nf-core module declares its own container, pinned to an exact version
tag, never `latest`. Most current nf-core modules use Seqera community images
(e.g. `community.wave.seqera.io/library/fastp:1.4.0--df415eab5e33cb11`); NanoPlot
still uses BioContainers (`quay.io/biocontainers/nanoplot:1.47.0--pyhdfd78af_0`).
`nf-core modules install` records each module's git commit in `modules.json`,
so the module code and its container version are frozen together.

| Tool | Version | nf-core module commit |
|---|---|---|
| seqkit | 2.13.0 | `8fff8cfa10` |
| fastp | 1.4.0 | `1e31555fd8` |
| FastQC | 0.13.0 | `18eb358065` |
| NanoPlot | 1.47.0 | `cc4d441624` |
| MultiQC | 1.35 | `98403d15b0` |

Local modules (our R scripts) use `docker.io/rocker/r-ver:4.5.2`, pinned by
**tag and digest** (`@sha256:fd4ccdd3...`). The tag names the R version; the
digest guarantees the exact image even if the tag is ever moved.

On Apple Silicon, amd64-only images run under Rosetta emulation via
`-profile docker,emulate_amd64` (see the README).

## Layer 4: the run-time record

Every module reports `[process, tool, version]` to a Nextflow `versions`
*topic channel* (the current nf-core pattern, replacing per-process
`versions.yml` files). The pipeline gathers them into
`results/pipeline_info/software_versions.tsv`. This checks layers 1-3: if a
container ever ships something other than its tag claims, this file shows
it. The end-to-end test fails if any tool is missing from it.

## Test data

`assets/samplesheet_test.csv` points at nf-core/test-datasets **by commit**
(`76179c1c...`), not by branch name, so the test inputs cannot change
underneath the pipeline.

## Recording a change

Any version bump is a deliberate commit that updates `environment.yml`,
regenerates the lockfile, and is noted in `CLAUDE.md`:

```bash
conda list -n seq-triage-dev --explicit --md5 > envs/seq-triage-dev.osx-arm64.lock.txt
```
