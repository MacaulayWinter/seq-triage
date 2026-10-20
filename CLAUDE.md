# CLAUDE.md: seq-triage

Working notes for building this pipeline with Claude Code: the plan, the
decisions made (with reasons), and where things stand. Updated at the end of
every session.

## Project

A Nextflow DSL2 pipeline that inspects incoming sequencing data, validates its
metadata, routes each sample to the right QC, and writes a plain-language
summary. A portfolio piece: every design choice must be defensible in
interview.

**Working agreement:** Macaulay runs things and makes the decisions. Claude
scaffolds, explains and debugs. Build step by step, never a finished repo in
one go.

## Hard rules

- No real research data in the repo, ever. Public test data only
  (nf-core/test-datasets, or small subsampled public ENA runs), referenced
  by URL.
- Reuse nf-core modules for every standard tool. Our own code is the
  detection, routing and reporting logic only.
- Detection reports a confidence level, and says `uncertain` rather than
  guessing.
- Pin every container and tool version.
- Keep this file current at the end of each session.

## v0.1 scope

1. Git repo + DSL2 skeleton (main.nf, nextflow.config, modules/, conf/,
   assets/, README).
2. Samplesheet `sample, fastq_1, fastq_2, platform?, sample_type?`,
   validated with nf-schema.
3. Platform and layout detection from `seqkit stats`, flagging any
   contradiction with the metadata.
4. Routed QC: Illumina -> fastp + FastQC; ONT -> NanoPlot; MultiQC over all.
5. Summary report: per-sample pass/warn/fail with reasons, detected vs
   declared.
6. `-profile test` on tiny public data + a GitHub Actions workflow running it.

**v0.2 (not yet):** isolate/metagenome/amplicon classification, host and
contaminant screen, and a recommendation step naming the downstream nf-core
pipeline with a ready-to-run params file.

## Architecture (agreed)

```
samplesheet ─► nf-schema ─► ch_samples [meta, reads]
     │                         ├─► SEQKIT_STATS ─► DETECT_PLATFORM ─┐
     └─► SCREEN_METADATA (1 task, whole sheet) ─► per-sample flags ─┤
                               └──── join by meta.id ◄──────────────┘
                                          ▼
              branch: empty | illumina | ont | unrouted
              illumina: FASTP ─► FASTQC     ont: NANOPLOT
                                          ▼
                     MULTIQC  +  TRIAGE_REPORT (one task, all samples)
```

- `[meta, reads]` tuples throughout. Every per-sample fact lives in `meta`.
- Detection is a process, because it needs seqkit output, and it runs per
  sample so it parallelises and resumes per sample.
- The join on `meta.id` is what makes `branch` possible: routing can only
  use values carried in the channel.
- An `empty` branch keeps zero-read files from crashing fastp/NanoPlot. They
  become FAIL lines in the report instead.
- Resource labels: nf-core standard (`process_single` ... `process_high`).
  Local modules use `process_single`. `resourceLimits` caps everything at the
  CI runner size.
- Errors: exit codes 104 and 130-145 (killed / out of memory) retry once with
  more resources. Anything else stops the run (`finish`). `uncertain` is a
  result, never an error.

## Decisions log

| # | Decision | Why |
|---|----------|-----|
| D1 | **Conservative routing.** No sample can PASS unless the reads independently confirm its metadata (high-confidence detection that agrees with what was declared). | User rule: when detection is in doubt, interpret conservatively. |
| D2 | Routing table: high + match -> detected, PASS possible. High + missing -> detected, at most WARN ("inferred"). High + contradiction -> detected, at most WARN ("review"). Medium/low + declared -> declared, at most WARN ("unconfirmed"). Medium/low/uncertain + nothing declared -> not routed, UNCERTAIN. Empty -> FAIL. | Puts D1 into practice. QC metrics can only lower this ceiling. |
| D3 | "High" confidence needs several independent signals to agree (read-length spread, max length, Q30, mate structure). One signal is never enough. | Avoids overconfident calls, e.g. PacBio HiFi looks long like ONT but has high Q30. |
| D4 | Detection and report logic in **R** (`bin/*.R`), base R only, pinned `r-base` 4.5.x container. | User is strong in R and must defend it line by line. Base R keeps the container small. |
| D5 | Two-layer metadata validation: **nf-schema** for hard schema errors (stop before reading data), then **SCREEN_METADATA** in R over the whole sheet (warnings, carried into the report). | Mirrors the user's own pattern: validate early and stop loudly, then screen for inconsistencies. Collisions only show up across rows, so the R screen runs once on the whole sheet. |
| D6 | SCREEN_METADATA adapts helpers from the user's R_Bible: `blank_na` (empty-value tokens -> NA), `normalise_identifier` plus match-tier recording (e.g. `nanopore` -> `ont`, reported as "normalised"), collision detection (`S1` vs `s1`), filename-noise detection (`_S1_L001_R1_001`), basename-vs-ID checks, mate checks (R1 and R2 must be different files, and read counts must match exactly once seqkit has run). | Reuses the user's established screening conventions. Generic helpers only, copied with a provenance comment. |
| D7 | Never silently rewrite IDs or values. Normalise for matching, but always report the original and how it was matched. | Deliberately unlike the old `sub("-", "", ...)` habit. Silent edits can't be audited. |
| D8 | No R_Bible *project* scripts are copied into the repo. | They contain real sample IDs, cohort codes and HPC paths, which would break the "no research data" rule. |
| D9 | Apple Silicon: an `arm` profile with `docker.runOptions = '--platform=linux/amd64'`, plus Docker Desktop Rosetta. CI runs natively on amd64. | Most BioContainers images are amd64-only. Explicit and documented beats implicit. |
| D10 | Commits use GitHub noreply email `108404294+MacaulayWinter@users.noreply.github.com`. | User does not want a private email in public history. |
| D11 | Commits made by Claude carry a `Co-Authored-By: Claude` line. Commits the user types do not. | Kept on the first commit (user pushed it as-is). Revisit if wanted. |
| D12 | The report explains every metric it shows, driven by `assets/metric_guide.tsv` (metric, applies_to, what it measures, how to read it, warn/fail thresholds, caveat, provenance). The report quotes value + threshold + meaning per sample, and ends with a "How to read these numbers" section covering only the metrics present. | User wants interpretation advice in the output. Keeping it as data makes thresholds reviewable, overridable (`--metric_guide`) and traceable. Wording adapted from R_Bible Parts 16.2 and 18.3. |
| D13 | Metrics are tagged by what they were measured on (`illumina_reads`, `ont_reads`, later `assembly`). R_Bible's N50 guidance (> 100 kb good) is for **assembly contigs** and is not applied to **read N50**. | Same formula, different scale. Applying assembly thresholds to reads would flag every healthy ONT run. Assembly guidance waits for v0.2. |
| D14 | **Reproducibility in four layers** (see `docs/reproducibility.md`): (1) launch env: `environment.yml` with exact pins plus a platform lockfile `envs/*.lock.txt` (every package, URL + md5); (2) `manifest.nextflowVersion` + `NXF_VER` in CI; (3) one pinned container per process, nf-core modules frozen in `modules.json`; (4) run-time `software_versions.yml` collected from every process. Analysis tools live only in containers, never in the conda env. | User requirement: every version recorded, and the environment installable by anyone. One source of truth per layer prevents drift. `environment.yml` verified to solve on linux-64, linux-aarch64, osx-64 and osx-arm64 (2026-10-20). |
| D15 | **Provisional thresholds drafted by Claude**, every row labelled `provisional pipeline default` (or "adapted from R_Bible ..."). Read-count thresholds are *sanity floors*, not adequacy thresholds. Depth-type metrics (total bases, GC) carry no threshold, because v0.1 has no genome size or species. Duplication and adapter content can WARN but never FAIL. | User chose option 2. Conservative: only flag what can be judged without knowing the organism. |

## Open questions

- Final allowed values for `sample_type` (proposed:
  `isolate | metagenome | amplicon | unknown`). Recorded and checked only in
  v0.1.
- Exact test-data URLs: confirm they exist in nf-core/test-datasets before
  use.
- Review the provisional thresholds in `assets/metric_guide.tsv` and replace
  them with lab-practice values over time (change `provenance` when you do).
- Verify on the first real run: the fastp JSON and NanoStats field names in
  `source_field`, and the claim that Q7 is below every MinKNOW default pass
  filter.
- Tiny test data will trip the read-count sanity floors, so CI reports will
  show WARN/FAIL. That is correct behaviour; decide whether to say so in the
  README or to choose larger test files.
- Add a linux-64 lockfile, generated from a real install in CI, so Linux users
  can do a byte-identical rebuild too.
- Depth of coverage needs a genome size, so it is not reported in v0.1. The
  guide states this explicitly.

## Environment (recorded 2026-10-20)

- Mac: Apple M1 Pro (arm64), macOS 26 (Darwin 25.1.0)
- Docker Desktop 29.2.1, Rosetta amd64 emulation confirmed working
- conda: Miniforge 25.11.0. Channels: conda-forge, bioconda (`defaults`
  removed)
- Dev env `seq-triage-dev`: OpenJDK 17.0.18, Nextflow 26.10.2,
  nf-core tools 4.1.0, nf-test 0.9.5
- Local R 4.5.2. git 2.50.1 (Apple) / 2.56.0 (in env). gh 2.102.0
- GitHub: `MacaulayWinter/seq-triage` (private), HTTPS auth via gh keyring

## Status

**Session 1 (2026-10-20)**
- Done: environment checked and installed, git identity set, repo
  initialised, `.gitignore` committed (`45a93d9`), private GitHub repo
  created and pushed, CLAUDE.md written. Reproducible launch env
  (`environment.yml` + osx-arm64 lockfile), `docs/reproducibility.md`, and
  provisional `assets/metric_guide.tsv` (13 metrics) drafted.
- Next: Nextflow skeleton (`nextflow.config`, `main.nf`, `conf/base.config`,
  README), then install the nf-core modules.
- Reminder: run `conda activate seq-triage-dev` before using
  `nextflow` or `nf-core`.
