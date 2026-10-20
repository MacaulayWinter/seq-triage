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
| D9 | Apple Silicon: an `emulate_amd64` profile (same name as the nf-core template's), used as `-profile docker,emulate_amd64`, with `docker.runOptions = '--platform=linux/amd64'`, plus Docker Desktop Rosetta. CI runs natively on amd64. | Most BioContainers images are amd64-only. Explicit and documented beats implicit. |
| D10 | Commits use GitHub noreply email `108404294+MacaulayWinter@users.noreply.github.com`. | User does not want a private email in public history. |
| D11 | Commits made by Claude carry a `Co-Authored-By: Claude` line. Commits the user types do not. | Kept on the first commit (user pushed it as-is). Revisit if wanted. |
| D12 | The report explains every metric it shows, driven by `assets/metric_guide.tsv` (metric, applies_to, what it measures, how to read it, warn/fail thresholds, caveat, provenance). The report quotes value + threshold + meaning per sample, and ends with a "How to read these numbers" section covering only the metrics present. | User wants interpretation advice in the output. Keeping it as data makes thresholds reviewable, overridable (`--metric_guide`) and traceable. Wording adapted from R_Bible Parts 16.2 and 18.3. |
| D13 | Metrics are tagged by what they were measured on (`illumina_reads`, `ont_reads`, later `assembly`). R_Bible's N50 guidance (> 100 kb good) is for **assembly contigs** and is not applied to **read N50**. | Same formula, different scale. Applying assembly thresholds to reads would flag every healthy ONT run. Assembly guidance waits for v0.2. |
| D14 | **Reproducibility in four layers** (see `docs/reproducibility.md`): (1) launch env: `environment.yml` with exact pins plus a platform lockfile `envs/*.lock.txt` (every package, URL + md5); (2) `manifest.nextflowVersion` + `NXF_VER` in CI; (3) one pinned container per process, nf-core modules frozen in `modules.json`; (4) run-time `pipeline_info/software_versions.tsv`, collected from every process via the `versions` topic channel. Analysis tools live only in containers, never in the conda env. | User requirement: every version recorded, and the environment installable by anyone. One source of truth per layer prevents drift. `environment.yml` verified to solve on linux-64, linux-aarch64, osx-64 and osx-arm64 (2026-10-20). |
| D15 | **Provisional thresholds drafted by Claude**, every row labelled `provisional pipeline default` (or "adapted from R_Bible ..."). Read-count thresholds are *sanity floors*, not adequacy thresholds. Depth-type metrics (total bases, GC) carry no threshold, because v0.1 has no genome size or species. Duplication and adapter content can WARN but never FAIL. | User chose option 2. Conservative: only flag what can be judged without knowing the organism. |
| D16 | Pin **nf-schema 2.8.0**, not 3.0.0. | 3.0.0 (2026-09-29) is a new major release with breaking changes to `paramsSummaryLog()`/`paramsHelp()`, and the nf-core 4.1 template still pins 2.5.1. The functions we rely on (`validateParameters`, `samplesheetToList`) are unchanged in 3.0, so upgrading later is a small change. |
| D17 | **Resource labels use nf-core names but laptop-sized values**: single 1 CPU/2 GB, low 2/4 GB, medium 4/8 GB, high 8/16 GB. Memory and time double on retry; CPUs do not. | Every tool here is lightweight QC. nf-core template values are HPC-sized (process_high = 72 GB). Containers on this Mac are capped by the Docker VM (11.7 GB), and private-repo GitHub runners are about 2 CPU / 7 GB. A kill for lack of memory is fixed by memory, not threads. |
| D18 | **Retry only exit 104 and 130-145, once; everything else `finish`.** Fusion Snapshot codes 175-177 from the nf-core template are deliberately left out. Verified 2026-10-20 with a scratch test: exit 137 retried with 4 GB -> 8 GB and succeeded; exit 1 got one attempt and the run stopped. | Retrying a genuine error only hides it. 175-177 cannot occur without Seqera Fusion. **Known interaction:** after a genuine error, `finish` submits no new tasks, *including pending retries elsewhere*. Accepted: the fix is to correct the error and `-resume`, which reruns the killed task anyway. |
| D19 | `validation.logging.unrecognisedParams` and `unrecognisedHeaders` set to `error` (nf-schema default is `warn`). Tested: `--outdr` typo -> exit 1; valid run -> exit 0. | Conservative rule: a misspelled samplesheet column (`fastq2`) would otherwise be dropped with a warning, silently turning paired data into single-end. The same options exist in nf-schema 3.0, so this survives an upgrade. |
| D20 | **MIT licence** (user choice, 2026-10-20), "for now". | Same licence as nf-core, so it matches the ecosystem the pipeline builds on. Without any licence the code is all rights reserved, even once public. |
| D21 | **Validation boundary:** anything unambiguously wrong in one row and checkable without reading the data is a hard error at parse time (`assets/schema_input.json`, plus Groovy for R1 == R2). Cross-row or merely suspicious issues (`S1` vs `s1`, filename-like IDs, the same FASTQ in two rows) are R-screen warnings. Tested 2026-10-20: 1 valid + 8 broken sheets (duplicate ID, `/` in ID, missing file, not gzipped, `PacBio`, R1 == R2, `fastq2` header, missing fastq_1); all 8 stop with a message naming row, field and value. | Fail fast on certainties, report on suspicions. R1 == R2 must be a hard error, because seqkit would see equal counts and call it paired with high confidence. JSON Schema cannot compare fields, hence the Groovy check. A shared FASTQ can be legitimate (our test sheet reuses files on purpose), so it is a warning, not an error. |
| D22 | `platform`/`sample_type` accept case variants and empty tokens (`N/A`, `.`, `none` ...) via a case-insensitive `(?i)` pattern. This was confirmed to work in nf-schema 2.8.0. Blank cells (returned by nf-schema as `[]`) become `null` in `main.nf`. Non-blank values reach the R screen **exactly as written**. | One meaning of "not declared" downstream. Normalisation happens in one place (the R screen), which reports it (D7). |
| D23 | **`bin/screen_metadata.R`** (base R) reads the RAW samplesheet and writes one row per sample: raw + normalised platform/type with a match tier (`blank`/`empty_token`/`exact`/`normalised`), plus warning flags `EMPTY_TOKEN`, `NORMALISED`, `WHITESPACE`, `ID_COLLISION`, `ID_LOOKS_LIKE_FILENAME`, `BASENAME_MISMATCH`, `MATES_POSSIBLY_SWAPPED`, `FASTQ_SHARED`, `LAYOUT_CONFLICT`. One root cause gives one flag (file noise is stripped from both ID and file name before comparing). | Raw input, because the screen must see what the user typed. Helpers adapted from R_Bible Parts 01 and 23. **Finding:** nf-schema silently trims whitespace (tested), so the `WHITESPACE` flag is the only record of that edit, and the output `sample` column is trimmed to match `meta.id`. Otherwise a padded ID would silently drop out of the join. |
| D24 | **Self-checking tests in the repo:** `tests/screen_metadata/` holds a synthetic 13-row samplesheet (no real files), the expected flag codes per row, and `test_screen_metadata.R`, which exits 1 on any mismatch. Verified it fails when an expectation is wrong. | Shows the logic is checked, not just written. Runnable in CI. Synthetic sheets contain no data, so the hard rule is kept. |
| D25 | **Test samplesheet (6 rows, ~17 MB, all public):** `BF_ILLUMINA_PE` (B. fragilis test2 R1+R2; illumina, isolate; expect PASS), `SC2_ONT` (SARS-CoV-2 test_2, real ONT; ont, amplicon; expect PASS), `SIM_ONT_Q40` (B. fragilis nanopore, placeholder Q40; ont; expect WARN, reads do not confirm), `CONTRADICTION` (B. fragilis R1 declared ont; expect WARN, review), `NO_METADATA` (SARS-CoV-2 test_2, nothing declared; expect WARN, inferred), `ILLUMINA_SE` (B. fragilis R1; illumina; expect PASS). | Every row has a known correct answer, so the end-to-end test can check verdicts automatically. Shared files also exercise `FASTQ_SHARED`. User approved 2026-10-20. |
| D26 | **Uninformative quality:** if Q20 % and Q30 % are both 100, quality is excluded from the detection evidence and flagged `QUALITY_UNINFORMATIVE`. Confidence comes from the remaining signals and is capped accordingly. | Placeholder or simulated quality (`SIM_ONT_Q40`: every base `I`) would otherwise count as evidence against ONT. Excluding it is the conservative, honest reading. User approved. |
| D27 | **`exempt_types` column in `metric_guide.tsv`:** a metric is skipped (and the report says why) when the sample's `sample_type` is listed. `read_n50` and `duplication_rate` are exempt for `amplicon`. This means `sample_type` now affects verdicts in v0.1 (README updated). | Short reads and high duplication are the design of amplicon libraries, not faults. Real SARS-CoV-2 ONT has a read N50 of 686 bp. User chose option (b) over a blanket warn-only rule. |
| D28 | **Walking skeleton first; nothing reaches `main` untested.** Each feature is built on its own branch. CI (`.github/workflows/ci.yml`) runs the R unit tests, installs from `environment.yml`, runs `-profile test` end to end, and checks the results with `tests/e2e/check_test_profile.R` (every row has a verdict, verdicts match `tests/e2e/expected_verdicts.tsv`, no FASTQ published, every tool versioned). Merge only when green. | User requirement (2026-10-20): test code end to end at the earliest version, not just publish it. The first end-to-end run immediately caught a bug the unit tests missed (false `MATES_POSSIBLY_SWAPPED` on single-file `test_2.fastq.gz`), now fixed with regression tests. |
| D29 | **`meta` carries only what downstream tasks need.** Screen flag counts and raw values stay in the screen TSV, which the report reads directly. | `meta` is part of every task's cache key. Tested: with the flag count in `meta`, a screen wording change re-ran NanoPlot. Without it, only SCREEN_METADATA and TRIAGE_REPORT re-run (20/22 cached). |
| D31 | **Own test data lives in a separate public repo, `MacaulayWinter/seq-triage-test-data`**, pinned by commit (`b18e39da`). It holds subsamples of public ENA runs only, with accession, cutting command and checksums in its README. First contents: *E. coli* ECE107 (PRJNA1287288), NovaSeq SRR36800568 (first 120,000 pairs, 19 MB) and R10.4.1 SRR34410547 (first 12,000 reads, 89 MB). Added to the test samplesheet as `ECE107_NOVASEQ` and `ECE107_R10`. | User chose option A (2026-10-20): CI must show a PASS under the real thresholds, not test-only ones. Hosted data is needed for `-profile test` to stay one command; a separate repo is the nf-core test-datasets model and keeps the "no data in this repo" rule. ECE107 was chosen over the bacass *Klebsiella* (no documented source) for traceable provenance, a current instrument, and the same isolate on both platforms. |
| D32 | **Keep the Illumina 100,000-read floor for now, counted as fastp counts it: reads, both mates included** (so 50,000 pairs for paired data). | User decisions, 2026-10-20. Matches the tool's own number, so the report quotes exactly what fastp says. Every real paired run in the survey clears it, including Portiera at 97,729 pairs (195,458 reads); only the 100-pair SARS-CoV-2 file fails. Depth of coverage is the better measure, but it needs a genome size (v0.2). |
| D30 | **Pin everything external by immutable reference:** test data by nf-core/test-datasets commit (`76179c1c`), R image by digest (`rocker/r-ver:4.5.2@sha256:fd4ccdd3...`), GitHub Actions by commit SHA (tag in a comment), Miniforge installer `26.7.2-0`, Nextflow via `NXF_VER`. Avoid releases younger than about two weeks (upload-artifact v7.0.1 chosen over v7.0.2, which was 13 days old). | Tags and branches can move; commits and digests cannot. SHA-pinned actions also close a known supply-chain risk. |

## Test data findings (2026-10-20, inspected with seqkit 2.13.0)

nf-core/test-datasets, `modules` branch, `data/genomics/...`:

| File | Reads | Length min/mean/max | N50 | Q30 % | AvgQual | What it really is |
|---|---|---|---|---|---|---|
| `prokaryotes/bacteroides_fragilis/illumina/fastq/test2_{1,2}` | 24,993 pairs | 151/151/151 | 151 | 96 / 90 | 27.5 / 23.7 | **Simulated** Illumina (headers `@CP009257.1_0_0/1`), NovaSeq-style binned quality (F, -, 8, #) |
| `prokaryotes/bacteroides_fragilis/nanopore/fastq/test` | 1,000 | 62/33,344/773,417 | 63,032 | **100** | **40.00** | **Simulated E. coli K-12** (header `NC-000913.3_..._unaligned_...`), *not* B. fragilis. **Every base quality is `I` (Q40): a placeholder.** |
| `sarscov2/nanopore/fastq/test_2` | 4,000 | 131/537/3,523 | 686 | 13 | 11.9 | **Real** ONT (MinKNOW UUID + `runid=` headers, 62 distinct quality values). Amplicon-length reads. |
| `sarscov2/nanopore/fastq/test` | 100 | 228/383/1,664 | 365 | 24 | 12.6 | Real ONT, very small |

Implications:
- Constant quality (Q20 = Q30 = 100%) must make quality an **excluded** signal (flagged), not evidence against ONT.
- Read headers can be supporting evidence when recognisable (MinKNOW `runid=`, Illumina `@instrument:run:flowcell:`), and are neutral otherwise. Test data headers are mostly simulated.
- Real amplicon ONT has read N50 well under 1 kb, so the provisional `read_n50` FAIL < 1 kb (D15) would wrongly fail it.
- seqkit/stats module container `community.wave.seqera.io/library/seqkit:2.13.0--3bacfb68370dd835` is **amd64 only** (digest `sha256:4f2a66d7...8567507`), so it needs `emulate_amd64` on the M1. A separate arm64 tag `...--2c3127bb2c908a0f` exists.
- nf-core modules now report versions through `topic: versions` channels, not `versions.yml`. `docs/reproducibility.md` layer 4 needs updating.

## Test data survey (2026-10-20, seqkit 2.13.0)

12 real public samples, kept outside the repo in `~/seqtriage-scratch/survey/` (manifest with source commit or accession, md5s, seqkit stats, screen output). Read counts per file; ONT quality is seqkit AvgQual.

| Sample (source) | Platform / instrument | Reads | N50 | Q30 % | AvgQual | Against D15 floors |
|---|---|---|---|---|---|---|
| ECE107 NovaSeq (ENA SRR36800568, first 120k pairs) | Illumina NovaSeq 6000, 4 binned quality values | 120,000 pairs | 151 | 87 / 88 | 23.1 | **PASS** |
| A1403KPN (bacass) | Illumina NextSeq 500 | 643,986 pairs | 150 | 96 / 95 | 30.3 | **PASS** |
| ERR044595 (bacass) | Illumina HiSeq 2000, 75 bp, old `#idx/1` headers | 1,000,000 pairs | 75 | 96 / 95 | 28.3 | **PASS** |
| ERR5766176 (taxprofiler) | Illumina MiSeq | 632,060 pairs | 151 | 91 / 85 | 24.6 | **PASS** |
| SRR2838702 Portiera (modules) | Illumina HiSeq 2500 via SRA | 97,729 pairs | 101 | 93 / 89 | 21.8 | **PASS** (195,458 reads; first written as WARN by counting pairs, corrected under D32) |
| sarscov2 test (modules) | Illumina NextSeq | 100 pairs | 150 | 88 | 24.9 | FAIL (reads) |
| ECE107 R10 (ENA SRR34410547, first 12k) | ONT R10.4.1 MinION | 12,000 | 15,463 | 72 | 16.7 | **PASS** |
| A1403KPN ONT (bacass) | ONT, older chemistry | 28,803 | 4,886 | 5 | 10.3 | WARN (N50) |
| subset15000 (bacass) | ONT | 3,750 | 21,650 | 0 | 10.3 | WARN (reads) |
| Portiera ONT (modules) | ONT via SRA, **no read UUID in headers** | 7,100 | 3,086 | 0 | 10.3 | WARN (reads, N50) |
| ERR3201952 (taxprofiler) | ONT via SRA, **headers end `/1`** | 10,000 | 1,942 | 0 | 6.5 | FAIL (quality) |
| NA12878 (nanoseq) | ONT, human | 334 | 22,913 | 7 | 9.1 | FAIL (reads) |
| SRR10382244 (modules) | **PacBio HiFi** | 204 | 18,494 | 97 | 27.3 | must be UNCERTAIN |

Findings:
- "No test sample can PASS" was a property of the small simulated CI set, not of the floors. Real runs pass on both platforms.
- mag's `test_minigut_R1/R2` are the same size as B. fragilis `test1_1/2`, so they are almost certainly the same simulated data. Not used.
- **Detection rules (D3):** SRA appends `/1` to ONT read names, so `/1` is not Illumina evidence. Some SRA ONT files carry no read UUID, so headers can be absent and detection must work from length and quality alone. Quality cannot separate platforms: PacBio HiFi (Q27) and R10 ONT (Q17, Q30 72%) overlap Illumina. **Read-length spread** separates them: Illumina lengths are near-constant (after trimming, at most 151), long-read lengths vary widely.
- ENA metadata can contradict the reads. SRR30615993 claims R9.4.1 + Guppy 2.0.10, but its reads look far better than that chemistry. Rejected as a reference.
- Arithmetic means of Phred scores overstate quality: SRR34410547 looked Q34 that way, but its error-rate-averaged quality is Q16.7. Thresholds must say which mean they use. Verify NanoPlot's definition against `mean_read_quality`.
- To verify: three different old-chemistry ONT files all give AvgQual 10.25–10.26. This could be coincidence or how seqkit averages; check with NanoPlot.
- **Screen on real file names:** 11 flags on 11 rows, all correct. No false `MATES_POSSIBLY_SWAPPED` on accession `_1`/`_2` names. One debatable flag: sample `A1403KPN_ont` with file `A1403KPN.fq.gz` gets `BASENAME_MISMATCH` (the ID has a suffix the file lacks). Decide whether that is noise.

**First end-to-end run with ECE107 (2026-10-20, local, emulate_amd64):** 29/29 tasks in 3 min, e2e checks 4/4. Tool output checked against `metric_guide.tsv`:
- **fastp: all 9 `source_field` paths exist** (ECE107_NOVASEQ: 236,248 reads after filtering, Q30 0.88, duplication 0.013, insert peak 209). But `total_reads` **counts both mates**: 120,000 pairs = 240,000 reads. The 100,000 floor therefore means 50,000 pairs, which the user confirmed is intended (D32).
- fastp single-end JSON has **no `adapter_cutting` or `insert_size` section**. The report must treat these as "not available", not as zero or an error.
- **NanoPlot: the guide's field names (`mean_qual`, `n50`, `number_of_reads` ...) do not appear in `NanoStats.txt`**, which uses labels ("Mean read quality:"). They are the keys of NanoPlot's `--tsv_stats` output. Fix when the report reads metrics: add `--tsv_stats` to `NANOPLOT` `ext.args`. ECE107_R10: mean read quality 16.6 (seqkit said 16.7), N50 15,463, 12,000 reads, so it clears every ONT floor.

## Open questions

- **BASENAME_MISMATCH when the sample ID is the file name plus a suffix** (`A1403KPN_ont` vs `A1403KPN.fq.gz`): warn (current) or accept? Found in the survey.
- **Do screen flags affect the verdict ceiling?** Currently they are listed as notes only. LAYOUT_CONFLICT, MATES_POSSIBLY_SWAPPED and ID_COLLISION arguably should cap at WARN.
- **Known gap:** an ONT-routed sample with two FASTQs makes NanoPlot fail, which stops the run under D18. Route it to `unrouted` (or join files) when detection lands.
- **Native arm64 option:** `nf-core modules install` generated `conf/containers_*_{amd64,arm64}.config` (not included by any config, so currently inert). `containers_docker_arm64.config` maps seqkit, fastp, FastQC and MultiQC to native arm64 images; NanoPlot has none (BioContainers, amd64 only). An `arm64` profile including it would speed up local M1 runs, with only NanoPlot emulated.
- Final allowed values for `sample_type` (proposed:
  `isolate | metagenome | amplicon | unknown`). Recorded and checked only in
  v0.1.
- Exact test-data URLs: confirm they exist in nf-core/test-datasets before
  use.
- Review the provisional thresholds in `assets/metric_guide.tsv` and replace
  them with lab-practice values over time (change `provenance` when you do).
- Verify the claim that Q7 is below every MinKNOW default pass filter.
  (Field names were checked on 2026-10-20: see the ECE107 run note.)
- Add a linux-64 lockfile, generated from a real install in CI, so Linux users
  can do a byte-identical rebuild too.
- A corrupt input (e.g. truncated gzip) makes seqkit exit non-zero, which
  stops the whole run under D18. Decide whether one bad file should halt the
  batch (current, conservative) or become a FAIL row in the report while the
  other samples continue.
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
- GitHub: `MacaulayWinter/seq-triage` (private). HTTPS auth via the gh keyring. `gh auth setup-git` puts the credential helper in `~/.gitconfig`, which is needed because the env ships its own git 2.56, and that git does not read Apple git's system `osxkeychain` helper.

## Status

**Session 1 (2026-10-20)**
- Done: environment checked and installed, git identity set, repo
  initialised, `.gitignore` committed (`45a93d9`), private GitHub repo
  created and pushed, CLAUDE.md written. Reproducible launch env
  (`environment.yml` + osx-arm64 lockfile), `docs/reproducibility.md`, and
  provisional `assets/metric_guide.tsv` (13 metrics) drafted.
- Skeleton so far: `nextflow.config` (manifest, params, nf-schema 2.8.0,
  strict validation, docker/emulate_amd64 profiles, run reports),
  `conf/base.config` (labels + retry policy, tested), `main.nf` +
  `nextflow_schema.json` (parameter validation, tested with 6 cases).
  Pushed up to `9fcdce3`.
- README + MIT LICENSE added. Samplesheet schema + parse-time checks (D21, D22).
  `bin/screen_metadata.R` + unit tests (D23, D24).
- **Walking skeleton runs end to end** on branch `walking-skeleton`: SCREEN_METADATA,
  SEQKIT_STATS, DETECT_PLATFORM (stub: all uncertain), join + routeSample + branch,
  FASTP/FASTQC/NANOPLOT, MULTIQC, TRIAGE_REPORT (real D2 ceiling), versions topic.
  Local run: 22/22 tasks, e2e checks 4/4. CI added (D28).
- Next: open a PR for `walking-skeleton`, get CI green, merge. Then replace the
  detection stub with real rules (D3, D26), updating `tests/e2e/expected_verdicts.tsv`
  to the final D25 answers.
- PR #1 merged after CI passed (`5bf5d13`). Test-data survey done (see "Test data survey").
  Public test-data repo created (D31), and `ECE107_NOVASEQ` + `ECE107_R10` added to the
  test samplesheet on branch `test-data-ece107` (stub-stage verdicts: WARN).
- Note: `~/Documents` is synced to iCloud ("Desktop & Documents"). The pipelines folder is
  set to "Keep Downloaded". Bulky data and Nextflow work dirs go in `~/seqtriage-scratch/`
  (not synced).
- Reminder: run `conda activate seq-triage-dev` before using
  `nextflow` or `nf-core`.

## CI fix log

Every automatic CI fix is recorded here in plain language: what failed, why, what changed, how it was verified (user rule, 2026-10-20). Fixes are never applied silently.

_No fixes yet._
