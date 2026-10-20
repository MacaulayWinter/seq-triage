# Validation

seq-triage is checked on data it was not tuned on, against answers known
independently of the pipeline. This page gives the method, the pass marks
(fixed before any result was seen) and the results. Decision record:
CLAUDE.md D37.

## Method

Three layers, each with a different source of truth:

| Layer | Question | Truth comes from | Pass mark |
|---|---|---|---|
| **A** | Is the platform call right? | The instrument platform recorded in ENA | **Zero wrong calls at high confidence.** *Uncertain* is never wrong. |
| **B** | Are the reported numbers right? | An independent recalculation from the FASTQ in base R, using none of the pipeline's tools | Counts exact; NanoPlot mean read quality within 0.1; fastp Q30 within 0.0001; report values within their display precision |
| **C** | Do damaged files get the verdict they should? | The damage itself: each fault is made deliberately, so its correct verdict is known in advance | Every expected verdict met and every expected flag present |

### Rules for choosing data

1. **Not tuned on.** Runs from the studies used to calibrate detection or to
   choose test data are excluded (listed in `validation/select_panel.R`).
2. **Chosen by rule.** For each category, ENA is queried for candidates in a
   stated date window. One run per study is kept, so a single large study
   cannot dominate, and runs are drawn at random with a fixed seed. Anyone
   re-running `select_panel.R` gets the same panel.
3. **No human samples.** Human genetic data is not needed to test triage.
   Metagenomes from human sites (e.g. gut) are included. They may contain
   some human reads, which matters for the v0.2 contamination screen.
4. **Subsamples.** The first 50,000 reads per file for short reads and the
   first 5,000 for long reads. This tests whether the pipeline treats the
   data correctly, not whether the original run was good: a subsample's
   verdict says nothing about the full run.

### Layer C cases

Built by `validation/make_faults.R` from the two known-good ECE107 samples in
[seq-triage-test-data](https://github.com/MacaulayWinter/seq-triage-test-data).
Each case changes one thing, and its expected verdict and flags are written
before the run, in `expected.tsv`. Two cases (`ILM_R2_SHORT`,
`ILM_TRUNCATED_GZ`) may stop a whole run by design, so each runs on its own.

### How to reproduce

```bash
Rscript validation/select_panel.R 1 panel_round1.tsv
bash validation/fetch_panel.sh panel_round1.tsv round1
Rscript validation/make_faults.R <seq-triage-test-data checkout> <another Illumina R2> round1_faults
nextflow run main.nf -profile docker --input round1/samplesheet.csv --outdir results_panel
nextflow run main.nf -profile docker --input round1_faults/samplesheet_main.csv --outdir results_faults
Rscript validation/independent_metrics.R round1/samplesheet.csv results_panel layer_B_panel.tsv
Rscript validation/independent_metrics.R round1_faults/samplesheet_main.csv results_faults layer_B_faults.tsv
Rscript validation/score_validation.R panel_round1.tsv results_panel round1_faults results_faults layer_B.tsv scores
```

Add `emulate_amd64` to `-profile` on Apple Silicon. The data stays outside
the repository.

## Results: round 1 (2026-10-20)

40 public runs from 18 categories, plus 29 damaged samples built from
ECE107. Pipeline at commit `c8e3b56` plus the fixes on the
`validation-round1` branch. Every table is in
[`validation/outcomes/round1/`](../validation/outcomes/round1/): `panel.tsv`
(what was chosen, from what pool), `fetch_md5.tsv` (checksums of every
subsample), `layer_A.tsv`, `layer_B.tsv`, `layer_C.tsv` and
`expected_faults.tsv` (written before the run).

| Layer | Result | Pass mark | |
|---|---|---|---|
| **A** Platform calls | **0 wrong at high confidence** (40 runs); layout right for 29 of 29 called runs | 0 | ✓ |
| **B** Independent recalculation | **526 / 526** checks agree | all | ✓ (see finding B1) |
| **C** Damaged samples | **28 / 29** as expected (26/27 main batch, 2/2 isolated) | all | ✗ one design question (C1) |

### Layer A: platform calls

| Category (2 runs each; ONT older / R10 era: 4) | Calls |
|---|---|
| Illumina NextSeq 2000, NovaSeq X, metagenome, RNA-seq, metatranscriptome | all correct, high |
| Illumina MiSeq, HiSeq 2500, NovaSeq 6000, amplicon | correct; 5 of 8 at medium (see A1) |
| ONT R10 era, ONT RNA | all correct, high |
| ONT older chemistry | correct, all medium (A1) |
| ONT metagenome | 1 correct (medium), 1 uncertain |
| PacBio HiFi (Revio) | 1 uncertain (right answer), 1 called ONT at **low** (A2, A3) |
| PacBio CLR (RS II) | both called ONT at **medium** (A4) |
| DNBSEQ, Element AVITI | all called illumina at **medium** |
| Ion Torrent | 1 uncertain (right answer), 1 called illumina at **medium** (A5) |

Every wrong call was at medium or low confidence. With no platform declared,
medium and low calls are not acted on: these samples got no QC and the
verdict UNCERTAIN. So no wrong call reached a PASS, which is the purpose of
the confidence ladder (D2, D33).

**Findings**

- **A1. DDBJ strips read headers.** DDBJ-hosted runs (`DRR…`) have names like
  `@DRR315135.1 1/1`, with no instrument information left. Detection then
  relies on read length alone and reports medium confidence. This is
  under-confident, never wrong. It is the main reason real Illumina and ONT
  runs fell short of high.
- **A2. PacBio headers are recognisable but not used.** Revio names reads
  `m84079_…/ccs`. Adding this pattern as evidence *against* Illumina and ONT
  would turn these calls into *uncertain*.
- **A3. Modern HiFi is not narrow.** Revio length spread was 0.43 and 0.59,
  not the 0.16 of the nf-core test file the "long but uniform = PacBio" rule
  was calibrated on. That rule does not generalise; A2 is the better signal.
  Revio also bins quality (7 values), which counts as Illumina-like evidence.
- **A4. PacBio CLR looks like nanopore** (variable lengths, headers stripped).
  It cannot be told apart from the reads alone. A known limit.
- **A5. One Ion Torrent run sits inside the Illumina limits** (longest read
  504 bp, spread 0.23). A known limit.

### Layer B: independent recalculation

All 526 comparisons agree: reads, bases and N counts are exact for every
sample and tool. NanoPlot read N50, median length and mean quality match.
fastp before-filtering reads, bases and Q30 match. The report shows the tools'
numbers faithfully.

- **B1. NanoPlot's "mean read quality" is a two-step average.** The first
  recalculation used the plain mean of per-read Q values. It disagreed with
  NanoPlot by 0.45-1.9 Q on all six R10-era and RNA runs, while agreeing on
  older chemistry. Testing four definitions on the actual files showed that
  NanoPlot averages each read's quality *as error probabilities* across
  reads. That matches its reported mean to the printed decimal; its median
  matches the per-read median exactly. The recalculation was changed to this
  definition, and **this change was made after seeing results**, so both are
  reported here. The pipeline itself was right. The metric guide's wording
  was ambiguous and has been made precise, because the definition changes
  how the Q7/Q10 limits read: NanoPlot's value is lower than a plain mean.

### Layer C: damaged samples

| Result | Cases |
|---|---|
| ✓ Expected verdict and flags | both controls; reads too few (4); N50 too low, including the amplicon exemption (3); Q30 too low (2); nanopore quality too low (2); duplication, with the RNA exemption (2); adapters; N bases (2); empty file; wrong or missing platform (2); mixed mates; swapped files; nanopore in two files; U bases with dna / rna (2); R2 shorter than R1 (isolated) |
| ✗ | **ONT_PLACEHOLDER_Q**: expected WARN, got PASS (C1) |
| ✓ but via a hidden bug | **ILM_TRUNCATED_GZ** (isolated): the run stops, as expected (C2) |

- **C1. Placeholder quality can PASS.** Every base had the same quality, so
  the quality metrics were rightly "not judged". But the read headers still
  confirmed nanopore at high confidence, and nothing else was wrong, so the
  sample passed. The expected answer, written before the run, was WARN:
  quality that cannot be assessed should not count as passing (D1).
  **Decision for the user.**
- **C2. seqkit fails silently on a truncated gzip.** It prints
  `unexpected EOF` but exits 0 and writes no row for the file. The run still
  stopped, because detection crashed on the empty table, but with an
  unhelpful R error. **Fixed:** detection now checks that seqkit returned a
  row for every file, and stops with "seqkit could not read <file>. The file
  is probably truncated or corrupt…". A regression test was added
  (`tests/detect_platform/`, case `SEQKIT_SILENT_FAILURE`).
- **Construction bugs caught before the run.** Two damaged samples first
  missed their targets: mean quality 9.5 instead of 5.5, and adapters in 15%
  of pairs instead of 40%. They were found by checking the values the
  builder achieved, fixed, and rebuilt before any pipeline result existed.
  The builder now stops if a case misses its target. Its quality targets use
  the plain per-read mean (see B1), so the nanopore-quality cases sit
  somewhat lower in NanoPlot's terms, still on the intended side of each
  limit.

### Not covered by round 1

- **Contamination.** seq-triage has no contamination screen yet (v0.2).
  Validating one will need samples with known contamination: reads of a
  public, openly consented human genome (e.g. Genome in a Bottle HG002) mixed
  into bacterial reads at known fractions (1%, 5%, 20%).
- **Round 2** (~150 runs) will widen the candidate pool per category: after
  the one-run-per-study rule, some categories had only 4-11 eligible runs.
