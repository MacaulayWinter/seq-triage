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

## Results: round 1 re-run with the D40 rules (2026-10-20)

Same 40 runs and 29 damaged samples, after the changes the user decided on
from round 1's findings: read-name markers for PacBio, MGI/DNBSEQ and
Element; an adapter checker; placeholder quality capped at WARN (CLAUDE.md
D40). Tables: [`validation/outcomes/round1_rerun/`](../validation/outcomes/round1_rerun/).

| | Round 1 | Re-run |
|---|---|---|
| **A** Wrong at high confidence | 0 | **0** ✓ (after one documented correction of ENA, below) |
| **A** Correct at high confidence | 19 / 40 | **29 / 40** |
| **A** Wrong at medium or low | 8 | **2** (Ion Torrent DRR191002; PacBio CLR DRR032990, read names stripped) |
| **A** Layout | 29 / 29 | 28 / 28 |
| **B** Independent checks | 526 / 526 | **542 / 542** ✓ |
| **C** Damaged samples | 28 / 29 | **29 / 29** ✓ |

- **PacBio, MGI/DNBSEQ and Element are now named** for every run that kept
  its read names (2 HiFi, 1 CLR, 2 DNBSEQ, 2 Element). Two older-chemistry
  nanopore runs rose from medium to high on ONT adapters.
- **Illumina runs with stripped names stay at medium.** This is by design:
  Illumina adapters do not confirm the instrument (see below).

### ENA was wrong about one run (adjudicated)

The first score of the re-run had **one wrong high-confidence call**:
SRR24724636, labelled "Illumina NovaSeq 6000" in ENA, called MGI/DNBSEQ.
The review by hand found the pipeline right and ENA wrong:

- every read name is an MGI flowcell name (`V350004723L1C001R003...`);
- the submitter's own design description says sequencing was run *"on
  DNBSEQ-G400RS (MGI, China)"* using MGI's *App-A* kit.

The correction is recorded with its evidence in
[`validation/adjudications.tsv`](../validation/adjudications.tsv). The scorer
applies it and keeps both truths in its output (`ena_truth`, `truth`,
`adjudicated`).

This changes the reading of round 1 in two ways:

1. **Round 1 scored this run as correct when it was wrong.** It was called
   Illumina at medium confidence, which matched ENA's mistaken label.
2. **It is the exact case behind the rule that Illumina adapters never raise
   confidence.** The library was made with an Illumina-compatible kit
   (NEBNext), then converted with App-A for MGI sequencing. Had Illumina
   adapters counted as evidence, data like this with stripped read names
   could be called Illumina at high confidence.

## Layer D: contamination screen accuracy

Mixtures built on purpose, so the true contamination is known. Full design,
run guide and scripts: [`validation/layer_D/`](../validation/layer_D/).

- **Sources:** 5 isolates sequenced on both platforms in one study
  (PRJNA1287288): *E. coli* ECE107 (base) and ECE170 (second strain),
  *K. pneumoniae*, *K. variicola*, *E. hormaechei*. Plus human HG002
  (Genome in a Bottle, openly consented), checked to match the bacterial
  reads in length and quality, so only contamination differs.
- **13 mixtures on both platforms (26 samples):**
  - controls;
  - human at 0.5, 2, 10 and 30%, plus 10% declared metagenome;
  - *K. pneumoniae* at 2, 10 and 30%;
  - *Enterobacter* at 10%;
  - *K. variicola* into *K. pneumoniae* (same genus);
  - a second *E. coli* strain at 30%. This one is expected to be **missed**,
    documenting the method's limit.
- **Pass marks (fixed before any real run):** estimate within 20% of the
  true value (floor 0.2 percentage points), and true contamination must be
  detected. Controls must show less than 5% off-target, and verdicts must be
  as listed in `mixtures.tsv`.
- **Rehearsed locally** end to end with nf-core's tiny test database (which
  makes the scores meaningless). The rehearsal found two problems, fixed
  before any real run:
  1. a script bug that stopped the run silently after the pipeline step;
  2. a flaw in the first tolerance. Its 2-point floor would have let a
     complete miss at 0.5% count as a success.
- **Run by the user on ISCA (2026-10-21)** with the prebuilt Kraken2
  **Standard-8** database (Oct 2025, 8 GB cap: archaea, bacteria, plasmid,
  viral, human, UniVec_Core), pipeline commit `0f5e2ff`, conda + SLURM. The
  source reads downloaded on ISCA had checksums identical to the local
  rehearsal's. Tables: [`validation/outcomes/layer_D_round1/`](../validation/outcomes/layer_D_round1/).

### Results, round 1 (pre-registered pass marks)

**19 of 24 scored samples passed, so the pass mark was not met.** The
second-strain cases behaved as expected: 30% of another *E. coli* strain read
as 1.0% (Illumina) and 2.7% (nanopore), confirming that species-level
screening cannot see it.

| | Illumina | Nanopore |
|---|---|---|
| **Human (host fraction)**, true 0.5 / 2 / 10 / 30% | 0.5 / 1.9 / 9.7 / 29.3% | 0.5 / 1.9 / 9.4 / 28.4% |
| Human scored | **6 / 6** | **6 / 6** |
| Second organism (off-target) scored | 4 / 6 | 3 / 6 |

All five failures were off-target estimates. Pure-culture controls show why:
a pure culture is never 0% off-target with this database.

| Pure culture | Illumina | Nanopore |
|---|---|---|
| *E. coli* ECE107 | 0.9% | 2.6% |
| *K. pneumoniae* KPE111 | 4.5% | **9.4% (falsely WARNed)** |

The background has three parts:
1. **A long tail of trace species:** 8 to 53 species, each under 1%, 2.4 to
   3.1% in total.
2. **Close relatives:** *K. variicola* and *K. quasipneumoniae* in the
   *K. pneumoniae* species complex, 1 to 3%.
3. **"*E. coli*" in pure *Klebsiella*:** 3.4% (nanopore) and 2.1% (Illumina).

### What the "*E. coli*" in pure *Klebsiella* is

Each read's true origin was found by aligning KPE111's reads (minimap2
2.30) to its own complete genome (GCF_051620725.1: 1 chromosome, 4
plasmids). That was compared with Kraken2's per-read call, from a separate
run on ISCA by the user.

| | Nanopore | Illumina |
|---|---|---|
| Reads Kraken2 itself called *E. coli* | 205 | 15 |
| from KPE111's chromosome | **160 (78%)** | 5 |
| from its plasmids | 45 (22%) | 10 |

- **The plasmid hypothesis is mostly wrong.** Plasmid reads are about four
  times more likely to be miscalled *E. coli* (5.7% vs 1.5% of reads), but
  most "*E. coli*" reads are ordinary *Klebsiella* chromosome reads. It is
  misclassification between related genera, worse with nanopore's higher
  error rate.
- **In Illumina, the 2.1% comes from Bracken, not Kraken2.** Kraken2 leaves
  85% of reads at genus or family level, and Bracken's species re-estimate
  then assigns some of them to *E. coli*.
- **So the background is a property of the classifier and database,** not
  of the sample.

### Change made: trace species no longer count (D47), and what it revealed

Off-target now counts only species at 1% or more of species-level reads.
The trace tail is reported but not counted. The rule was chosen after seeing
these controls, so **this round cannot test it**. Applied to the same tables
only as a preview:

- **Fixed:** `KPN_2` on both platforms and `KVAR_10` nanopore. The *E. coli*
  controls drop to 0%.
- **Still over:** pure *K. pneumoniae* nanopore, 9.4 → 6.4%. This is the
  *E. coli* background above, which a 1% floor cannot remove.
- **Revealed:** `ENT_10` nanopore reads **4.9%** for a true 10%. Only about
  half the *Enterobacter* reads get a species call with nanopore and this
  database; the old estimate was inflated by background. **The screen
  under-estimates *Enterobacter* contamination in nanopore data with
  Standard-8.**

**Next: round 2, pre-registered (D49).** The same 13 mixtures and pass marks
with different isolates ([`validation/layer_D/round2/`](../validation/layer_D/round2/)),
run against both Standard-8 and the full Standard database. Round 2 is what
tests the D47 rule, and the comparison between the two databases separates
the effect of the 8 GB cap from the classifier itself.

### Not covered by round 1

- **Contamination.** seq-triage has no contamination screen yet (v0.2).
  Validating one will need samples with known contamination: reads of a
  public, openly consented human genome (e.g. Genome in a Bottle HG002) mixed
  into bacterial reads at known fractions (1%, 5%, 20%).
- **Round 2** (~150 runs) will widen the candidate pool per category: after
  the one-run-per-study rule, some categories had only 4-11 eligible runs.
