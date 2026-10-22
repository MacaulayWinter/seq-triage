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

## Round 2 of layers A-C: design (pre-registered 2026-10-21)

Fixed before any round 2 data was fetched or run (CLAUDE.md D51). Round 2
repeats round 1 on new data, so the pass marks are unchanged.

**What stays the same:** the three layers and their pass marks, the 18
strata, the subsample sizes (50,000 short reads, 5,000 long reads per file),
and the rule that the samplesheet never declares the platform. The round 2
run counts (115 in total) were set in `select_panel.R` before round 1 ran.
The contamination screen is off: these runs have no known contamination to
score against.

**What changes, and why:**

1. **Every round 1 study is excluded.** The D40 read-name markers were
   designed after looking at round 1 runs (the MGI flowcell name came from
   SRR24724636), so testing them on the same studies would be circular.
2. **Whole-pool draw.** Round 1 drew from the first 500 runs ENA returned
   per stratum. Keeping one run per study then left pools as small as 4
   runs, and 19 of 40 picks were DDBJ runs, which often lose their read
   names. Round 2 fetches every matching run's accession and study, shuffles
   the studies with a fixed seed, and checks one random run per study until
   the stratum is full. A stratum that cannot fill is reported as short.
3. **New pass mark for layer A:** no run of a platform seq-triage does not
   handle (PacBio, DNBSEQ, Element, Ion Torrent) may get a PASS or WARN
   verdict. This is the safety property behind the D40 markers. (Ion Torrent
   called Illumina at high confidence is already a wrong high-confidence
   call.)
4. **Layer C from a different base:** *K. pneumoniae* KPE111 (100,000
   NovaSeq pairs, 12,000 R10 reads) instead of *E. coli* ECE107. The N50
   cases aim at the same N50s from the new base. Two new cases:
   `ILM_MALFORMED` (valid gzip, one quality line a character short; expected
   FAIL with `FILE_UNREADABLE`) and `ONT_READS_10K_EXACT` (exactly on the
   10,000-read warn limit; expected PASS, because a value at a limit
   passes). `ILM_TRUNCATED_GZ` runs in the main batch, where since D45 it
   must be a FAIL row while every other sample finishes.

**The panel** ([`validation/outcomes/round2/panel.tsv`](../validation/outcomes/round2/panel.tsv)),
drawn 2026-10-21 and committed before any reads were fetched: 115 runs and
36 backups, each from a different study. Pools after exclusions held 182
(DNBSEQ) to 17,874 (Illumina amplicon) studies per stratum, against 4 to 121
eligible runs in round 1. DDBJ runs are 4 of 115 (round 1: 19 of 40). Each
pool's size and checksum is in
[`pools.tsv`](../validation/outcomes/round2/pools.tsv). The pools themselves
(69 MB of accessions) are not stored; re-running `select_panel.R` fetches
them again, and ENA's contents change over time, so a later draw can differ.
Two strata needed many runs checked to fill (Illumina RNA-seq 184, PacBio
CLR 314), because most of their runs have no FASTQ in ENA.

**Reported, not judged:** the share of correct high-confidence calls per
stratum (how often the pipeline is useful, not only safe), and confidence by
archive (the effect of DDBJ's stripped read names, open question A1).

**How to reproduce:**

```bash
Rscript validation/select_panel.R 2 panel_round2.tsv
bash validation/fetch_panel.sh panel_round2.tsv round2
Rscript validation/make_faults.R KPE111 KPE111_ILM_R1.fastq.gz KPE111_ILM_R2.fastq.gz KPE111_ONT.fastq.gz ECE107_ILM_R2.fastq.gz round2_faults 2
```

The KPE111 and ECE107 files are the layer D sources
(`validation/layer_D/fetch_sources.sh`). The pipeline and scoring commands
are the same as for round 1.

## Results: round 2 of layers A-C (2026-10-21)

115 public runs from 18 strata and 31 damaged samples, scored against the
pre-registered pass marks above. Tables:
[`validation/outcomes/round2/`](../validation/outcomes/round2/). Run on the
Mac (`-profile docker,arm64`, three batches), pipeline at `main` `e25edb0`
plus the round 2 validation scripts.

| Layer | Pass mark | Result |
|---|---|---|
| **A** Wrong at high confidence | 0 | **2: FAIL as registered** (both under review, below) |
| **A** Unhandled platform given PASS or WARN (new) | 0 | **1: FAIL as registered** (the same run as one of the two) |
| **A** Layout | all | 86 / 86 |
| **B** Independent checks | all | **1158 / 1158** ✓ |
| **C** Damaged samples, main batch | all | **30 / 30** ✓ |
| **C** `ILM_R2_SHORT` (isolated run) | WARN + `MATE_COUNT_MISMATCH` | **WARN + `MATE_COUNT_MISMATCH`** ✓ |

**As registered, round 2 fails layer A.** The two wrong high-confidence
calls are set out with their evidence below. ENA was wrong once in round 1
(SRR24724636), so each was reviewed by hand. Neither has been adjudicated:
correcting a recorded truth is the user's decision, and the evidence differs
in strength between the two.

### The two wrong high-confidence calls

**SRR28793840: ENA says Illumina HiSeq 2500; the pipeline says MGI/DNBSEQ.**

- Every read name is an MGI flowcell name (`E100056810L1C001R00200000064`),
  a format no Illumina instrument writes. A sibling run in the same study,
  also labelled HiSeq 2500 (SRR28793836), carries names from the same MGI
  flowcell (`E100056810L1C001R002...`).
- The submitter's description adds nothing ("RNA"). Quality values (34
  distinct, up to Q40) do not separate the two platforms.
- The read names are the only evidence, and they are the same signal the
  pipeline used. Correcting ENA on that alone is partly circular, so it is
  left to the user. If corrected, this call is right.

**SRR14318705: ENA says BGISEQ-500; the pipeline says Illumina (high).**

- Read names are stripped by the archive (`1/1`), so only lengths and
  quality are left. The qualities use exactly 7 values, `)-7<AFJ`, Illumina's
  HiSeq X / 4000 binning. Every MGI run in the panel that kept its names uses
  35-37 values, or one placeholder value (SRA Lite), never this set.
- The library is named `ill`, and the study's only other run is PacBio
  (a hybrid assembly, which is usually Illumina + long reads).
- **Against:** the submitter's own design description says "BGISEQ
  sequencing".
- If ENA is right, this is a **real error with a clear cause**: with read
  names gone, binned quality was the only confirming signal, and it raised
  the call to high. Because the run was called Illumina, it was routed to
  QC and got WARN, which is also the one failure of the new
  unhandled-platform mark.

**How much rests on binned quality alone:** 8 runs were called Illumina at
high confidence with binned quality as the only confirming signal (no
Illumina read names). 7 were right; the 8th is SRR14318705. All 41 calls
confirmed by Illumina read names were right.

**Scores under each reading (for the user's decision):**

| Truth used | Wrong at high | Unhandled given PASS/WARN | Layer A |
|---|---|---|---|
| ENA as recorded | 2 | 1 | FAIL |
| SRR28793840 corrected to dnbseq | 1 | 1 | FAIL |
| both corrected | 0 | 0 | pass |

Options if SRR14318705 is BGISEQ (none applied): (a) binned quality alone
gives at most medium confidence: the 7 correct calls above would become
medium, so unrouted and UNCERTAIN unless the platform is declared; (b) a
second confirming signal is required for high when read names are missing;
(c) no change, documenting the limit.

### Other findings

- **Usefulness:** 72 of 90 runs on handled platforms were called correctly at
  high confidence (round 1 re-run: 29 of 40, all platforms). Medium: 14.
  Uncertain: 3 (two 2x301 Illumina runs served already trimmed, so their
  length spread was 0.36 and 0.45; one nanopore amplicon-like metagenome).
  Uncertain is never wrong.
- **Unhandled platforms:** 24 of 25 ended UNCERTAIN. Read-name markers named
  4 of 5 DNBSEQ, 5 of 5 Element and 3 of 5 PacBio HiFi runs at high
  confidence. All 5 PacBio CLR runs, 1 HiFi run and 4 of 5 Ion Torrent runs
  were called ONT or Illumina at medium, so not routed: safe, as designed.
- **Archive:** DDBJ runs were 4 of 115 (round 1: 19 of 40). Most medium
  calls are SRA runs (13 of 77) whose read names were replaced by numbers.
- **Layer B:** every recalculated number matched, across 1,158 checks.
  Samples the integrity check rejected are skipped in layer B, because they
  have no QC numbers; their handling is tested in layer C.
- **Layer C:** all three new cases behaved as registered: the malformed and
  truncated files became FAIL rows while the batch finished, and exactly
  10,000 reads, on the limit, passed.

### Decisions after round 2 (2026-10-22, D58)

Made by the user on Claude's recommendation, after seeing the results. Round
2 therefore cannot test them, and **round 2 stays recorded as failing layer A
as registered.**

1. **SRR28793840 corrected to DNBSEQ** in
   [`adjudications.tsv`](../validation/adjudications.tsv): an MGI flowcell
   name in every read, in this run and a sibling run, is a fact about the
   file that no Illumina instrument produces.
2. **SRR14318705 not corrected**: the submitter's own description says
   "BGISEQ sequencing", so it stays a wrong call.
3. **Binned quality alone gives at most medium confidence** (it was the only
   signal behind the SRR14318705 error). The other cures listed above came
   to the same thing: without read names, no second signal is available.

Re-scored with detection re-run under the new rule
([`layer_A_D58.tsv`](../validation/outcomes/round2/layer_A_D58.tsv)). The
re-run was checked first: the old code reproduces all 115 published calls
exactly. Panel samples declare no platform, so only high-confidence
Illumina/ONT calls are routed to QC.

| | As registered | After D58 |
|---|---|---|
| Wrong at high confidence | 2 | **0** |
| Unhandled platform routed to QC | 1 | **0** of 26 |
| Correct at high confidence (handled platforms) | 72 / 90 | **65 / 89** |

The cost is the 7 correct Illumina calls that rested on binned quality alone:
they are now medium, so UNCERTAIN unless the platform is declared.
SRR14318705 becomes "wrong, medium": not routed, UNCERTAIN.

**Found by the change, then decided (D60):** the test-profile sample
`CONTRADICTION` (Illumina reads declared `ont`) is now called Illumina at
medium confidence. Under the old routing rule (D2) a medium call deferred to
the declared platform, so the Illumina reads went through nanopore QC and
FAILED on read N50: a safe verdict with a misleading reason. Before D58 a
high-confidence call overrode the label (WARN, "review"). Decided by the
user on Claude's recommendation: a medium or low call that contradicts the
declared platform is not routed. The sample is now UNCERTAIN: "Metadata says
ont but the reads suggest illumina (detection confidence: medium) ... Review
this sample." No round 2 panel sample declares a platform, so the re-score
above is unchanged.

## Round 3, layer A: testing D58 and D60 (pre-registered 2026-10-22)

Fixed before any round 3 run was selected or fetched (CLAUDE.md D61). D58
(binned quality alone gives at most medium confidence) and D60 (a medium or
low call that contradicts the declared platform is not routed) were chosen
after seeing round 2, so round 2 cannot test them. Round 3 does, on studies
neither rule has seen.

- **Exclusions:** every study in the round 1 and round 2 panel files, plus
  the calibration studies.
- **Draw:** the round 2 method (whole pool, studies in seeded random order,
  one random run per study), from round 2's saved pool snapshot
  (2026-10-21; checksums in `round2/pools.tsv`), seed offset for round 3.
- **54 runs**, weighted to where the rules act: Illumina strata whose runs
  often lose their read names (NovaSeq 6000 6, NovaSeq X 4, metagenome 4,
  metatranscriptome 4, RNA-seq 4), other Illumina (NextSeq, MiSeq, HiSeq,
  3 each), DNBSEQ 8, Element 4, Ion Torrent 3, PacBio HiFi 2 and CLR 2,
  nanopore 4 (2 older, 2 R10-era). Amplicon and nanopore RNA/metagenome
  strata are left out: the rules do not act there.
- **Declared platforms (new):** each run declares ENA's platform where it is
  illumina or ont (`DECLARE=ena`), as a user would; other platforms are left
  blank. A mislabelled run becomes a real contradiction, which D60 governs.
- **Same subsampling** (50,000 short reads, 5,000 long reads per file),
  contamination screen off, no layer C.

**Pass marks:**

| Mark | Rule |
|---|---|
| **A1** | Zero wrong calls at high confidence (truth: ENA, corrected only by a documented hand review). |
| **A2** | No run of a platform seq-triage does not handle (PacBio, DNBSEQ, Element, Ion Torrent) is routed to QC: no PASS, WARN or FAIL from QC. |
| **B** | Every independently recalculated number within its tolerance, as before. |

**Reported, not judged:** correct high-confidence calls per stratum (the
cost of D58), every contradiction and how D60 settled it, confidence by
archive.

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

### Results, round 2: new isolates, two databases (2026-10-22)

Pre-registered as D49: the same 13 mixtures, expected verdicts and pass
marks as round 1, with all-new isolates, run against the 8 GB-capped
Standard-8 and the full Standard database (Oct 2025, 95 GB). Both runs
used the same mixtures (checked: identical `truth.tsv`) and the same tool
versions; only the database differs. Both databases hold the same
libraries; Standard-8 is the same database capped to 8 GB. Tables:
[`layer_D_round2_std8/`](../validation/outcomes/layer_D_round2_std8/),
[`layer_D_round2_full/`](../validation/outcomes/layer_D_round2_full/),
side by side in
[`layer_D_round2_compare.tsv`](../validation/outcomes/layer_D_round2_compare.tsv).

| | Standard-8 | Full Standard |
|---|---|---|
| Scored samples PASS (pass mark: all 24) | **20 / 24: not met** | **22 / 24: not met** |
| Human host (10 samples) | 10 / 10 | 10 / 10 |
| Off-target (14 samples) | 10 / 14 | 12 / 14 |
| Strain limit (2) | confirmed on both | confirmed on Illumina; **not seen** on nanopore |
| Unclassified reads, Illumina / nanopore | 7.6-9.9% / 1.3-3.0% | 0.37-0.66% / 0.17-0.38% |

**Host DNA is right on both databases**, as in round 1: every human fraction
from 0.5% to 30% within tolerance, on both platforms.

**The cap explains most of the Illumina off-target error.** With Standard-8,
Bracken overestimated *K. pneumoniae* in all three Illumina spikes (true
2 / 10 / 30%: 4.3 / 14.9 / 41.2%). With the full database: 2.4 / 11.1 /
32.1%, all within tolerance. Reads that only match the capped database at
genus level or above are redistributed by Bracken to the species present,
and the capped database leaves more of them. Nanopore estimates, which do
not use Bracken, were within tolerance on both databases.

**The two misses that remain on the full database:**

- **K. pneumoniae 2%, nanopore:** 2.5% for a tolerance of +/- 0.4 points,
  just outside (Standard-8: 2.4%, inside).
- **K. variicola 10% in K. pneumoniae, Illumina:** 6.6% (Standard-8: 4.9%).
  Short reads assign part of *K. variicola* to its sister species *K.
  pneumoniae* on both databases; nanopore reads get the same mixture right
  (10.6% and 11.2%). Same-genus contamination is therefore under-estimated
  on Illumina, which is the conservative direction for a verdict only when
  it still crosses the limit (it did not here: PASS on Standard-8, WARN on
  the full database, expected WARN).

**Enterobacter is right on both databases** (Illumina 11.8 / 8.9%,
nanopore 9.3 / 9.7%). Round 1's nanopore under-detection (4.9% for 10%, on
Standard-8) was therefore a property of that isolate, CPO025, not of the
cap.

**The 1% trace-species rule (D47), on data it was not designed on:** every
control stayed below 5% off-target on both databases (0.0-2.7%). The full
database reports far longer trace tails on Illumina (18-36 species per
sample, against 4-9 on Standard-8; nanopore is similar on both, 38-95), so
without the floor more pure Illumina samples would have warned.

**The strain mixture on nanopore, full database: a database effect, not
contamination.** 30% of a second *E. coli* strain (AH21H052) should look
like pure *E. coli*, and on Standard-8 it does (0.0% off-target). On the
full database, 6.4% of species-level reads are called *K. pneumoniae*, so
the sample WARNs and the expected limit is not seen. Checked on the Mac
(D59): the AH21H052 nanopore reads (byte-identical to those used on ISCA)
were mapped competitively against AH21H052's own complete genome
(GCF_051563005.1) and *K. pneumoniae* HS11286 (GCF_000240185.1, chromosome
and plasmids). Of the 3,600 reads in the mixture, 3,362 fit their own
genome best (median identity 94.6%, 97.6% of each read aligned); only 5
(0.14%) align well to *K. pneumoniae*. 6.4% would need about a fifth of the
AH21H052 reads. So the source is not contaminated: the full database
assigns genuine *E. coli* reads of this lineage to *K. pneumoniae*, most
likely through a reference genome filed under *K. pneumoniae* that holds
*E. coli* sequence. Confirming which reads, and so which lineage, needs
Kraken2's per-read calls on ISCA
([`per_read_kraken2.sh`](../validation/layer_D/per_read_kraken2.sh)).

**What this means for use:** with the full database, the screen's host
estimates are dependable on both platforms and its species estimates are
dependable for different-genus contamination; same-genus contamination on
Illumina is under-estimated. With Standard-8, Illumina species estimates
are inflated by Bracken's redistribution and should be read as upper
bounds. The pass mark (every sample within tolerance) is not met by either
database, and is not changed after seeing these results.

### Not covered yet

- Fungal or protozoan contamination: neither database holds those
  libraries.
- Strain-level contamination: a separate method (ConFindr-style), planned.
