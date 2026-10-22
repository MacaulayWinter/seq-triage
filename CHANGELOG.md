# Changelog

All notable changes to seq-triage. Versions follow
[semantic versioning](https://semver.org/); each entry names the decisions
(`D` numbers) recorded in [CLAUDE.md](CLAUDE.md).

## v0.2.0dev (unreleased)

### Added

- PacBio HiFi QC (D64): reads named PacBio AND marked HiFi (read names
  ending `/ccs`, or a BAM written by `ccs`) are routed to NanoPlot and judged
  with HiFi limits (mean read quality warn below Q20, fail below Q10; read
  count and N50 as for nanopore); Kraken2 screens them like nanopore.
  `platform` accepts `pacbio`; a PacBio declaration alone never routes.
  Other PacBio data (CLR, subreads) is still recognised and not QC'd.

- Sample type from the reads (D65): isolate / metagenome / amplicon /
  uncertain, with confidence; contradicting a declared type caps at WARN; an
  inferred type never relaxes a verdict. Validation layer E pre-registered.
- Downstream recommendation (D66): for each ready sample, the next nf-core
  pipeline at a pinned version (bacass 2.6.1, taxprofiler 2.0.1, ampliseq
  2.18.0, viralrecon 3.0.0, rnaseq 3.27.0; rules in
  `assets/downstream_pipelines.tsv`), with a samplesheet in that pipeline's
  format and a params file. CI validates every generated samplesheet
  against the pipeline's own schema at the pinned commit.
- After validation layer E (D67): isolate at high confidence also needs at
  most 25% of reads unclassified; a virus-dominated sample is isolate at
  most at medium confidence; new information metric `viral_fraction`.

- Input formats (D69): plain FASTQ and `.bz2`/`.xz`/`.zst` (rewritten as
  `.fastq.gz` first, `NORMALISE_FASTQ`; a file that will not decompress is a
  FAIL row); SAM and CRAM in the `bam` column (`--cram_fasta`, or EBI's
  reference registry; a missing reference is a FAIL row, not a crash).
- MGI/DNBSEQ and Element are QC'd (D69): short-read QC like Illumina, own
  guide rows, MGI adapters given to fastp; `platform` accepts `dnbseq`/`mgi`/
  `bgiseq` and `element`/`aviti`; the downstream rules accept them where the
  next pipeline does.

- Several files per sample (D75): repeat the sample ID on several rows
  (lanes, top-up runs) or give a folder in the new `fastq_dir` column (a
  nanopore barcode folder). Every file is checked; one damaged file fails
  the whole sample. The platform is detected per lane and combined (new
  local step `COMBINE_DETECTION`), then the files are joined with nf-core
  `cat/fastq` before QC.

### Changed

- A contradiction between detected and declared platform (D60) is now short
  reads vs long reads: a medium "illumina" call no longer contradicts a
  declared MGI or Element sample (D69).
- `CHECK_BAM` decodes every base (`samtools fastq`), not only counts reads:
  counting passed a CRAM whose reference was missing (D69).
- The report takes each metric's limits from the guide row it was judged
  by, so one metric can have different limits per platform (D64).
- After validation layer E round 2 (D73): a metagenome call is high
  confidence only if the most common genus holds under 50% of genus-level
  reads; a new phage or poorly-known bacterium split over sister species is
  now metagenome at medium confidence. Layer E round 2 failed as registered
  (7 wrong at high; 2 with D73).

### Fixed

- The virus rule (D67) matched Viruses at rank D; the Oct 2025 databases
  report it as R1, so it never fired on real data. Now matched by taxid (D70).
- `validation/fetch_panel.sh` kept no checksums when a re-run reused files
  (D73).

## v0.1.0 (2026-10-22)

The first release. It also contains two features first planned for v0.2,
the contamination screen and BAM input; both are tested and documented
below, so they are released rather than held back.

### Added

- Samplesheet validation with nf-schema: hard errors before any data is read
  (D19, D21), then a metadata screen in R across the whole sheet (D23).
- Platform and layout detection from the reads, with a confidence level and
  `uncertain` rather than a guess (D33). Read-name markers name PacBio,
  MGI/DNBSEQ and Element runs, which are not routed (D40).
- Conservative routing: no sample can PASS unless the reads confirm its
  metadata (D1, D2).
- Routed QC: fastp + FastQC for Illumina, NanoPlot for nanopore, MultiQC over
  all, from pinned nf-core modules.
- Plain-language report: per-sample PASS / WARN / FAIL / UNCERTAIN with every
  reason, each metric explained from `assets/metric_guide.tsv` (D12, D34).
- `molecule` column (dna / rna) and RNA-aware judgements (D36).
- Integrity check: a damaged file becomes a FAIL row and the batch carries on
  (D45).
- Contamination screen at species level with Kraken2 and Bracken, off unless
  `--kraken2_db` is given (D41); trace species under 1% do not count as
  off-target (D47); `--kraken2_batch` loads a large database once (D50).
- Profiles: `docker`, `emulate_amd64` and native `arm64` for Apple Silicon
  (D9, D44), `conda` and `slurm` for HPCs without containers (D42).
- BAM input (`bam` column), aligned or unaligned: checked from header to end
  (`CHECK_BAM`), converted with nf-core `samtools/collatefastq`; layout from
  the reads' READ1/READ2 flags; a damaged BAM becomes a FAIL row (D55).
- BAM header: `@PG` software of one platform counts as detection evidence;
  `@RG PL` stands in for a blank samplesheet platform, and a disagreement
  with the samplesheet is flagged (D57).
- `--help` generated from the parameter schema (D53).
- Tests: R unit tests, nf-test module and pipeline tests (D46), end-to-end
  runs with Docker and conda, lockfile rebuild, repository guard (D52).
- Repository guard: pre-commit hook and CI job against sequencing data, HPC
  details and non-noreply emails (D52); nf-core lint clean in CI (D53).

### Changed after validation round 2 (decided after seeing the data)

- Binned quality alone gives at most medium confidence (D58).
- A medium or low call that contradicts the declared platform is not routed;
  the sample is UNCERTAIN, for review (D60).

### Validation ([docs/validation.md](docs/validation.md))

- Round 1 (40 runs, 29 damaged samples): every pass mark met (D39, D40).
- Round 2 (115 runs from a whole-ENA draw, 31 damaged samples): layers B and
  C met; **layer A failed as registered** (2 wrong at high confidence; D54).
- Round 3 (54 runs from unseen studies, pre-registered to test D58): **every
  pass mark met**; D58 stopped one wrong call it was not designed on (D62).
- Contamination (layer D): host DNA right in every case; the all-samples
  pass mark was not met in either round (round 1 19/24; round 2 20/24 with
  Standard-8, 22/24 with the full Standard database; D47, D56).

### Known limitations

- PacBio, MGI/DNBSEQ, Element and Ion Torrent are recognised (or left
  uncertain) but not QC'd.
- D60 has not yet met a real mislabelled run; it is tested by the test
  profile and a unit test only.
- Contamination: another strain of the same species cannot be detected;
  same-genus contamination is under-estimated on Illumina; with Standard-8,
  Illumina species estimates are inflated (read them as upper bounds); the
  full database assigns some E. coli of one lineage to K. pneumoniae (D59).
- CRAM input is not accepted yet.
