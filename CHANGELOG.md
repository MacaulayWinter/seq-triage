# Changelog

All notable changes to seq-triage. Versions follow
[semantic versioning](https://semver.org/); each entry names the decisions
(`D` numbers) recorded in [CLAUDE.md](CLAUDE.md).

## v0.1.0 (unreleased)

Released once validation round 2 of layers A-C (D51) has run.

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
- `--help` generated from the parameter schema (D53).
- Tests: R unit tests, nf-test module and pipeline tests (D46), end-to-end
  runs with Docker and conda, lockfile rebuild, repository guard (D52).
- Validation: round 1 of layers A-C (D39), layer D contamination rounds 1-2
  (D43, D47, D49).
