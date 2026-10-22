#!/usr/bin/env Rscript
# Test for bin/recommend.R (D66)
#
# Synthetic triage results: one sample per rule in assets/downstream_pipelines.tsv
# and one per "no recommendation" reason. Checks each sample's recommendation,
# that every recommended pipeline gets samplesheet.csv, params.yaml and
# README.md, and that values only the user knows are FILL_IN. The outputs are
# kept in <out_dir> (argument, optional) so CI can validate each samplesheet
# against that pipeline's own schema (tests/recommend/check_schemas.py).
#
# Run from the repository root:
#   Rscript tests/recommend/test_recommend.R [out_dir]

args <- commandArgs(trailingOnly = TRUE)
out  <- if (length(args)) args[1] else file.path(tempdir(), "recommendations")
tmp  <- tempfile("rec_"); dir.create(tmp)
wt   <- function(x, f) write.table(x, f, sep = "\t", quote = FALSE, row.names = FALSE, na = "")

#        sample        verdict  routed_by   platform  decl_type    det_type     conf    molecule flags                 dominant                                                  bam  expected
cases <- read.table(header = TRUE, sep = "|", strip.white = TRUE, stringsAsFactors = FALSE, na.strings = "", text = "
sample        | verdict   | routed_by | platform | decl_type  | det_type   | conf      | molecule | flags            | dominant                                               | bam | expected
ILM_ISOLATE   | PASS      | detected  | illumina | isolate    | isolate    | high      | dna      |                  | Escherichia coli (99.0%)                               |     | nf-core/bacass
ONT_ISOLATE   | WARN      | detected  | ont      | isolate    | isolate    | medium    |          |                  | Escherichia coli (99.0%)                               |     | nf-core/bacass
ILM_META      | WARN      | detected  | illumina |            | metagenome | high      | dna      |                  | Bacteroides uniformis (20.0%)                          |     | nf-core/taxprofiler
ONT_META      | WARN      | declared  | ont      | metagenome | metagenome | medium    | dna      |                  | Bacteroides uniformis (20.0%)                          |     | nf-core/taxprofiler
ILM_16S       | PASS      | detected  | illumina | amplicon   | amplicon   | high      | dna      |                  | Escherichia coli / Shigella (30.0%)                    |     | nf-core/ampliseq
ILM_SC2       | WARN      | detected  | illumina | amplicon   | amplicon   | high      | rna      |                  | Severe acute respiratory syndrome coronavirus 2 (99%)  |     | nf-core/viralrecon
ILM_RNA       | PASS      | detected  | illumina | isolate    | uncertain  | uncertain | rna      |                  | Pseudomonas aeruginosa (98.0%)                         |     | nf-core/rnaseq
ILM_SC2_NEWNAME | WARN    | detected  | illumina | amplicon   | amplicon   | high      | dna      |                  | Betacoronavirus pandemicum (100.0% of species-level reads, Bracken) |     | nf-core/viralrecon
ILM_MUMPS     | WARN      | detected  | illumina | amplicon   | amplicon   | high      | rna      |                  | Orthorubulavirus parotitidis (100.0% of species-level reads, Bracken) |   | none
DNB_ISOLATE   | PASS      | detected  | dnbseq   | isolate    | isolate    | high      | dna      |                  | Escherichia coli (99.0%)                               |     | nf-core/bacass
DNB_META      | WARN      | detected  | dnbseq   | metagenome | metagenome | high      | dna      |                  | Bacteroides uniformis (20.0%)                          |     | nf-core/taxprofiler
ELEMENT_META  | WARN      | detected  | element  | metagenome | metagenome | high      | dna      |                  | Bacteroides uniformis (20.0%)                          |     | none
FAILED        | FAIL      | detected  | illumina | isolate    | isolate    | high      | dna      |                  |                                                        |     | none
BAM_INPUT     | PASS      | detected  | illumina | isolate    | isolate    | high      | dna      |                  |                                                        | yes | none
CONFLICT      | WARN      | detected  | illumina | isolate    | amplicon   | high      | dna      | TYPE_CONTRADICTS |                                                        |     | none
NO_TYPE       | WARN      | detected  | illumina |            | uncertain  | uncertain | dna      |                  |                                                        |     | none
ONT_AMPLICON  | WARN      | detected  | ont      | amplicon   | isolate    | medium    | dna      |                  | Severe acute respiratory syndrome coronavirus 2 (99%)  |     | none
AMP_NO_SCREEN | PASS      | detected  | illumina | amplicon   | amplicon   | high      | dna      |                  |                                                        |     | none
PB_ISOLATE    | WARN      | detected  | pacbio   | isolate    | isolate    | medium    | dna      |                  | Escherichia coli (99.0%)                               |     | none
META_RNA      | WARN      | detected  | illumina | metagenome | metagenome | high      | rna      |                  | Bacteroides uniformis (20.0%)                          |     | none
LANES_ISOLATE | PASS      | detected  | illumina | isolate    | isolate    | high      | dna      |                  | Escherichia coli (99.0%)                               |     | nf-core/bacass
LANES_META    | WARN      | detected  | illumina | metagenome | metagenome | high      | dna      |                  | Bacteroides uniformis (20.0%)                          |     | nf-core/taxprofiler
LANES_16S     | PASS      | detected  | illumina | amplicon   | amplicon   | high      | dna      |                  | Escherichia coli / Shigella (30.0%)                    |     | none
DIR_ONT       | WARN      | detected  | ont      | isolate    | isolate    | medium    | dna      |                  | Escherichia coli (99.0%)                               |     | none
HUMAN_WGS     | PASS      | detected  | illumina | isolate    | isolate    | high      | dna      |                  | Homo sapiens (99.0%)                                   |     | nf-core/sarek
LANES_HUMAN   | PASS      | detected  | illumina | isolate    | isolate    | high      | dna      |                  | Homo sapiens (99.0%)                                   |     | nf-core/sarek
YEAST_ISOLATE | WARN      | detected  | illumina | isolate    | uncertain  | uncertain | dna      |                  | Escherichia coli (60.0%)                               |     | none
SC_10X        | PASS      | detected  | illumina |            | uncertain  | uncertain |          | BARCODE_READ     |                                                        |     | nf-core/scrnaseq
METHYL_WGBS   | PASS      | detected  | illumina |            | uncertain  | uncertain | dna      | BISULFITE        |                                                        |     | nf-core/methylseq
INTERLEAVED_ISO | PASS    | detected  | illumina | isolate    | isolate    | high      | dna      | INTERLEAVED      | Escherichia coli (99.0%)                               |     | none
INDEX_MATE_ISO | WARN     | detected  | illumina | isolate    | isolate    | high      | dna      | INDEX_READ_AS_MATE | Escherichia coli (99.0%)                             |     | none
")
summ <- data.frame(sample = cases$sample, verdict = cases$verdict, routed_by = cases$routed_by,
                   detected_platform = cases$platform, declared_platform = cases$platform,
                   declared_type = cases$decl_type, detected_type = cases$det_type, type_confidence = cases$conf,
                   declared_molecule = cases$molecule, detection_flags = cases$flags,
                   detected_layout = ifelse(cases$sample == "INTERLEAVED_ISO", "interleaved", "paired"),
                   declared_target = ifelse(cases$sample %in% c("HUMAN_WGS", "LANES_HUMAN"), "9606", ""))
has_dom <- !is.na(cases$dominant)
viral   <- grepl("coronavirus|Betacoronavirus|Orthorubulavirus", cases$dominant)
mets <- rbind(data.frame(sample = cases$sample[has_dom], metric = "dominant_species", value = cases$dominant[has_dom]),
              data.frame(sample = cases$sample[has_dom], metric = "viral_fraction",
                         value = ifelse(viral[has_dom], "100.0%", "0.10%")),
              # D76: a yeast-like isolate: few classified reads are bacterial.
              data.frame(sample = "YEAST_ISOLATE", metric = "bacterial_fraction", value = "2.0%"))
paired <- cases$platform %in% c("illumina", "dnbseq", "element")
sheet <- data.frame(sample = cases$sample,
                    fastq_1 = ifelse(is.na(cases$bam), sprintf("/data/%s_R1.fastq.gz", cases$sample), ""),
                    fastq_2 = ifelse(is.na(cases$bam) & paired, sprintf("/data/%s_R2.fastq.gz", cases$sample), ""),
                    bam = ifelse(is.na(cases$bam), "", sprintf("/data/%s.bam", cases$sample)), fastq_dir = "")
# D75: LANES_* samples have a second lane (a second row); DIR_ONT is a folder.
lanes <- sheet[startsWith(sheet$sample, "LANES_"), ]
lanes$fastq_1 <- sub("_R1", "_L002_R1", lanes$fastq_1); lanes$fastq_2 <- sub("_R2", "_L002_R2", lanes$fastq_2)
sheet <- rbind(sheet, lanes)
sheet[sheet$sample == "DIR_ONT", c("fastq_1", "fastq_dir")] <- c("", "/data/fastq_pass/barcode01")
wt(summ, file.path(tmp, "summary.tsv")); wt(mets, file.path(tmp, "metrics.tsv"))
write.csv(sheet, file.path(tmp, "samplesheet.csv"), row.names = FALSE, quote = FALSE)
unlink(out, recursive = TRUE)
status <- system2("Rscript", c("--vanilla", "bin/recommend.R", file.path(tmp, "summary.tsv"), file.path(tmp, "metrics.tsv"),
                               file.path(tmp, "samplesheet.csv"), "assets/downstream_pipelines.tsv", out))
if (status != 0) stop("recommend.R exited with status ", status)

rec <- read.delim(file.path(out, "recommendations.tsv"), colClasses = "character", na.strings = "")
res <- data.frame(case = cases$sample, got = rec$pipeline[match(cases$sample, rec$sample)], expected = cases$expected)
# Every recommended pipeline: three files, FILL_IN where the user must decide.
for (rl in unique(rec$rule[rec$pipeline != "none"])) {
  d <- file.path(out, rl)
  ok <- all(file.exists(file.path(d, c("samplesheet.csv", "params.yaml", "README.md"))))
  res <- rbind(res, data.frame(case = paste0("FILES_", rl), got = if (ok) "3 files" else "missing", expected = "3 files"))
}
# One folder per rule, so short- and long-read bacass each run as they are.
for (k in c("short", "long")) {
  y <- readLines(file.path(out, paste0("bacass_", k), "params.yaml"))
  res <- rbind(res, data.frame(case = paste0("BACASS_", toupper(k), "_ASSEMBLY_TYPE"),
                               got = if (any(y == sprintf("assembly_type: %s", k))) k else "wrong", expected = k))
}
fill <- function(p, key) any(grepl(sprintf("^%s: FILL_IN$", key), readLines(file.path(out, p, "params.yaml"))))
# Free-text values must NOT be written as a live FILL_IN (ampliseq would run
# with the "primer" FILL_IN); they appear commented out, so ampliseq stops.
for (key in c("FW_primer", "RV_primer")) {
  y <- readLines(file.path(out, "ampliseq", "params.yaml"))
  ok <- !any(grepl(sprintf("^%s:", key), y)) && any(grepl(sprintf("^# %s: ", key), y))
  res <- rbind(res, data.frame(case = sprintf("COMMENTED_ampliseq_%s", key), got = if (ok) "commented" else "live", expected = "commented"))
}
checks <- list(c("viralrecon", "primer_set"),
               c("rnaseq", "fasta"), c("taxprofiler", "databases"))
for (k in checks) res <- rbind(res, data.frame(case = sprintf("FILL_IN_%s_%s", k[1], k[2]),
                                               got = if (fill(k[1], k[2])) "FILL_IN" else "set", expected = "FILL_IN"))
# D75: a sample in two lanes is two rows in bacass's and taxprofiler's
# samplesheets (taxprofiler: one run_accession per row).
bs <- read.csv(file.path(out, "bacass_short", "samplesheet.csv"), colClasses = "character")
res <- rbind(res, data.frame(case = "LANES_BACASS_TWO_ROWS", got = as.character(sum(bs$ID == "LANES_ISOLATE")), expected = "2"))
tp <- read.csv(file.path(out, "taxprofiler", "samplesheet.csv"), colClasses = "character")
res <- rbind(res, data.frame(case = "LANES_TAXPROFILER_RUNS", got = paste(tp$run_accession[tp$sample == "LANES_META"], collapse = ","),
                             expected = "LANES_META_1,LANES_META_2"))
am <- read.csv(file.path(out, "ampliseq", "samplesheet.csv"), colClasses = "character")
res <- rbind(res, data.frame(case = "LANES_NOT_IN_AMPLISEQ", got = as.character(sum(am$sample == "LANES_16S")), expected = "0"))
# D76: sarek gets one row per lane, lanes distinct (its uniqueEntries).
sk <- read.csv(file.path(out, "sarek", "samplesheet.csv"), colClasses = "character")
res <- rbind(res, data.frame(case = "SAREK_LANES", got = paste(sk$lane[sk$sample == "LANES_HUMAN"], collapse = ","), expected = "L1,L2"))
res <- rbind(res, data.frame(case = "SAREK_UNIQUE_ENTRIES", got = as.character(!anyDuplicated(sk[, c("patient", "sample", "lane")])), expected = "TRUE"))
for (k in list(c("scrnaseq", "fasta"), c("scrnaseq", "gtf"), c("methylseq", "fasta")))
  res <- rbind(res, data.frame(case = sprintf("FILL_IN_%s_%s", k[1], k[2]), got = if (fill(k[1], k[2])) "FILL_IN" else "set", expected = "FILL_IN"))
res <- rbind(res, data.frame(case = "SC_10X_ASSUMES_RNA", got = as.character(grepl("RNA assumed", rec$reason[rec$sample == "SC_10X"])), expected = "TRUE"))
res <- rbind(res, data.frame(case = "FILL_IN_sarek_tools", got = if (fill("sarek", "tools")) "FILL_IN" else "set", expected = "FILL_IN"))
# A reason is always given.
res <- rbind(res, data.frame(case = "EVERY_ROW_HAS_A_REASON", got = if (all(nzchar(rec$reason))) "yes" else "no", expected = "yes"))
res$result <- ifelse(res$got == res$expected, "PASS", "FAIL")
print(res, row.names = FALSE, right = FALSE)
cat(sprintf("\n%d / %d cases as expected\n", sum(res$result == "PASS"), nrow(res)))
if (any(res$result == "FAIL")) quit(status = 1)
