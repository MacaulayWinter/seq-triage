#!/usr/bin/env python3
"""kmer_redundancy.py: how many genomes do these reads come from? (D86 calibration)

WHAT IT DOES  Reads up to BUDGET bases from the start of each FASTQ and counts
              a fixed sample of its 21-mers: the canonical k-mers (the smaller
              of a k-mer and its reverse complement) that start with "ACG",
              about 1 in 32. The sample is chosen by sequence, not position,
              so the same k-mer is kept wherever and on whichever strand it
              occurs, and repeats are counted correctly.
WHY           Reads of one genome (an isolate, even one no database knows)
              cover it more than once, so its k-mers recur; reads of a
              community come from hundreds of megabases, so almost every k-mer
              is seen once. No reference database is involved.
REPORTS       per file: bases used, sampled distinct k-mers, the share of them
              seen at least twice (frac_repeated), and a rough genome size:
              distinct k-mers seen at least twice x 32 (only meaningful when
              the genome is covered at least about twice).
NEEDS         Python 3 standard library only.
Usage:        python3 kmer_redundancy.py <out.tsv> <reads.fastq.gz> [more ...]
"""
import gzip
import re
import sys
from collections import Counter

K = 21
BUDGET = 10_000_000          # bases read per file, so files are compared at equal depth
PREFIX = "ACG"               # canonical k-mers starting with this are kept (~1/32)
SAMPLE_FACTOR = 32
COMP = str.maketrans("ACGT", "TGCA")
FWD = re.compile("(?=ACG)")  # a forward k-mer starting here begins with ACG
REV = re.compile("(?=CGT)")  # a k-mer ending here in CGT has a reverse complement starting with ACG


def sampled_kmers(seq):
    """Yield the canonical sampled k-mers of one read (ACGT only)."""
    n = len(seq)
    for m in FWD.finditer(seq):
        i = m.start()
        if i + K <= n:
            kmer = seq[i:i + K]
            rc = kmer.translate(COMP)[::-1]
            if kmer <= rc and "N" not in kmer:
                yield kmer
    for m in REV.finditer(seq):
        j = m.start() + 3 - K          # start of the k-mer that ends with CGT
        if j >= 0:
            kmer = seq[j:j + K]
            rc = kmer.translate(COMP)[::-1]
            if rc < kmer and "N" not in kmer:
                yield rc


def scan(path):
    counts = Counter()
    used = 0
    with gzip.open(path, "rt", errors="replace") as fh:
        for i, line in enumerate(fh):
            if i % 4 != 1:
                continue
            seq = line.strip().upper()
            if len(seq) < K:
                continue
            counts.update(sampled_kmers(seq))
            used += len(seq)
            if used >= BUDGET:
                break
    distinct = len(counts)
    repeated = sum(1 for c in counts.values() if c >= 2)
    return used, distinct, (repeated / distinct if distinct else float("nan")), repeated * SAMPLE_FACTOR / 1e6


def main():
    out, files = sys.argv[1], sys.argv[2:]
    with open(out, "w") as fh:
        fh.write("file\tbases\tsampled_distinct\tfrac_repeated\test_genome_mb\n")
        for f in files:
            try:
                used, distinct, frac, size = scan(f)
            except (OSError, EOFError) as err:
                print(f"skipped {f}: {err}", file=sys.stderr)
                continue
            fh.write(f"{f.rsplit('/', 1)[-1]}\t{used}\t{distinct}\t{frac:.4f}\t{size:.2f}\n")
            fh.flush()


if __name__ == "__main__":
    main()
