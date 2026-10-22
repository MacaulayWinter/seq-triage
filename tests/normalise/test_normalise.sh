#!/usr/bin/env bash
# Test for bin/normalise_fastq.sh (D69). Run inside the pinned samtools image
# (gzip, bzip2, xz, zstd, bgzip), from the repository root:
#   docker run --rm -v "$PWD":/w -w /w <samtools image> bash tests/normalise/test_normalise.sh
# Builds the same two reads in every format, plus a truncated .bz2 and .xz,
# and checks: every good file comes out as identical .fastq.gz content, and
# every damaged one is reported "unreadable" (never a silently short file).
set -u
tmp=$(mktemp -d); trap 'rm -rf -- "$tmp"' EXIT
script=$(pwd)/bin/normalise_fastq.sh
cd "$tmp"
printf '@r1\nACGT\n+\nIIII\n@r2\nGGCC\n+\nIIII\n' > ref.fastq
cp ref.fastq plain.fastq; cp ref.fastq plain.fq
cp ref.fastq b.fastq && bzip2 b.fastq
cp ref.fastq x.fastq && xz x.fastq
cp ref.fastq z.fastq && zstd -q --rm z.fastq
gzip -c ref.fastq > g.fastq.gz
head -c 30 b.fastq.bz2 > bad.fastq.bz2
head -c 40 x.fastq.xz > bad.fastq.xz
pass=0; fail=0
check() { if [ "$2" = "$3" ]; then pass=$((pass + 1)); printf 'PASS  %-16s %s\n' "$1" "$2"; else fail=$((fail + 1)); printf 'FAIL  %-16s got %s, expected %s\n' "$1" "$2" "$3"; fi; }
for f in plain.fastq plain.fq b.fastq.bz2 x.fastq.xz z.fastq.zst g.fastq.gz; do
    s=S_${f//./_}; bash "$script" "$s" 1 "$f" >/dev/null 2>&1
    st=$(tail -1 "$s.normalise.tsv" | cut -f2)
    same=$([ -e "${s}_1.fastq.gz" ] && gzip -dc "${s}_1.fastq.gz" | cmp -s - ref.fastq && echo identical || echo different)
    check "$f" "$st/$same" "ok/identical"
done
for f in bad.fastq.bz2 bad.fastq.xz; do
    s=S_${f//./_}; bash "$script" "$s" 1 "$f" >/dev/null 2>&1
    check "$f" "$(tail -1 "$s.normalise.tsv" | cut -f2)/$([ -e "${s}_1.fastq.gz" ] && echo output || echo none)" "unreadable/none"
done
# A pair: both mates come out, in order.
cp b.fastq.bz2 m_R1.fastq.bz2; cp x.fastq.xz m_R2.fastq.xz
bash "$script" PAIR 1 m_R1.fastq.bz2 m_R2.fastq.xz >/dev/null 2>&1
check "pair" "$(ls PAIR_*.fastq.gz | wc -l | tr -d ' ') files/$(cut -f2 PAIR.normalise.tsv | tail -n +2 | sort -u)" "2 files/ok"
# D79: FASTA reads (multi-line, Windows line ends, blank line) become FASTQ
# with a placeholder quality, in every compression.
printf '>r1 first read\nAC\nGT\n>r2\r\nGGCC\r\n\n' > reads.fasta
printf '@r1 first read\nACGT\n+\nIIII\n@r2\nGGCC\n+\nIIII\n' > fasta_ref.fastq
cp reads.fasta reads.fa; cp reads.fasta reads.fna
gzip -c reads.fasta > reads.fasta.gz; bzip2 -c reads.fasta > reads.fa.bz2; xz -c reads.fasta > reads.fna.xz
for f in reads.fasta reads.fa reads.fna reads.fasta.gz reads.fa.bz2 reads.fna.xz; do
    s=F_${f//./_}; bash "$script" "$s" 1 "$f" >/dev/null 2>&1
    st=$(tail -1 "$s.normalise.tsv" | cut -f2)
    same=$([ -e "${s}_1.fastq.gz" ] && gzip -dc "${s}_1.fastq.gz" | cmp -s - fasta_ref.fastq && echo identical || echo different)
    check "$f" "$st/$same" "ok/identical"
done
head -c 25 reads.fasta.gz > bad.fasta.gz
bash "$script" BADFA 1 bad.fasta.gz >/dev/null 2>&1
check "bad.fasta.gz" "$(tail -1 BADFA.normalise.tsv | cut -f2)/$([ -e BADFA_1.fastq.gz ] && echo output || echo none)" "unreadable/none"
echo; echo "$pass / $((pass + fail)) cases as expected"
[ "$fail" -eq 0 ]
