#!/usr/bin/env bash
#*******************************************************************************
# check_repo.sh: refuse anything the hard rules forbid in this repository
#                (CLAUDE.md D52)
#
# WHAT IT DOES  Checks files for three things the repo must never hold:
#               1. sequencing data (FASTQ, BAM, CRAM, SAM, POD5, FAST5, SRA),
#                  by extension and by content (a FASTQ renamed .txt is
#                  still caught: line 1 starts "@", line 3 starts "+");
#               2. files over 1 MB, or 10 MB under validation/outcomes/
#                  (published result tables: layer D round 1's qc_values.tsv
#                  is 5.7 MB; everything else tracked is under 60 KB);
#               3. HPC details: generic patterns (SLURM accounts like
#                  research_project-..., /lustre/ and /gpfs/ paths), plus
#                  your own private strings, one per line, in
#                  .githooks/private_patterns.txt (git-ignored, so the
#                  strings themselves never enter the repo).
#               With --staged (the pre-commit hook) it also checks that the
#               commit author email is a GitHub noreply address (D10).
# MODES         --staged  the files staged for this commit (pre-commit hook)
#               --all     every tracked file (CI, and a one-off full check)
# EXIT          0 clean, 1 something found (each problem is printed).
#
# Usage:  bash .githooks/check_repo.sh --staged | --all
#*******************************************************************************
set -u

mode=${1:-}
case $mode in --staged|--all) ;; *) echo "Usage: check_repo.sh --staged | --all" >&2; exit 2 ;; esac

root=$(git rev-parse --show-toplevel)
cd "$root" || exit 2
MAX_BYTES=1048576
problems=0
problem() { echo "  BLOCKED: $*" >&2; problems=$((problems + 1)); }

# Read a file as committed: from the index (--staged) or the working tree (--all).
content() { if [ "$mode" = --staged ]; then git show ":$1" 2>/dev/null; else cat "$1"; fi; }

if [ "$mode" = --staged ]; then
    files=$(git diff --cached --name-only --diff-filter=ACMR)
else
    files=$(git ls-files)
fi

patterns='research_project-[A-Za-z0-9]+|/lustre/|/gpfs/'
private="$root/.githooks/private_patterns.txt"
private_list=()
if [ -f "$private" ]; then
    while IFS= read -r p; do
        case $p in ''|'#'*) continue ;; esac
        private_list+=("$p")
    done < "$private"
fi

while IFS= read -r f; do
    [ -n "$f" ] || continue
    # The rules themselves name the patterns; this script is exempt from 3.
    [ "$f" = .githooks/check_repo.sh ] && continue
    case $f in
        *.fastq|*.fastq.gz|*.fq|*.fq.gz|*.bam|*.bai|*.cram|*.crai|*.sam|*.pod5|*.fast5|*.sra)
            problem "$f: sequencing data file (hard rule: no research data in this repo)"; continue ;;
    esac
    size=$(content "$f" | wc -c | tr -d ' ')
    limit=$MAX_BYTES; case $f in validation/outcomes/*) limit=$((10 * MAX_BYTES)) ;; esac
    if [ "$size" -gt "$limit" ]; then
        problem "$f: $size bytes, over the $((limit / MAX_BYTES)) MB limit"; continue
    fi
    # Binary files (images, zips) skip the text checks.
    # (A NUL byte in the first 8 KB marks binary. Counted with tr, because
    # bash cannot hold a NUL in a string: grep $'\x00' would match everything.)
    head_bytes=$(content "$f" | head -c 8000 | wc -c)
    text_bytes=$(content "$f" | head -c 8000 | tr -d '\000' | wc -c)
    [ "$text_bytes" -lt "$head_bytes" ] && continue
    # FASTQ by content, whatever the extension: @header, sequence, +, quality.
    if content "$f" | head -n 4 | awk 'NR==1 && /^@/ {a=1} NR==3 && /^\+/ && a {f=1} END {exit !f}'; then
        problem "$f: looks like FASTQ reads (line 1 '@', line 3 '+')"
    fi
    hit=$(content "$f" | LC_ALL=C grep -n -E "$patterns" | head -n 3)
    [ -n "$hit" ] && problem "$f: HPC account or path (keep it in a git-ignored settings file):"$'\n'"$hit"
    for p in ${private_list[@]+"${private_list[@]}"}; do
        if content "$f" | LC_ALL=C grep -q -F -- "$p"; then
            problem "$f: contains a string from .githooks/private_patterns.txt"
        fi
    done
done <<< "$files"

if [ "$mode" = --staged ]; then
    email=$(git var GIT_AUTHOR_IDENT | sed -E 's/.*<([^>]*)>.*/\1/')
    case $email in
        *@users.noreply.github.com) ;;
        *) problem "author email $email is not a GitHub noreply address (CLAUDE.md D10)" ;;
    esac
fi

if [ "$problems" -gt 0 ]; then
    echo "check_repo.sh: $problems problem(s); nothing was committed." >&2
    echo "If a file is wrongly blocked, say why in CLAUDE.md before changing this check." >&2
    exit 1
fi
echo "check_repo.sh: clean ($(printf '%s\n' "$files" | grep -c .) files checked)"
