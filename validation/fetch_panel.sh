#!/usr/bin/env bash
#*******************************************************************************
# fetch_panel.sh: download the start of every panel run (CLAUDE.md D37)
#
# WHAT IT DOES  For each run with role "panel" in panel.tsv, streams the FIRST
#               N reads of each FASTQ from ENA (short reads: 50,000 per file;
#               long reads: 5,000) and stops; the rest is never downloaded.
#               If a run fails (server error, empty or broken file, mates of
#               different length), the next backup run of the SAME stratum is
#               used instead, and the swap is logged.
# WRITES        <out_dir>/fastq/*.fastq.gz, <out_dir>/fetch_log.tsv (every
#               attempt, with read counts and md5 of the uncompressed reads),
#               and <out_dir>/samplesheet.csv for the pipeline. The samplesheet
#               never declares the platform (the reads must decide), unless
#               DECLARE=ena is set (round 3, D61): then it declares ENA's
#               platform where that is illumina or ont, as a user would, so
#               mislabelled runs become real contradictions.
#
# Usage:  bash validation/fetch_panel.sh <panel.tsv> <out_dir>
#*******************************************************************************
set -u

panel=$1
out=$2
mkdir -p "$out/fastq"
log="$out/fetch_log.tsv"
sheet="$out/samplesheet.csv"
printf 'stratum\trun\tresult\tfile\treads\tmd5_uncompressed\tnote\n' > "$log"
printf 'sample,fastq_1,fastq_2,platform,sample_type,molecule\n' > "$sheet"

# Reads to keep per file: short-read platforms 50,000; long-read 5,000.
reads_for() { case $1 in ont|pacbio) echo 5000 ;; *) echo 50000 ;; esac; }

# Stream one file: first $2 reads of URL $1 into $3. Returns the read count.
# curl stops with "error 23" once head has enough; that is expected.
fetch_one() {
    # Already fetched on an earlier attempt: reuse it.
    if [ -s "$3" ] && [ -s "$3.done" ]; then cat "$3.done"; return; fi
    curl -sS --retry 3 --max-time 900 "https://$1" 2>/dev/null | gunzip -c 2>/dev/null \
        | head -n $(( $2 * 4 )) > "$3.tmp"
    local lines; lines=$(wc -l < "$3.tmp" | tr -d ' ')
    if [ "$lines" -eq 0 ] || [ $(( lines % 4 )) -ne 0 ]; then rm -f "$3.tmp"; echo 0; return; fi
    # md5 of the uncompressed reads goes to a side file: this function runs
    # in a subshell, so a variable set here would be lost.
    if command -v md5sum >/dev/null; then md5sum "$3.tmp" | cut -d' ' -f1; else md5 -q "$3.tmp"; fi > "$3.md5"
    gzip -n -c "$3.tmp" > "$3" && rm -f "$3.tmp"
    echo $(( lines / 4 )) | tee "$3.done"
}

col() { awk -F'\t' -v c="$1" 'NR==1{for(i=1;i<=NF;i++) h[$i]=i; next} {print $h[c]}' "$panel"; }

# Strata in panel order; for each, take panel runs, falling back to backups.
for stratum in $(col stratum | awk '!seen[$0]++'); do
    need=$(awk -F'\t' -v s="$stratum" 'NR>1 && $1==s && $2=="panel"' "$panel" | wc -l | tr -d ' ')
    got=0
    # Fields are split on '|', not tab: bash merges consecutive tabs (tab is
    # whitespace to IFS), so an empty field would shift every later one.
    while IFS='|' read -r run truth layout stype mol ftp; do
        [ "$got" -ge "$need" ] && break
        n=$(reads_for "$truth")
        IFS=';' read -r -a urls <<< "$ftp"
        id="${stratum}_${run}"
        ok=1; counts=()
        for i in "${!urls[@]}"; do
            f="$out/fastq/${id}_R$((i + 1)).fastq.gz"
            c=$(fetch_one "${urls[$i]}" "$n" "$f")
            counts+=("$c")
            printf '%s\t%s\t%s\t%s\t%s\t%s\t\n' "$stratum" "$run" "$([ "$c" -gt 0 ] && echo ok || echo failed)" \
                "$(basename "$f")" "$c" "$(cat "$f.md5" 2>/dev/null)" >> "$log"
            [ -s "$f.md5" ] && rm -f "$f.md5"
            [ "$c" -gt 0 ] || ok=0
        done
        if [ "$ok" -eq 1 ] && [ "${#counts[@]}" -eq 2 ] && [ "${counts[0]}" -ne "${counts[1]}" ]; then
            ok=0; printf '%s\t%s\tfailed\t\t\t\tmates have different read counts\n' "$stratum" "$run" >> "$log"
        fi
        if [ "$ok" -eq 0 ]; then
            rm -f "$out/fastq/${id}"_R*.fastq.gz
            printf '%s\t%s\tswapped_out\t\t\t\tusing next backup\n' "$stratum" "$run" >> "$log"
            continue
        fi
        r1="$out/fastq/${id}_R1.fastq.gz"
        r2=""; [ "${#urls[@]}" -eq 2 ] && r2="$out/fastq/${id}_R2.fastq.gz"
        decl=""
        if [ "${DECLARE:-}" = ena ]; then case $truth in illumina|ont) decl=$truth ;; esac; fi
        printf '%s,%s,%s,%s,%s,%s\n' "$id" "$r1" "$r2" "$decl" "$stype" "$mol" >> "$sheet"
        got=$((got + 1))
        echo "  $id: ${counts[*]} reads"
    done < <(awk -F'\t' -v s="$stratum" 'NR==1{for(i=1;i<=NF;i++) h[$i]=i; next}
               $h["stratum"]==s {print $h["run_accession"]"|"$h["truth"]"|"$h["library_layout"]"|"$h["sample_type"]"|"$h["molecule"]"|"$h["fastq_ftp"]}' "$panel")
    [ "$got" -lt "$need" ] && echo "WARNING: $stratum: only $got of $need runs fetched" >&2
done
echo "Done: $(($(wc -l < "$sheet") - 1)) samples in $sheet"
