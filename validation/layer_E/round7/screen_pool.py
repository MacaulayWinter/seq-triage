#!/usr/bin/env python3
"""Layer E round 7 (CLAUDE.md D94): keep the G-dominated runs of a random pool.

WHAT IT DOES  Reads the pool panel (select_panel.R, its rows already in a
              seeded random order) and, in that order, streams the first 2,000
              reads of each run's first FASTQ (R1, the read detection checks)
              from ENA, measuring the share of A, C, G and T. The first 12 runs
              with G >= 50% become role "panel", the next 4 "backup" (used only
              if a download fails); every other run is "screened_out". No other
              choice is made, so no run is picked by hand.
WRITES        <out_panel.tsv>: the pool's rows and columns, the new roles, and
              columns a_share, c_share, g_share, t_share, probe_reads.

Usage: python3 validation/layer_E/round7/screen_pool.py <pool_panel.tsv> <out_panel.tsv>
"""
import csv, gzip, io, sys, urllib.request

PANEL_N, BACKUP_N, G_MIN, PROBE = 12, 4, 0.50, 2000

def shares(url):
    try:
        with urllib.request.urlopen("https://" + url.removeprefix("ftp://"), timeout=120) as r:
            seqs = []
            with gzip.GzipFile(fileobj=r) as g:
                for i, line in enumerate(io.TextIOWrapper(g, errors="replace")):
                    if i % 4 == 1:
                        seqs.append(line.strip().upper()[:1_000_000])   # the D92 cap
                    if i >= 4 * PROBE - 1:
                        break
    except Exception:
        return None
    s = "".join(seqs)
    if not s:
        return None
    return [s.count(b) / len(s) for b in "ACGT"] + [len(seqs)]

rows = list(csv.DictReader(open(sys.argv[1]), delimiter="\t"))
kept = 0
for r in rows:
    r.update(role="screened_out", a_share="", c_share="", g_share="", t_share="", probe_reads="")
    if kept >= PANEL_N + BACKUP_N or not r.get("fastq_ftp"):
        continue
    v = shares(r["fastq_ftp"].split(";")[0])
    if v is None:
        continue
    r.update(a_share=f"{v[0]:.4f}", c_share=f"{v[1]:.4f}", g_share=f"{v[2]:.4f}", t_share=f"{v[3]:.4f}", probe_reads=str(v[4]))
    if v[2] >= G_MIN:
        r["role"] = "panel" if kept < PANEL_N else "backup"
        kept += 1
cols = list(rows[0].keys())
with open(sys.argv[2], "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(rows)
screened = sum(1 for r in rows if r["g_share"])
print(f"screened {screened} runs; G >= {G_MIN:.0%}: {kept} kept ({min(kept, PANEL_N)} panel, {max(0, kept - PANEL_N)} backup)")
