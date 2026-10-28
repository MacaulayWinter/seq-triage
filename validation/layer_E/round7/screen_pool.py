#!/usr/bin/env python3
"""Layer E round 7 (CLAUDE.md D94): keep the G-dominated runs of a random pool.

WHAT IT DOES  Reads the pool panel (select_panel.R, its rows already in a
              seeded random order) and, in that order, streams the first 2,000
              reads of each run's first FASTQ (R1, the read detection checks)
              from ENA, measuring the share of A, C, G and T. The first 12 runs
              with G >= 50% become role "panel", the next 4 "backup" (used only
              if a download fails); every other run is "screened_out". No other
              choice is made, so no run is picked by hand.
WIDENED       (D94, pool 1 gave 7 of the 8 needed) Several pools, screened in
              the order given. A run whose shares are already recorded (an
              earlier screen's output) is not downloaded again, so pool 1's
              result stands as screened. New runs are fetched 12 at a time,
              with 2 retries; the choice still follows pool order.
WRITES        <out_panel.tsv>: the pools' rows and columns, the new roles, and
              columns a_share, c_share, g_share, t_share, probe_reads, pool.

Usage: python3 validation/layer_E/round7/screen_pool.py <out_panel.tsv> <pool.tsv> [<pool.tsv> ...]
"""
import csv, gzip, io, sys, urllib.request
from concurrent.futures import ThreadPoolExecutor

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

def shares_retry(url):
    for _ in range(3):
        v = shares(url)
        if v is not None:
            return v
    return None

rows = []
known = []
for k, f in enumerate(sys.argv[2:], start=1):
    rd = csv.DictReader(open(f), delimiter="\t")
    screened = "g_share" in (rd.fieldnames or [])   # an earlier screen's output, failures included
    for r in rd:
        r.setdefault("pool", str(k))
        rows.append(r); known.append(screened)
todo = [i for i, r in enumerate(rows) if not known[i] and r.get("fastq_ftp")]
with ThreadPoolExecutor(12) as ex:
    got = dict(zip(todo, ex.map(lambda i: shares_retry(rows[i]["fastq_ftp"].split(";")[0]), todo)))
kept = 0
for i, r in enumerate(rows):
    if not known[i]:
        v = got.get(i)
        r.update(a_share="", c_share="", g_share="", t_share="", probe_reads="")
        if v is not None:
            r.update(a_share=f"{v[0]:.4f}", c_share=f"{v[1]:.4f}", g_share=f"{v[2]:.4f}", t_share=f"{v[3]:.4f}", probe_reads=str(v[4]))
    r["role"] = "screened_out"
    if r["g_share"] and float(r["g_share"]) >= G_MIN and kept < PANEL_N + BACKUP_N:
        r["role"] = "panel" if kept < PANEL_N else "backup"
        r["stratum"] = "E7_ILM_TINY_G"   # one group over both pools, so a backup covers either
        kept += 1
cols = list(dict.fromkeys(c for r in rows for c in r))
with open(sys.argv[1], "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n", restval=""); w.writeheader(); w.writerows(rows)
screened = sum(1 for r in rows if r["g_share"])
print(f"screened {screened} runs; G >= {G_MIN:.0%}: {kept} kept ({min(kept, PANEL_N)} panel, {max(0, kept - PANEL_N)} backup)")
