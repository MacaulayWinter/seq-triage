#!/usr/bin/env python3
"""Validate each samplesheet recommend.R wrote against that pipeline's OWN
samplesheet schema (assets/schema_input.json), fetched at the commit pinned in
assets/downstream_pipelines.tsv (CLAUDE.md D66). A format mismatch fails here,
in CI, instead of on the user's HPC.

Empty cells are left out of each row, as nf-schema does. File existence is not
checked: the paths are synthetic.

Usage: python3 tests/recommend/check_schemas.py <recommendations_dir> [assets/downstream_pipelines.tsv]
"""
import csv, json, sys, urllib.request
from pathlib import Path
from jsonschema import Draft202012Validator

out = Path(sys.argv[1])
table = sys.argv[2] if len(sys.argv) > 2 else "assets/downstream_pipelines.tsv"
rules = {}
with open(table) as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        rules[row["rule"]] = (row["pipeline"], row["commit"])

failures, checked = 0, 0
for rule, (pipeline, commit) in sorted(rules.items()):
    sheet = out / rule / "samplesheet.csv"
    if not sheet.exists():
        continue
    url = f"https://raw.githubusercontent.com/nf-core/{pipeline}/{commit}/assets/schema_input.json"
    schema = json.load(urllib.request.urlopen(url, timeout=60))
    with open(sheet) as fh:
        rows = [{k: v for k, v in r.items() if v != ""} for r in csv.DictReader(fh)]
    errors = list(Draft202012Validator(schema).iter_errors(rows))
    checked += 1
    if errors:
        failures += 1
        print(f"FAIL {rule} (nf-core/{pipeline}@{commit[:7]}): {len(errors)} schema error(s)")
        for e in errors[:5]:
            print(f"     {list(e.absolute_path)}: {e.message}")
    else:
        print(f"PASS {rule} (nf-core/{pipeline}@{commit[:7]}): {len(rows)} row(s) valid")
print(f"\n{checked - failures} / {checked} samplesheets valid against their pipeline's schema")
sys.exit(1 if failures or checked == 0 else 0)
