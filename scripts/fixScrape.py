#!/usr/bin/env python3
import pandas as pd
import glob
import os

# 1) Load the “gold-standard” header from the FANCA file
ref_path = "data/scraped/FANCA_variants_full.tsv"
ref_cols = pd.read_csv(ref_path, sep="\t", nrows=0).columns.tolist()

# 2) Prepare input and output dirs
in_dir  = "data/scraped"
out_dir = "scraped_adj"
os.makedirs(out_dir, exist_ok=True)

# 3) Process each TSV
for path in glob.glob(f"{in_dir}/*.tsv"):
    fname = os.path.basename(path)
    # Skip the reference itself (it already matches!)
    if os.path.abspath(path) == os.path.abspath(ref_path):
        print(f"Skipping reference file {fname}")
        continue

    # Read file (preserve empties as "")
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)

    # If columns exactly match, do nothing
    if df.columns.tolist() == ref_cols:
        print(f"{fname}: already matches, skipping")
        continue

    # Otherwise: add any missing, drop extras, reorder
    for col in ref_cols:
        if col not in df.columns:
            df[col] = ""
    df = df[ref_cols]

    # Write adjusted file
    out_path = os.path.join(out_dir, fname)
    df.to_csv(out_path, sep="\t", index=False)
    print(f"{fname}: written adjusted → {out_path}")
