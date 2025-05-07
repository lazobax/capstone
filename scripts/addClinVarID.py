import os
import pandas as pd
import sys

input_dir = "data/ClinVar_vcf/"  # ← change this to your directory

file = sys.argv[1]

if file.endswith(".vcf"):
    file_path = os.path.join(input_dir, file)
    df = pd.read_csv(file_path, sep="\t", header=None)  # Keep header

    if df.shape[1] != 7:
        print(f"Skipping {file}: does not have 7 columns")
        sys.exit(0)

    gene = file.rsplit(".", 1)[0]
    ids = [f"CV.{gene}{str(i+1).zfill(5)}" for i in range(len(df))]

    df.insert(5, "ID", ids)  # Insert after the 5th column (0-indexed)

    # Define new header
    header = ["#CHROM", "POS", "REF", "ALT", "TYP", "ID", "CONSEQ", "REPORTED"]

    df.to_csv(file_path, sep="\t", header=header, index=False)
    print(f"✅ Updated: {file}")
