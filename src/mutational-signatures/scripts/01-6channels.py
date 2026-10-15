import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import argparse
from pathlib import Path
import os

# get arguments from shell call
parser = argparse.ArgumentParser(description="Process input data.")

parser.add_argument("--input", required=True, help="Path to input file")
parser.add_argument("--output", required=True, help="Output directory")

args = parser.parse_args()

input_path = Path(args.input)
output_dir = Path(args.output)
output_dir.mkdir(parents=True, exist_ok=True)

print(f"Input: {input_path}")
print(f"Output: {output_dir}")

filename_pic_O = output_dir / f"{input_path.stem}_6channels.svg"

# Read either a prepared TSV or a VCF file.
if input_path.name.lower().endswith((".vcf", ".vcf.gz")):
    vcf_columns = [
        "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
        "FORMAT", "SAMPLE",
    ]
    df = pd.read_csv(
        input_path,
        sep="\t",
        comment="#",
        header=None,
        names=vcf_columns,
    )
else:
    df = pd.read_csv(input_path, sep="\t")

# Keep only SNVs
df = df[
    df["REF"].str.len().eq(1, fill_value=False)
    & df["ALT"].str.len().eq(1, fill_value=False)
]

# Define the 6 substitution classes
subs = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
bases = ["A", "C", "G", "T"]

# Build the 96-channel labels
labels = []
for sub in subs:
    labels.append(f"{sub}")

counts = pd.Series(0, index=labels)

# Complement map for normalization
comp = {"A":"T", "T":"A", "C":"G", "G":"C"}

# Iterate mutations
for _, row in df.iterrows():
    ref, alt = row["REF"].upper(), row["ALT"].upper()

    # Normalize to pyrimidine convention
    if ref in ["C", "T"]:
        ref_norm, alt_norm = ref, alt
    else:
        ref_norm, alt_norm = comp[ref], comp[alt]

    sub = f"{ref_norm}>{alt_norm}"
    label = f"{sub}"

    if label in counts:
        counts[label] += 1

# Plotting colors by substitution type
colors = {
    "C>A": "#1f77b4",
    "C>G": "#ff7f0e",
    "C>T": "#2ca02c",
    "T>A": "#d62728",
    "T>C": "#9467bd",
    "T>G": "#8c564b"
}

plt.figure(figsize=(18,6))
counts.plot(kind="bar", color=[colors[label] for label in counts.index], width=1)

plt.ylabel("Mutation Counts")
plt.xlabel("1-base Context")
plt.title("6-channel Mutational Spectrum")
plt.xticks(rotation=90, fontsize=6)
plt.tight_layout()

plt.savefig(filename_pic_O)