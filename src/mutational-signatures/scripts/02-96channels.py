import pandas as pd
import matplotlib.pyplot as plt
import os

mutation_file_I = "MUTATIONS/DNA2406740.vcf_mut.tsv"
mutation_file_I = "MUTATIONS/DNA2406784.vcf_mut.tsv"
mutation_file_I = "MUTATIONS/DNA2408066.vcf_mut.tsv"
mutation_file_I = "MUTATIONS/DNA2412161.vcf_mut.tsv"
mutation_file_I = "MUTATIONS/DNA2505059.vcf_mut.tsv"
mutation_file_I = "MUTATIONS/DNA2505926.vcf_mut.tsv"
mutation_file_I = "MUTATIONS/DNA2506590.vcf_mut.tsv"
mutation_file_I = "MUTATIONS/DNA2508977.vcf_mut.tsv"
mutation_file_I = "MUTATIONS/DNA2510817.vcf_mut.tsv"

filename_pic_O = mutation_file_I.replace('.tsv','_96channels.svg')
filename_table_O = mutation_file_I.replace('.tsv','_table.tsv')

# Load prepared file
df = pd.read_csv(mutation_file_I, sep="\t")

# Keep only SNVs
df = df[(df["REF"].str.len() == 1) & (df["ALT"].str.len() == 1) & (df["TRINUC"].str.len() == 3)]

# Define the 6 substitution classes
subs = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
bases = ["A", "C", "G", "T"]

# Build the 96-channel labels
labels = []
for sub in subs:
    ref, _ = sub.split(">")
    for left in bases:
        for right in bases:
            labels.append(f"{left}[{sub}]{right}")

counts = pd.Series(0, index=labels)

# Complement map for normalization
comp = {"A":"T", "T":"A", "C":"G", "G":"C"}

# Iterate mutations
for _, row in df.iterrows():
    ref, alt, ctx = row["REF"].upper(), row["ALT"].upper(), row["TRINUC"].upper()

    # Center base in trinuc must match REF
    if ctx[1] != ref:
        continue

    # Normalize to pyrimidine convention
    if ref in ["C", "T"]:
        ref_norm, alt_norm, left, right = ref, alt, ctx[0], ctx[2]
    else:
        ref_norm, alt_norm = comp[ref], comp[alt]
        left, right = comp[ctx[2]], comp[ctx[0]]  # reverse complement

    sub = f"{ref_norm}>{alt_norm}"
    label = f"{left}[{sub}]{right}"

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
counts.plot(kind="bar", color=[colors[label[2:5]] for label in counts.index], width=1)

plt.ylabel("Mutation Counts")
plt.xlabel("Trinucleotide Context")
plt.title("96-channel Mutational Spectrum")
plt.xticks(rotation=90, fontsize=6)
plt.tight_layout()
#plt.show()

plt.savefig(filename_pic_O)

counts.name = 'sample'

counts.to_csv(filename_table_O,sep='\t',index=True)