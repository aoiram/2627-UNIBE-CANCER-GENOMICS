import matplotlib.pyplot as plt
import pandas as pd

# 1. Load and filter variant data
df = pd.read_csv("06_mut_signatures/somatic_vcf_mut.tsv", sep="\t")
df = df[
    (df["REF"].str.len() == 1)
    & (df["ALT"].str.len() == 1)
    & (df["TRINUC"].str.len() == 3)
]

subs = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
bases = ["A", "C", "G", "T"]
comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
colors = {
    "C>A": "#1f77b4",
    "C>G": "#ff7f0e",
    "C>T": "#2ca02c",
    "T>A": "#d62728",
    "T>C": "#9467bd",
    "T>G": "#8c564b",
}

# 2. Build 6-channel and 96-channel series
counts_6 = pd.Series(0, index=subs)
labels_96 = [
    f"{left}[{sub}]{right}"
    for sub in subs
    for left in bases
    for right in bases
]
counts_96 = pd.Series(0, index=labels_96)

for _, row in df.iterrows():
    ref, alt, ctx = row["REF"].upper(), row["ALT"].upper(), row["TRINUC"].upper()
    if ctx[1] != ref:
        continue
    if ref in ["C", "T"]:
        ref_norm, alt_norm, left, right = ref, alt, ctx[0], ctx[2]
    else:
        ref_norm, alt_norm = comp[ref], comp[alt]
        left, right = comp[ctx[2]], comp[ctx[0]]

    sub = f"{ref_norm}>{alt_norm}"
    label_96 = f"{left}[{sub}]{right}"

    if sub in counts_6:
        counts_6[sub] += 1
    if label_96 in counts_96:
        counts_96[label_96] += 1

# 3. Save 6-channel plot
plt.figure(figsize=(10, 6))
counts_6.plot(kind="bar", color=[colors[k] for k in counts_6.index], width=0.8)
plt.ylabel("Mutation Counts")
plt.xlabel("Substitution Type")
plt.title("6-channel Mutational Spectrum")
plt.xticks(rotation=0)
plt.tight_layout()
plt.savefig("06_mut_signatures/somatic_6channels.svg")
plt.close()

# 4. Save 96-channel plot
plt.figure(figsize=(18, 6))
counts_96.plot(
    kind="bar", color=[colors[k[2:5]] for k in counts_96.index], width=1.0
)
plt.ylabel("Mutation Counts")
plt.xlabel("Trinucleotide Context")
plt.title("96-channel Mutational Spectrum")
plt.xticks(rotation=90, fontsize=6)
plt.tight_layout()
plt.savefig("06_mut_signatures/somatic_96channels.svg")
plt.close()

# 5. Save 96-channel table for MutationalPatterns in R
counts_96.name = "tumor_sample"
counts_96.to_csv("06_mut_signatures/somatic_mut_table.tsv", sep="\t", index=True)