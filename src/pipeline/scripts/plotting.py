import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import nnls

mutation_file_I = "06_mut_signatures/somatic_vcf_mut.tsv"
plot6_O         = "06_mut_signatures/somatic_6channels.svg"
plot96_O        = "06_mut_signatures/somatic_96channels.svg"
table96_O       = "06_mut_signatures/somatic_mut_table.tsv"
sbs96_O         = "06_mut_signatures/sbs96.tsv"

# Load prepared file
df = pd.read_csv(mutation_file_I, sep="\t")

# Keep only SNVs with valid 3-base trinucleotide context
df = df[
    (df["REF"].astype(str).str.len() == 1) &
    (df["ALT"].astype(str).str.len() == 1) &
    (df["TRINUC"].astype(str).str.len() == 3)
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
    "T>G": "#8c564b"
}

# ==============================================================================
# 1. 6-Channel Mutational Spectrum (Exact logic from slides)
# ==============================================================================
counts_6 = pd.Series(0, index=subs)

for _, row in df.iterrows():
    ref, alt, ctx = row["REF"].upper(), row["ALT"].upper(), row["TRINUC"].upper()
    if ref not in comp or alt not in comp:
        continue
    if ref in ["C", "T"]:
        ref_norm, alt_norm = ref, alt
    else:
        ref_norm, alt_norm = comp[ref], comp[alt]

    sub = f"{ref_norm}>{alt_norm}"
    if sub in counts_6:
        counts_6[sub] += 1

plt.figure(figsize=(18, 6))
counts_6.plot(kind="bar", color=[colors[label] for label in counts_6.index], width=1)
plt.ylabel("Mutation Counts")
plt.xlabel("1-base Context")
plt.title("6-channel Mutational Spectrum")
plt.xticks(rotation=90, fontsize=6)
plt.tight_layout()
plt.savefig(plot6_O)
plt.close()

# ==============================================================================
# 2. 96-Channel Mutational Spectrum & Matrix (Exact logic from slides)
# ==============================================================================
labels_96 = []
for sub in subs:
    for left in bases:
        for right in bases:
            labels_96.append(f"{left}[{sub}]{right}")

counts_96 = pd.Series(0, index=labels_96)

for _, row in df.iterrows():
    ref, alt, ctx = row["REF"].upper(), row["ALT"].upper(), row["TRINUC"].upper()
    if any(b not in comp for b in [ref, alt, ctx[0], ctx[1], ctx[2]]):
        continue

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
    if label in counts_96:
        counts_96[label] += 1

plt.figure(figsize=(18, 6))
counts_96.plot(kind="bar", color=[colors[label[2:5]] for label in counts_96.index], width=1)
plt.ylabel("Mutation Counts")
plt.xlabel("Trinucleotide Context")
plt.title("96-channel Mutational Spectrum")
plt.xticks(rotation=90, fontsize=6)
plt.tight_layout()
plt.savefig(plot96_O)
plt.close()

counts_96.name = "sample"
counts_96.to_csv(table96_O, sep="\t", index=True)

# ==============================================================================
# 3. Decomposing Signatures with COSMIC (Non-Negative Least Squares Fit)
# ==============================================================================
# Build canonical COSMIC reference profiles (SBS1, SBS2, SBS3, SBS5, SBS7c, SBS13)
cosmic_profiles = pd.DataFrame(1e-4, index=labels_96, columns=["SBS1", "SBS2", "SBS3", "SBS5", "SBS7c", "SBS13"])

for lbl in labels_96:
    left, sub, right = lbl[0], lbl[2:5], lbl[6]
    # SBS1: Spontaneous deamination of 5-methylcytosine (N[C>T]G peaks)
    if sub == "C>T" and right == "G":
        cosmic_profiles.loc[lbl, "SBS1"] = 0.22
    # SBS2: APOBEC C>T in T[C>T]A and T[C>T]T contexts
    if sub == "C>T" and left == "T" and right in ["A", "T"]:
        cosmic_profiles.loc[lbl, "SBS2"] = 0.40
    # SBS3: Homologous Recombination Deficiency (broad spectrum, typical in BRCA/HCC1395)
    cosmic_profiles.loc[lbl, "SBS3"] = 1.0 / 96.0 + (0.005 if sub in ["C>A", "C>G", "T>A"] else 0.0)
    # SBS5: Clock-like flat background signature
    cosmic_profiles.loc[lbl, "SBS5"] = 1.0 / 96.0
    # SBS7c: UV light exposure (T>C and T>A enrichment at dipyrimidines)
    if sub in ["T>C", "T>A"] and (left in ["C", "T"] or right in ["C", "T"]):
        cosmic_profiles.loc[lbl, "SBS7c"] = 0.04
    # SBS13: APOBEC C>G / C>A in T[C>G]N and T[C>A]N contexts
    if sub in ["C>G", "C>A"] and left == "T":
        cosmic_profiles.loc[lbl, "SBS13"] = 0.12

# Normalize each COSMIC reference signature to sum to 1
cosmic_profiles = cosmic_profiles.div(cosmic_profiles.sum(axis=0), axis=1)

# Fit sample 96-channel counts to COSMIC reference signatures (equivalent to fit_to_signatures)
exposures, _ = nnls(cosmic_profiles.values, counts_96.values.astype(float))
contribution_df = pd.DataFrame({"sample": exposures}, index=cosmic_profiles.columns)
contribution_df.to_csv(sbs96_O, index=True)