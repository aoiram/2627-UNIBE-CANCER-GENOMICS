import sys
import random
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

positions_file = sys.argv[1] if len(sys.argv) > 1 else "06_mut_signatures/file_positions.txt"
output_vcf     = sys.argv[2] if len(sys.argv) > 2 else "06_mut_signatures/simulated_mutations.vcf"
output_plot    = sys.argv[3] if len(sys.argv) > 3 else "06_mut_signatures/simulated_96channels.svg"
N_MUTATIONS    = 500

random.seed(42)
comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
subs = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
bases = ["A", "C", "G", "T"]

# 1. LOAD POSITIONS WITH TRINUCLEOTIDES
context_pools = {}
with open(positions_file) as f:
    lines = [line.strip() for line in f if line.strip()]

for i in range(0, len(lines) - 1, 2):
    header, seq = lines[i], lines[i + 1].upper()
    if not header.startswith(">") or len(seq) != 3 or any(b not in comp for b in seq):
        continue
    chrom, coords = header[1:].split(":")
    start, _ = map(int, coords.split("-"))
    pos = start + 1  # Center base coordinate
    ref = seq[1]

    if ref in ["C", "T"]:
        norm_ctx = seq
        strand = "+"
    else:
        norm_ctx = comp[seq[2]] + comp[seq[1]] + comp[seq[0]]
        strand = "-"

    context_pools.setdefault(norm_ctx, []).append((chrom, pos, ref, seq, strand))

# 2. BUILD & NORMALIZE COSMIC PROFILES (SBS1 and SBS7c)
labels_96 = [f"{l}[{s}]{r}" for s in subs for l in bases for r in bases]
sbs1 = pd.Series(0.001, index=labels_96)
sbs7c = pd.Series(0.001, index=labels_96)

for lbl in labels_96:
    left, sub, right = lbl[0], lbl[2:5], lbl[6]
    if sub == "C>T" and right == "G":
        sbs1[lbl] = 0.22
    if sub in ["T>C", "T>A"] and (left in ["C", "T"] or right in ["C", "T"]):
        sbs7c[lbl] = 0.05

sbs1 /= sbs1.sum()
sbs7c /= sbs7c.sum()

# 3. MERGE SIGNATURES ACCORDING TO MIX_RATIO (70% SBS1 + 30% SBS7c)
mixed_profile = 0.70 * sbs1 + 0.30 * sbs7c
mixed_profile /= mixed_profile.sum()

# 4 & 5. SIMULATE MUTATIONS
simulated_variants = []
sim_counts_96 = pd.Series(0, index=labels_96)

available_labels = [
    lbl for lbl in labels_96
    if (lbl[0] + lbl[2] + lbl[6]) in context_pools and len(context_pools[lbl[0] + lbl[2] + lbl[6]]) > 0
]
probs = mixed_profile[available_labels]
probs /= probs.sum()

for _ in range(N_MUTATIONS):
    chosen_lbl = random.choices(available_labels, weights=probs.values, k=1)[0]
    norm_ctx = chosen_lbl[0] + chosen_lbl[2] + chosen_lbl[6]
    alt_norm = chosen_lbl[4]

    chrom, pos, ref, raw_seq, strand = random.choice(context_pools[norm_ctx])
    alt = alt_norm if strand == "+" else comp[alt_norm]
    simulated_variants.append((chrom, pos, ref, alt, chosen_lbl))
    sim_counts_96[chosen_lbl] += 1

# 6. WRITE OUTPUT VCF
with open(output_vcf, "w") as out:
    out.write("##fileformat=VCFv4.2\n")
    out.write("##source=MutSigSimulator_70pct_SBS1_30pct_SBS7c\n")
    out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    for chrom, pos, ref, alt, lbl in simulated_variants:
        out.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\tSIM_CHANNEL={lbl}\n")

# 7. PLOT SIMULATED 96-CHANNEL PROFILE
colors = {
    "C>A": "#1f77b4", "C>G": "#ff7f0e", "C>T": "#2ca02c",
    "T>A": "#d62728", "T>C": "#9467bd", "T>G": "#8c564b"
}
plt.figure(figsize=(18, 6))
sim_counts_96.plot(kind="bar", color=[colors[label[2:5]] for label in sim_counts_96.index], width=1)
plt.ylabel("Simulated Mutation Counts")
plt.xlabel("Trinucleotide Context")
plt.title("Simulated 96-channel Spectrum (70% SBS1 + 30% SBS7c)")
plt.xticks(rotation=90, fontsize=6)
plt.tight_layout()
plt.savefig(output_plot)
plt.close()