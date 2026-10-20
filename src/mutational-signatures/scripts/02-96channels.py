import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pandas as pd


parser = argparse.ArgumentParser(description="Create a 96-channel mutational spectrum.")
parser.add_argument("--input", required=True, help="Path to an input VCF file")
parser.add_argument("--output", required=True, help="Output directory")
parser.add_argument("--reference", required=True, help="Indexed reference FASTA file")
args = parser.parse_args()

input_path = Path(args.input)
output_dir = Path(args.output)
reference_path = Path(args.reference)
output_dir.mkdir(parents=True, exist_ok=True)

if not reference_path.with_suffix(reference_path.suffix + ".fai").is_file():
    raise FileNotFoundError(
        f"Reference index not found: {reference_path}.fai. "
        "Create it with samtools faidx before running the workflow."
    )


class IndexedFasta:
    def __init__(self, fasta_path):
        self.fasta_path = fasta_path
        self.contigs = {}
        with fasta_path.with_suffix(fasta_path.suffix + ".fai").open() as index:
            for line in index:
                name, length, offset, line_bases, line_width = line.rstrip().split("\t")
                self.contigs[name] = (
                    int(length),
                    int(offset),
                    int(line_bases),
                    int(line_width),
                )

    def sequence(self, contig, start, end):
        if contig not in self.contigs:
            return None

        length, offset, line_bases, line_width = self.contigs[contig]
        if start < 0 or end > length or start >= end:
            return None

        sequence = bytearray()
        with self.fasta_path.open("rb") as fasta:
            position = start
            while position < end:
                line_number, column = divmod(position, line_bases)
                read_length = min(end - position, line_bases - column)
                fasta.seek(offset + line_number * line_width + column)
                sequence.extend(fasta.read(read_length))
                position += read_length
        return sequence.decode("ascii").upper()


reference = IndexedFasta(reference_path)
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

df["POS"] = pd.to_numeric(df["POS"], errors="coerce")
df = df.dropna(subset=["POS"])
df["POS"] = df["POS"].astype(int)
df["REF"] = df["REF"].astype(str).str.upper()
df["ALT"] = df["ALT"].astype(str).str.split(",").str[0].str.upper()
df = df[
    df["REF"].str.len().eq(1)
    & df["ALT"].str.len().eq(1)
    & df["REF"].isin(["A", "C", "G", "T"])
    & df["ALT"].isin(["A", "C", "G", "T"])
]

subs = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
bases = ["A", "C", "G", "T"]
labels = [
    f"{left}[{sub}]{right}"
    for sub in subs
    for left in bases
    for right in bases
]
counts = pd.Series(0, index=labels, dtype=int)
comp = {"A": "T", "T": "A", "C": "G", "G": "C"}

for _, row in df.iterrows():
    ref, alt, position = row["REF"], row["ALT"], row["POS"]
    context = reference.sequence(row["CHROM"], position - 2, position + 1)
    if context is None or len(context) != 3 or context[1] != ref:
        continue

    if ref in ["C", "T"]:
        ref_norm, alt_norm = ref, alt
        left, right = context[0], context[2]
    else:
        ref_norm, alt_norm = comp[ref], comp[alt]
        left, right = comp[context[2]], comp[context[0]]

    label = f"{left}[{ref_norm}>{alt_norm}]{right}"
    if label in counts:
        counts[label] += 1

colors = {
    "C>A": "#1f77b4",
    "C>G": "#ff7f0e",
    "C>T": "#2ca02c",
    "T>A": "#d62728",
    "T>C": "#9467bd",
    "T>G": "#8c564b",
}

plt.figure(figsize=(18, 6))
counts.plot(
    kind="bar",
    color=[colors[label[2:5]] for label in counts.index],
    width=1,
)
plt.ylabel("Mutation Counts")
plt.xlabel("Trinucleotide Context")
plt.title("96-channel Mutational Spectrum")
plt.xticks(rotation=90, fontsize=6)
plt.tight_layout()

plot_path = output_dir / f"{input_path.stem}_96channels.svg"
table_path = output_dir / f"{input_path.stem}_table.tsv"
plt.savefig(plot_path)
counts.name = input_path.stem
counts.to_csv(table_path, sep="\t", index=True)
