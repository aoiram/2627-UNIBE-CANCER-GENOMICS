import sys
import os
import random

bed_file    = sys.argv[1] if len(sys.argv) > 1 else "exome_regions.bed"
output_file = sys.argv[2] if len(sys.argv) > 2 else "06_mut_signatures/random_positions.txt"
N           = int(sys.argv[3]) if len(sys.argv) > 3 else 1000

random.seed(42)

# Load exonic regions (fallback to interval_list if .bed is not directly present)
regions = []
if os.path.exists(bed_file):
    with open(bed_file) as f:
        for line in f:
            if line.startswith("#") or line.startswith("@") or not line.strip():
                continue
            parts = line.strip().split()
            regions.append((parts[0], int(parts[1]), int(parts[2])))
elif os.path.exists(bed_file + ".interval_list"):
    with open(bed_file + ".interval_list") as f:
        for line in f:
            if line.startswith("@") or not line.strip():
                continue
            parts = line.strip().split()
            regions.append((parts[0], int(parts[1]), int(parts[2])))
else:
    regions = [("chr6", 60000, 170000000), ("chr17", 60000, 83000000)]

# Randomly sample positions
coords = []
while len(coords) < N:
    chrom, start, end = random.choice(regions)
    if end - 1 > start + 1:
        pos = random.randint(start + 1, end - 1)
        coords.append((chrom, pos))

# Write output in exact format from slides
with open(output_file, "w") as out:
    for chrom, pos in coords:
        out.write(f">{chrom}:{pos}-{pos}\n")
        out.write("N\n")