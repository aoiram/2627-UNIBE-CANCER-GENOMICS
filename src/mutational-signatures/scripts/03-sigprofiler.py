#!/usr/bin/env python3
"""Run SigProfiler on a single-sample VCF: SBS96 matrix generation + COSMIC fitting.

Usage:
    03-sigprofiler.py --input sample.vcf --output out_dir --genome GRCh38 [--threads 4]

Results are written under <output>/ (exposures in
Assignment_Solution/Activities/Assignment_Solution_Activities.txt).
"""
import argparse
import os
import shutil

from SigProfilerAssignment import Analyzer as Analyze


p = argparse.ArgumentParser()
p.add_argument("--input", required=True, help="Single-sample VCF")
p.add_argument("--output", required=True, help="Output directory for this sample")
p.add_argument("--genome", default="GRCh38", help="Genome build (GRCh37, GRCh38, mm10...)")
p.add_argument("--cosmic_version", type=float, default=3.4)
p.add_argument("--threads", type=int, default=1)
args = p.parse_args()

out = os.path.abspath(args.output)
# SigProfiler expects a folder of VCFs, so stage the sample in its own folder
vcf_dir = os.path.join(out, "input_vcf")
os.makedirs(vcf_dir, exist_ok=True)
shutil.copy(args.input, os.path.join(vcf_dir, os.path.basename(args.input)))

Analyze.cosmic_fit(
    samples=vcf_dir,
    output=out,
    input_type="vcf",
    context_type="96",
    genome_build=args.genome,
    cosmic_version=args.cosmic_version,
    exome=False,
    make_plots=True,
    sample_reconstruction_plots="none",  # "pdf" needs to write COSMIC templates inside the read-only container
    verbose=False,
)