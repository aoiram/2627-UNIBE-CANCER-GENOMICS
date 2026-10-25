# Cancer Genomics: Somatic Variant Calling and Mutational Signatures

**Group 7**

- Mario Amos
- Isaac Ambrogetti
- Marcos Valdivia

This repository contains Group 7's cancer genomics course project. It implements
a Snakemake workflow for paired normal/tumor short-read data, from read quality
control through somatic variant calling, annotation, and mutational-signature
analysis.

The project also includes an in-silico mutational-signature simulation based on
exonic positions from the course reference resources and a mixture of the
COSMIC SBS1 and SBS7c profiles.

## Project status

The workflow is written is launched through Slurm.
The input data and supporting resources are provided through the course and
cluster directories referenced by the `Snakefile`; they are not stored in this
repository. The workflow therefore requires access to those paths, as well as
the configured annotation caches and Apptainer/Singularity images.

## Repository layout

```text
.
├── src/
│   ├── README.md
│   └── pipeline/
│       ├── Snakefile
│       ├── pipeline.slurm
│       └── scripts/
│           ├── 01_generate_random_coords.py
│           ├── 02_apply_cosmic_signatures.py
│           └── plotting.py
└── README.md
```

The `Snakefile` defines the authoritative locations for the paired
`normal`/`tumor` FASTQs, the reference genome, interval files, known-sites
resources, germline resource, panel of normals, ClinVar data, and VEP cache.
The repository contains the workflow code and scripts, not those input data
files.

## Workflow overview

The pipeline is organized into numbered output directories:

1. **Read QC and trimming (`01_cleaned_data`, `02_qc`)**
   - Runs FastQC on both read pairs.
   - Trims/adapts reads with fastp.
   - Aggregates alignment and preprocessing metrics with MultiQC.

2. **Alignment (`03_alignment`)**
   - Aligns cleaned reads with BWA-MEM.
   - Converts, sorts, and indexes BAM files with Samtools.

3. **BAM preprocessing (`04_preprocessing`)**
   - Marks duplicates with GATK.
   - Runs GATK base-quality score recalibration using dbSNP and known indels.
   - Produces recalibrated, indexed BAM files.

4. **Post-alignment QC (`02_qc`)**
   - Calculates Samtools flagstat metrics.
   - Collects hybrid-selection metrics.
   - Creates `multiqc_report.html`.

5. **Somatic variant calling (`05_variants`)**
   - Calls somatic variants with Mutect2 using the recalibrated tumor and normal
     BAMs.
   - Uses a germline resource and panel of normals.
   - Filters calls with FilterMutectCalls.
   - Retains `PASS` variants with bcftools.

6. **Variant annotation and interpretation (`05_variants`)**
   - Annotates variants with offline Ensembl VEP using the GRCh38 cache.
   - Produces a standard VEP annotation and a ClinVar-enriched annotation.
   - Filters for high/moderate-impact, damaging-prediction, or ClinVar
     pathogenic variants.
   - Writes a simple mutation-type summary.

7. **Observed mutational signatures (`06_mut_signatures`)**
   - Extracts trinucleotide contexts around passing SNVs.
   - Generates six-channel and 96-channel mutation spectra.
   - Fits the observed 96-channel counts to the signatures represented in
     `plotting.py` using non-negative least squares.

8. **In-silico signature simulation (`06_mut_signatures`)**
   - Samples 2,000 random positions from the configured exome regions.
   - Extracts their trinucleotide contexts.
   - Simulates 500 mutations from a 70% SBS1 / 30% SBS7c mixture.
   - Writes a simulated VCF and a simulated 96-channel spectrum.

## Expected outputs

The `all` rule in [`src/pipeline/Snakefile`](src/pipeline/Snakefile) declares
these final outputs:

```text
02_qc/multiqc_report.html
05_variants/chr6ch17_mutect2_annotated.txt
05_variants/chr6ch17_mutect2_annotated_filtered.txt
05_variants/chr6ch17_mutect2_annotated_with_clinVar.txt
05_variants/chr6ch17_mutect2_annotated_with_clinVar_filtered.txt
05_variants/summary.txt
06_mut_signatures/somatic_6channels.svg
06_mut_signatures/somatic_96channels.svg
06_mut_signatures/somatic_mut_table.tsv
06_mut_signatures/sbs96.tsv
06_mut_signatures/simulated_mutations.vcf
06_mut_signatures/simulated_96channels.svg
```

Intermediate files include cleaned FASTQs, sorted and recalibrated BAMs and
indexes, duplicate-marking metrics, recalibration tables, VCFs, and additional
QC reports.

## Requirements

The cluster launch script assumes:

- Slurm with the `pibu_el8` partition.
- Snakemake 6.6.1, loaded through the cluster module system.
- Apptainer (the script supplies a Singularity compatibility wrapper).
- The following tools available through the configured container images:
  FastQC, fastp, BWA, Samtools, GATK, MultiQC, bcftools, and VEP.
- Python 3 for plotting and
  signature fitting.
- The course reference and resource files listed near the top of the
  `Snakefile`, including the GRCh38 reference, exome intervals, dbSNP, known
  indels, germline resource, panel of normals, ClinVar, and VEP cache.


## Running on the IBU cluster

1. Confirm that the shared-resource and container paths at the top of
   [`src/pipeline/Snakefile`](src/pipeline/Snakefile) are available on the
   target cluster.
2. Confirm that the input FASTQs are available at the paths expected by the
   `fastqc` and `fastp` rules:

   ```text
   /data/courses/cancergenomics/VAR_CALLING/fastq/normal.R1.fq.gz
   /data/courses/cancergenomics/VAR_CALLING/fastq/normal.R2.fq.gz
   /data/courses/cancergenomics/VAR_CALLING/fastq/tumor.R1.fq.gz
   /data/courses/cancergenomics/VAR_CALLING/fastq/tumor.R2.fq.gz
   ```

3. Check the Slurm account, reservation, email address, partition, memory, and
   time limit in [`src/pipeline/pipeline.slurm`](src/pipeline/pipeline.slurm).
4. Submit the workflow , check the readme in `./src/`

The script creates the numbered output directories, unlocks the Snakemake working directory, and runs:

```bash
snakemake --cores 8 --use-singularity --rerun-incomplete
```

For a dry run or a local inspection of the dependency graph, run Snakemake
from `src/pipeline` with the appropriate cluster/module environment, for
example:

```bash
snakemake -n
snakemake --dag | dot -Tsvg > dag.svg
```

## Configuration notes and limitations

- **Input location:** All input data and reference resources are taken from
  the course and cluster directories defined in the `Snakefile`. No local
  `data/` directory is required.
- **Reference build:** The workflow is configured for the GRCh38 reference
  and VEP annotation cache supplied through the course resources. Reference
  and annotation resources must use compatible genome builds.
- **Exome intervals:** The workflow is configured around the course exome
  intervals and currently calls Mutect2 on those intervals. The output naming
  `chr6ch17` reflects the course data/resource setup and should be revised if
  the interval set changes.
- **Reproducibility:** The simulation scripts seed Python's random number
  generator with `42`. The observed signature fitting uses fixed, hand-built
  reference profiles in `plotting.py`, not a downloaded COSMIC profile matrix.