# Cancer Genomics

## Group 7

- Mario Amos
- Isaac Ambrogetti
- Marcos Valdivia


# Mutational signature analysis

### Setup
Install reference genome and index it
```
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit
http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
gunzip hg19.fa.gz
```
Now index the fasta in a srun
```
srun --cpus-per-task=2 --mem-per-cpu=10GB --time=00:30:00 --partition=pshort_el8 \
    apptainer exec \
    --bind /data \
    /containers/apptainer/samtools-1.19.sif \
    samtools faidx \
    /data/users/iambrogetti/2526-unibe-cancer-genomics/data/MUT_SIGNATURES/hg19.fa
```

## Pipeline
Start by checking the `src/mutational-signatures/config_main.yaml` since it's there where the session information need to be change to secure a smooth run of the pipeline.

Then run the pipeline with this command:
```
python src/mutational-signature/main-pipeline.py
```

remember to load the snakemake environment (do not use conda env at least for the second command)
