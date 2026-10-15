#!/bin/bash
#SBATCH --job-name=class_assignment
#SBATCH --partition=pibu_el8
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mem=48G
#SBATCH --time=06:00:00
##SBATCH --account=class-407009-ws-2026 # Comment if out of reservation schedule
##SBATCH --reservation=class-407009-2026-10-08 # Comment if out of reservation schedule
#SBATCH --mail-user=mario.amos@students.unibe.ch
#SBATCH --mail-type=end

# Clean previous runs
find . -name "slurm-*.out" ! -name "slurm-${SLURM_JOB_ID}.out" -delete

sleep 3

module load snakemake/6.6.1-foss-2021a

# --- APPTAINER VERSION WORKAROUND ---
# Create a fake 'singularity' command that returns a version Snakemake 6 expects
mkdir -p ~/.tmp_bin
cat << 'EOF' > ~/.tmp_bin/singularity
#!/bin/bash
if [[ "$1" == "--version" || "$1" == "version" ]]; then
    echo "singularity version 3.8.0"
else
    apptainer "$@"
fi
EOF
chmod +x ~/.tmp_bin/singularity

# Prioritize our fake singularity in the system PATH
export PATH=~/.tmp_bin:$PATH
# ------------------------------------

# Pre-create the sequentially numbered directories to ensure clean execution
mkdir -p 01_cleaned_data 02_qc 03_alignment 04_preprocessing 05_variants

# Note : There can be some folder locking issues when a job get interrupted
#        (ig by the end of the course reservation schedule),
#        just add --unlock and rerun to unlock and then rerun normally
snakemake --cores 12 --default-resources "threads=12" \
          --use-singularity \
          --singularity-args "--bind /data,/containers,$PWD --pwd $PWD" \
          --rerun-incomplete