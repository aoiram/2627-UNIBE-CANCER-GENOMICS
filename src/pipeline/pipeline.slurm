#!/bin/bash
#SBATCH --job-name=class_assignment
#SBATCH --partition=pibu_el8
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=48G
#SBATCH --time=06:00:00
##SBATCH --account=class-407009-ws-2026 # Comment if out of reservation schedule
##SBATCH --reservation=class-407009-2026-10-22 # Comment if out of reservation schedule
#SBATCH --mail-user=mario.amos@students.unibe.ch
#SBATCH --mail-type=end

# Clean previous runs
find . -name "slurm-*.out" ! -name "slurm-${SLURM_JOB_ID}.out" -delete

sleep 2

module load snakemake/6.6.1-foss-2021a

# --- APPTAINER VERSION WORKAROUND ---
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
export PATH=~/.tmp_bin:$PATH
# ------------------------------------

# Pre-create sequentially numbered directories and scripts directory
mkdir -p 01_cleaned_data 02_qc 03_alignment 04_preprocessing 05_variants 06_mut_signatures scripts

# Ensure Python virtual environment has required packages
if [ ! -d "$HOME/.venvs/plotting" ]; then
    python3 -m venv ~/.venvs/plotting
    env PYTHONPATH="" ~/.venvs/plotting/bin/pip install --no-cache-dir --ignore-installed pandas matplotlib numpy scipy
else
    env PYTHONPATH="" ~/.venvs/plotting/bin/pip install --quiet pandas matplotlib numpy scipy
fi

env PYTHONPATH="" ~/.venvs/plotting/bin/python -c "import pandas, matplotlib, numpy, scipy; print('Venv ready!')"

# Unlock in case previous run was interrupted
snakemake --unlock

# Run Snakemake (keep --bind /data,/home so /data/courses and /data/users are accessible)
snakemake --cores 8 --default-resources "threads=8" \
          --use-singularity \
          --singularity-args "--bind /data,/home --writable-tmpfs --fakeroot --pwd $PWD" \
          --rerun-incomplete --show-failed-logs