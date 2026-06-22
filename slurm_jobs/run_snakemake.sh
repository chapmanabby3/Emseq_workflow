#!/bin/bash
set -euo pipefail

############################################################
# Snakemake cluster submission script (generic template)
############################################################

# Load Snakemake module
module reset
module load snakemake/8.4.2-foss-2023a

echo "Starting Snakemake workflow (cluster mode)..."

# Create local log directory for Slurm outputs
mkdir -p logs/slurm

############################################################
# Initial unlock (safe to keep)
############################################################
snakemake --unlock

############################################################
# Run workflow on cluster
############################################################
snakemake \
    --use-singularity \
    --singularity-args "--cleanenv -B /cluster" \
    --jobs 50 \
    --executor cluster-generic \
    --cluster-generic-submit-cmd "sbatch --account=ACCOUNT_ID --parsable --output=logs/slurm/%x_%j.out --error=logs/slurm/%x_%j.err" \
    --default-resources mem_mb=8192 threads=1 \
    --printshellcmds \
    --rerun-incomplete
