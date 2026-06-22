#!/bin/bash

set -euo pipefail

# Load the module (adjust version if needed)
module reset
module load snakemake/8.4.2-foss-2023a

echo "Running Snakemake in Cluster Mode..."

# Create log directory for Slurm output files
mkdir -p logs/slurm

# --cluster "sbatch ..." tells Snakemake to submit each rule as a Slurm job
# --jobs 50 limits the max concurrent jobs to 50
# --default-resources sets defaults for rules that don't have them defined in Snakefile
snakemake  --unlock --cores 1
snakemake --use-singularity \
    --singularity-args "--cleanenv -B /cluster" \
    --jobs 50 \
    --executor cluster-generic \
    --cluster-generic-submit-cmd "sbatch --account=nn12038k --parsable --output=logs/slurm/slurm-%j.out --error=logs/slurm/slurm-%j.err --time=96:00:00 --cpus-per-task={threads} --mem={resources.mem_mb}M" \
    --default-resources mem_mb=8192 threads=1 \
    --printshellcmds \
    --rerun-incomplete \
    --latency-wait 120
