#!/bin/bash
#
#=========================================================
# SLURM job: DSS pairwise methylation analysis
#
# This job runs the DSS differential methylation script.
#
# Before running:
#   - Update R script name
#   - Update log paths
#   - Ensure BSseq object exists
#=========================================================

#SBATCH --job-name=DSS_pairwise
#SBATCH --account=your_account
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=150G
#SBATCH --time=24:00:00

#-----------------------------
# Log files
#-----------------------------
#SBATCH --output=/path/to/logs/dss_%j.out
#SBATCH --error=/path/to/logs/dss_%j.err

set -euo pipefail

#-----------------------------
# Create log directory
#-----------------------------
mkdir -p /path/to/logs

#-----------------------------
# Load environment
#-----------------------------
module load R/4.5.2-gfbf-2025b

#-----------------------------
# Run DSS analysis
#-----------------------------
Rscript DSS_pairwise.R
