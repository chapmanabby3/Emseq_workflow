#!/bin/bash
#SBATCH --job-name=annotate_pie_volcano
#SBATCH --account=your_account
#SBATCH --output=plots_annotate.out
#SBATCH --error=plots_annotate.err
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

module load R/4.5.2-gfbf-2025b

Rscript volc_pie_annotate.R
