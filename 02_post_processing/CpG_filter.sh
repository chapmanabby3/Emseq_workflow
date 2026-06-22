#!/bin/bash
#
# ------------------------------------------------------------
# CpG_filter.sh
#
# Filters Bismark CX reports to retain only CpG methylation calls.
#
# Input:
#   *.CX_report.txt.gz
#
# Output:
#   *.CpG.txt
#
# Before running:
#   - Update the input and output directories.
#   - Update the sample list.
#   - Set the correct SLURM array size.
# ------------------------------------------------------------

#SBATCH --job-name=CpG_filter
#SBATCH --account=<account_name>

# Update log file locations as needed
#SBATCH --output=/path/to/logs/CpG_filter_%A_%a.out
#SBATCH --error=/path/to/logs/CpG_filter_%A_%a.err

# Set this to the number of samples listed below.
#SBATCH --array=1-<number_of_samples>

#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=16G

# Directory containing Bismark CX reports
indir="/path/to/methylation/results"

# Directory for filtered CpG output
outdir="/path/to/output_directory"

mkdir -p "$outdir"

# List sample IDs exactly as they appear in the Bismark output filenames.
SAMPLES=(
    "sample1_treatment_rep1"
    "sample2_control_rep1"
    "sample3_treatment_rep2"
)

echo "Starting CpG site filtering for selected samples..."

# SLURM array is 1-based; Bash arrays are 0-based.
sample_index=$((SLURM_ARRAY_TASK_ID - 1))
sample="${SAMPLES[$sample_index]}"

if [[ -z "$sample" ]]; then
    echo "Error: No sample mapped for SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID"
    exit 1
fi

# Construct the full path to the input file
file="$indir/${sample}_deduplicated.deduplicated.CX_report.txt.gz"

# Check if the file exists before processing
if [[ -f "$file" ]]; then

    # Get base name for output
    base="${sample}_deduplicated.deduplicated"

    echo "Processing $base..."

    zcat "$file" | awk '$6 == "CG"' > "$outdir/${base}.CpG.txt"

else

    echo "Error: Input file not found for sample $sample: $file"
    exit 1

fi

echo "Finished filtering all specified CpG sites."
echo "Output written to:"
echo "$outdir"
