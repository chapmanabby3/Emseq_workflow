```bash
#!/bin/bash
#
# ------------------------------------------------------------
# zipper.sh
#
# Compresses all uncompressed *.CpG.txt files in a directory
# using gzip. Existing compressed (.gz) files are skipped.
#
# Before running:
#   - Update the SLURM account.
#   - Update the log file locations if desired.
#   - Set the target directory containing the CpG files.
# ------------------------------------------------------------

#SBATCH --job-name=zipper
#SBATCH --account=your_account

# Update log file locations as needed
#SBATCH --output=/path/to/logs/zipper_%A.out
#SBATCH --error=/path/to/logs/zipper_%A.err

#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=4

# Directory containing the uncompressed *.CpG.txt files
target_dir="/path/to/CpG_directory"

echo "Starting compression of uncompressed .CpG.txt files in:"
echo "$target_dir"

# Navigate to the target directory
cd "$target_dir" || {
    echo "Error: Could not change to directory $target_dir."
    exit 1
}

# Compress each uncompressed CpG file
for file in *.CpG.txt; do

    # Skip if no matching files are found
    [[ -e "$file" ]] || {
        echo "No .CpG.txt files found."
        break
    }

    if [[ ! -f "${file}.gz" ]]; then

        echo "Compressing $file..."

        gzip "$file" &

    else

        echo "Skipping $file (already compressed)."

    fi

done

# Wait for all background gzip jobs to finish
wait

echo "Compression complete."
```
