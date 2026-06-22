# Snakemake EM-seq Pipeline

This is a Snakemake pipeline implementation for processing EM-seq data on Saga.
It handles QC, trimming (manual Cutadapt + Trim Galore), alignment (Bismark),
deduplication, and methylation extraction.

## Directory Structure

```
snakemake_test/
├── Snakefile           # Workflow definition
├── config.yaml         # Configuration (sample list, paths, container images)
├── read_sequences.csv  # CSV file containing sample-specific Cutadapt sequences
├── run_snakemake.sh    # SLURM submission script
├── containers/         # Singularity sif images
├── genome/             # Reference genome and indices
├── data/               # Input FASTQ files (symlinked or placed here)
└── results/            # Output directory
```

## Configuration

### 1. Samples & Adapters (`read_sequences.csv`)
This pipeline uses a CSV file to manage sample-specific adapter sequences for the manual `cutadapt` step.
The file `read_sequences.csv` must contain the following columns:
*   `sample_id`: The unique identifier for the sample (must match the beginning of the filename, e.g., `11` for `11_male_...`).
*   `forward_read`: The adapter sequence for the forward read (passed to `-g`).
*   `reverse_read`: The adapter sequence for the reverse read (passed to `-G`).

**Example:**
```csv
sample_id,forward_read,reverse_read
1,CGAATACG,CGTATTCG
11,TCCTCATG,CATGAGGA
```

### 2. Main Config (`config.yaml`)
*   **samples**: A list of full sample basenames (e.g., `1_male_gonad_high_rep1`) to process.
*   **paths**: Directories for input, output, and genome.
*   **containers**: Paths to Singularity images.
*   **params**: Global settings for Trim Galore.

## Adding New Samples

1.  **Data**: Place paired-end FASTQ files in the `data/` directory.
    *   Naming convention: `{sample_name}_R1.fastq.gz` and `{sample_name}_R2.fastq.gz`.
2.  **Adapters**: Add the sample's ID and specific adapter sequences to `read_sequences.csv`.
3.  **Config**: Add the full `{sample_name}` to the `samples` list in `config.yaml`.

## Running the Pipeline

To submit the workflow to the cluster (SLURM), use the provided wrapper script:

```bash
sbatch run_snakemake.sh
```

This script:
*   Loads the Snakemake module.
*   Runs Snakemake with Singularity support.
*   Requests 64GB memory resource constraint (to prevent Bismark from crashing the node).

## Outputs
All outputs are stored in `results/`, organized by step:
*   `qc_raw/`: FastQC on raw reads.
*   `trimmed/`: Output from Cutadapt and Trim Galore.
*   `qc_trimmed/`: FastQC on trimmed reads.
*   `aligned/`: Bismark alignment BAM files.
*   `deduplicated/`: Deduplicated BAMs and reports.
*   `methylation/`: Methylation extraction reports (M-bias, CX report, etc.).
*   `qc_multiqc/`: Aggregate MultiQC report.
