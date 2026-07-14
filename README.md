# EM-seq Bioinformatics Workflow (Snakemake + DSS + R Analysis)

This repository contains a complete workflow for processing EM-seq / bisulfite sequencing data, from raw reads to differential methylation analysis and visualization.

The pipeline is designed for execution on an HPC system (e.g. Saga) using Snakemake, Slurm, and R-based downstream analysis.
This pipeline was created for DNA methylation analysis of a non model teleost fish (Boreogadus saida).
---

## Workflow Overview

The workflow is organized into four main stages:

### 01_data_prep
Initial setup and data preparation:
- Downloading raw sequencing data
- Preparing reference genome (Bismark indexing)
- Preparing control genomes (e.g. lambda, pUC19)
- Organizing input structure for Snakemake

---

### 02_post_processing
Processing of methylation call outputs:
- CpG site filtering
- Zipping/compression of processed files
- Preparing data for downstream DSS analysis

---

### 03_dss_input
Preparation of DSS input objects:
- Merging CpG sites across samples
- Creating BSseq / HDF5-backed methylation matrices
- Generating DSS-compatible input objects for statistical testing

---

### 04_differential_methylation
Downstream statistical analysis and visualization:
- Pairwise differential methylation analysis using DSS
- Identification of DMCs and DMRs
- Genomic annotation of CpGs/DMRs
- Generation of plots:
  - Volcano plots
  - Feature distribution pie charts
  - Summary tables for interpretation

---

### 05_functional_annotation
-Annotating DMCs/DMRs to genomic features using UniProtKB/Swiss-Prot
-Assigning CpGs/DMRs to genes
-Gene-level summaries
-Ortholog mapping

---

### 06_enrichment_analysis
    - GO enrichment
    - KEGG enrichment
    - Reactome (if applicable)
    - Visualization (dotplots, cnetplots, enrichment maps)
    - Biological interpretation

---

## Snakemake Pipelines

Two Snakemake workflows are used:

### Main workflow
Located in:

snakemake_main/


Handles:
- QC
- trimming
- alignment
- methylation calling
- general EM-seq processing pipeline

---

### Control workflows
Located in:

snakemake_controls/


Separate pipelines for:
- Lambda control genome
- pUC19 control genome

These are used for quality control and conversion efficiency estimation.

---

## Slurm Execution

All HPC submission scripts are located in:


slurm_jobs/


To run a job:

```bash
sbatch script_name.sh

Snakemake workflows are also launched via Slurm submission scripts in this directory.

General Workflow Order
Prepare reference and control genomes (01_data_prep)
Run Snakemake main pipeline (snakemake_main)
Run control pipelines (snakemake_controls)
Post-process CpG files (02_post_processing)
Build DSS input objects (03_dss_input)
Run pairwise DSS analysis + visualization (04_differential_methylation)
Gene annotation & enrichment analysis (05_functional_annotation)
Enrichment analysis with GO/KEGG (06_enrichment_analysis)
Software Requirements
Snakemake (HPC execution mode)
Bismark
DSS (R package)
bsseq
HDF5Array
GenomicRanges
ChIPseeker
R (>= 4.2 recommended)
Slurm workload manager
Notes
All scripts are designed for HPC execution (Saga or similar systems).
Paths should be adapted using user-specific config files.
Sample names are placeholders in this repository and must be replaced with actual dataset identifiers.
DSS analyses assume a minimum coverage filter (default: ≥5 reads per CpG site).
