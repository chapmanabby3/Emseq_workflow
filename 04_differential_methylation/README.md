Overview

This directory contains R scripts for downstream analysis of EM-seq / methylation data after the DSS pipeline has been executed.

The analysis includes:

Pairwise differential methylation analysis outputs (DSS results)
Annotation of DMCs and DMRs
Generation of volcano plots
Genomic feature enrichment visualizations (pie charts)
Scripts in this directory
1. DSS_pairwise.R

This script performs downstream processing of the DSS results generated from the EM-seq pipeline.

It:

Loads the processed BSseq object
Applies sample metadata
Filters CpG sites based on coverage thresholds
Performs pairwise comparisons between experimental groups
Generates:
DML (differential methylation loci)
DMC (differentially methylated CpGs)
DMR (differentially methylated regions)
Saves results in RDS and CSV formats for downstream analysis
2. volc_pie_annotation.R

This script performs downstream biological interpretation and visualization of DSS results.

It:

Loads DML/DMC/DMR results from pairwise comparisons
Annotates genomic features using a reference GTF file
Generates:
Volcano plots of methylation differences
Annotated gene-level summaries
Pie charts of genomic feature distribution (Promoter, Exon, Intron, etc.)
Produces figures for all defined pairwise comparisons
Required inputs

Before running these scripts, ensure that the following outputs exist:

DSS pairwise results:

DSS_pairwise_complete/
    high_vs_control/
    low_vs_control/
    high_vs_low/

HDF5-backed BSseq object:

emseq_hdf5/

Genome annotation file (GTF):

*.gtf.gz
Execution workflow

These scripts are intended to be run in sequence:

Step 1 – Run DSS analysis

Submit via SLURM:

sbatch slurm_jobs/dss_pairwise.sh

This generates:

DML test results
DMC results
DMR results
Step 2 – Run annotation and plotting

Once DSS results are available, run:

sbatch slurm_jobs/volc_pie_annotation.sh

This produces:

Volcano plots
Annotated DMC/DMR tables
Pie charts of genomic feature distribution
