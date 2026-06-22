# Create DSS input

This script converts filtered Bismark CpG reports into an HDF5-backed BSseq object for downstream differential methylation analysis with DSS.

## Before running

Edit the following variables in `make_DSS_table.R`:

- working directory
- output directory
- input file list
- sample metadata

## Run

Load R:

module load R/<version>

Then execute:

Rscript make_DSS_table.R > make_DSS_table.log 2>&1

The output is an HDF5-backed BSseq object that can be used for DSS analysis.
