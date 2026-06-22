# Post-processing

This directory contains scripts used after Bismark methylation extraction.

Workflow:

1. Run `CpG_filter.sh` to retain only CpG methylation calls from Bismark CX reports.

2. Run `zipper.sh` to compress the resulting `.CpG.txt` files prior to downstream analyses.
