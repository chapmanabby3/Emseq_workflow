#!/bin/bash
#
# ------------------------------------------------------------
# prep_control_genomes.sh
#
# Prepares Bismark indices for the control genomes used to
# assess bisulfite conversion efficiency.
#
# Required files:
#   - pUC19.fa
#   - lambda.fa
#
# Place both FASTA files in the current directory before
# running this script.
#
# Example:
#   bash prep_control_genomes.sh
# ------------------------------------------------------------

echo "Preparing Bismark indices for control genomes..."

bismark_genome_preparation --verbose .
bismark_genome_preparation --bowtie2 .

echo "Finished preparing control genome indices."
