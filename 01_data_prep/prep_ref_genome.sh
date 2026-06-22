#!/bin/bash

# Load the current Bismark module available on Saga
# module avail Bismark
# module load Bismark/<version>

# Download the current polar cod assembly

wget ...

gunzip ...

mv ...

# Prepare the genome for Bismark

bismark_genome_preparation --verbose .
bismark_genome_preparation --bowtie2 .
