#!/bin/bash

###############################################################################
# Functional Annotation
# Step 1: Download annotation resources and perform BLASTP
#
# Description:
#   Downloads the predicted protein FASTA, genome feature table,
#   and UniProtKB/Swiss-Prot database. Creates a local BLAST database
#   and identifies the best Swiss-Prot match for each predicted protein.
#
# Requirements:
#   - BLAST+
#   - wget
#
# Inputs:
#   - Protein FASTA (.faa)
#   - Genome feature table
#
# Outputs:
#   - blast_results.txt
#   - Swiss-Prot BLAST database (sprot_db.*)
###############################################################################

# Create working directory
mkdir -p annotation_blast
cd annotation_blast

###############################################################################
# Download annotation files
###############################################################################

# Predicted protein FASTA
wget <PROTEIN_FASTA_URL>

# Genome feature table
wget <FEATURE_TABLE_URL>

# UniProtKB/Swiss-Prot database
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz

###############################################################################
# Uncompress downloaded files
###############################################################################

gunzip *.gz

###############################################################################
# Load BLAST
###############################################################################

module load BLAST+/2.17.0-gompi-2025b

###############################################################################
# Create local BLAST database
###############################################################################

makeblastdb \
    -in uniprot_sprot.fasta \
    -dbtype prot \
    -out sprot_db

###############################################################################
# Run BLASTP
###############################################################################

blastp \
    -query <PROTEIN_FASTA_FILE>.faa \
    -db sprot_db \
    -out blast_results.txt \
    -outfmt "6 qseqid sseqid pident evalue stitle" \
    -evalue 1e-5 \
    -num_threads 8 \
    -max_target_seqs 1

###############################################################################