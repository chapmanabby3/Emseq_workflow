#!/bin/bash
#SBATCH --job-name=eggnog_polar_cod_pfg
#SBATCH --account=<ACCNT NO>
#SBATCH --time=24:00:00
#SBATCH --mem=100G
#SBATCH --cpus-per-task=32
#SBATCH --output=eggnog_polar_cod_pfg_%j.out
#SBATCH --error=eggnog_polar_cod_pfg_%j.err

###############################################################################
# EDIT THESE SETTINGS
###############################################################################
# script expects wd to be set at <PATH TO PREVIOUS WD>/eggnog_mapper #
WORKDIR=<PATH/TO/WD>

# Location where eggNOG databases will be stored.
# Make sure this is in a directory with enough space.
# Create a data directory such as working directory/eggnog_data #
EGGNOG_DATA_DIR=<PATH/TO/EGGNOG_DATA_DIRECTORY>

CPU=32

###############################################################################

set -euo pipefail

module purge
module load Miniconda3/25.7.0-2

# Activate or create an eggNOG environment.
# Use only conda-forge/bioconda and override the default Anaconda channels,
# whose Terms of Service cannot be accepted non-interactively in a batch job.
if ! conda env list | grep -q "^eggnog_env "; then
    conda create -y -n eggnog_env \
        --override-channels \
        -c conda-forge -c bioconda \
        eggnog-mapper diamond
fi

source activate eggnog_env

mkdir -p "$WORKDIR" "$EGGNOG_DATA_DIR"
cd "$WORKDIR"

###############################################################################
# Download eggNOG databases once
#
# NOTE: download_eggnog_data.py uses the hardcoded host eggnogdb.embl.de,
# which no longer resolves in DNS (NXDOMAIN). We therefore download the three
# required files directly from the working mirror eggnog5.embl.de instead.
###############################################################################

EGGNOG_BASE_URL="http://eggnog5.embl.de/download/emapperdb-5.0.2"

download_if_absent() {
    local file="$1"
    local dest="$EGGNOG_DATA_DIR/$file"
    if [ ! -s "$dest" ]; then
        echo "Downloading $file ..."
        wget -q -O "$dest" "$EGGNOG_BASE_URL/$file"
    fi
}

if [ ! -f "$EGGNOG_DATA_DIR/eggnog.db" ]; then
    download_if_absent eggnog.db.gz
    download_if_absent eggnog.taxa.tar.gz
    download_if_absent eggnog_proteins.dmnd.gz

    echo "Decompressing eggNOG databases ..."
    gunzip -f "$EGGNOG_DATA_DIR/eggnog.db.gz"
    gunzip -f "$EGGNOG_DATA_DIR/eggnog_proteins.dmnd.gz"
    tar -xzf "$EGGNOG_DATA_DIR/eggnog.taxa.tar.gz" -C "$EGGNOG_DATA_DIR"
fi

###############################################################################
# Run eggNOG-mapper
###############################################################################

if [ ! -s unresolved_polar_cod_proteins_for_eggnog.faa ]; then
    echo "ERROR: query FASTA unresolved_polar_cod_proteins_for_eggnog.faa not found in $WORKDIR" >&2
    echo "Run the prep script (prep_eggNOG_quer_pfg.R) first to generate it." >&2
    exit 1
fi

emapper.py \
    -i unresolved_polar_cod_proteins_for_eggnog.faa \
    -o polar_cod_eggnog \
    --output_dir "$WORKDIR" \
    --data_dir "$EGGNOG_DATA_DIR" \
    -m diamond \
    --itype proteins \
    --cpu "$CPU" \
    --override

echo "eggNOG-mapper finished."
