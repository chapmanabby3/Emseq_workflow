#!/bin/bash
#SBATCH --job-name=rfam_unresolved_dmrs
#SBATCH --account=<ACCT_NO>
#SBATCH --time=2-00:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
#SBATCH --output=rfam_unresolved_dmrs_%j.out
#SBATCH --error=rfam_unresolved_dmrs_%j.err

###############################################################################
# EDIT THESE SETTINGS
###############################################################################
# Create a new working directory inside atl_blastp_pfg for example rfam_unresolved #
WORKDIR="<INCLUDE_PATH_TO/NAME_YOUR_WORKING_DIRECTORY>"

GENOME_FASTA="<PATH/TO/polar_cod.fa>"

RFAM_DIR="<INCLUDE_PATH_TO/NAME_YOUR_RFAM_DIRECTORY>"
RFAM_CM="${RFAM_DIR}/Rfam.cm"
RFAM_CLANIN="${RFAM_DIR}/Rfam.clanin"
RFAM_BASE_URL="https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT"

# Input files will be same files from "annotated_file" from 06_merge_eggNOG.R #
INPUT_FILES=(
  "<PATH/TO/INPUT_FILES_FROM_EGGNOG_GROUP1_VS_GROUP2.csv>"
  "<PATH/TO/INPUT_FILES_FROM_EGGNOG_GROUP1_VS_GROUP3.csv>"
  "<PATH/TO/INPUT_FILES_FROM_EGGNOG_GROUP2_VS_GROUP3.csv>"
)

FLANK=2000
THREADS=16

###############################################################################

set -euo pipefail

mkdir -p "${WORKDIR}" "${RFAM_DIR}"
cd "${WORKDIR}"

module purge
module load R/4.5.2-gfbf-2025b
module load BEDTools/2.31.1-GCC-14.3.0
module load Infernal/1.1.5-gompi-2025b
module load SAMtools/1.22.1-GCC-14.3.0

###############################################################################
# Prepare a local genome copy with its FASTA index
#
# The project genome directory is at its file-count quota, so bedtools cannot
# create polar_cod.fa.fai there. Symlink the genome into the work directory
# and build the index locally instead.
###############################################################################

LOCAL_GENOME="${WORKDIR}/polar_cod.fa"

if [ ! -e "${LOCAL_GENOME}" ]; then
    ln -s "${GENOME_FASTA}" "${LOCAL_GENOME}"
fi

if [ ! -s "${LOCAL_GENOME}.fai" ]; then
    echo "Building FASTA index for ${LOCAL_GENOME} ..."
    samtools faidx "${LOCAL_GENOME}"
fi

###############################################################################
# Download and prepare the Rfam database once
###############################################################################

if [ ! -s "${RFAM_CM}" ]; then
    echo "Downloading Rfam.cm ..."
    wget -q -O "${RFAM_CM}.gz" "${RFAM_BASE_URL}/Rfam.cm.gz"
    gunzip -f "${RFAM_CM}.gz"
fi

if [ ! -s "${RFAM_CLANIN}" ]; then
    echo "Downloading Rfam.clanin ..."
    wget -q -O "${RFAM_CLANIN}" "${RFAM_BASE_URL}/Rfam.clanin"
fi

if [ ! -s "${RFAM_CM}.i1m" ]; then
    echo "Running cmpress on Rfam.cm ..."
    cmpress "${RFAM_CM}"
fi

for INPUT in "${INPUT_FILES[@]}"; do
  
  COMPARISON=$(basename "$(dirname "${INPUT}")")
  OUTDIR="${WORKDIR}/${COMPARISON}"
  mkdir -p "${OUTDIR}"
  
  echo "Processing ${COMPARISON}"
  
  Rscript - <<EOF
  
library(data.table)
library(tidyverse)

dmr <- fread("${INPUT}", na.strings = c("", "NA", "NaN"))

required <- c(
  "seqnames", "start", "end", "dmr_id",
  "final_annotation_tier", "feature_class"
)

missing <- setdiff(required, names(dmr))

if (length(missing) > 0) {
  stop(
    "Missing required columns in ${INPUT}:\n",
    paste(missing, collapse = ", ")
  )
}

unresolved <- dmr %>%
  filter(final_annotation_tier == "Unresolved") %>%
  mutate(
    chrom = as.character(seqnames),
    start = as.integer(start),
    end = as.integer(end),
    bed_start = pmax(0, start - 1 - ${FLANK}),
    bed_end = end + ${FLANK},
    region_id = dmr_id
  ) %>%
  distinct(chrom, bed_start, bed_end, region_id)

fwrite(
  unresolved[, .(chrom, bed_start, bed_end, region_id)],
  "${OUTDIR}/unresolved_dmrs.bed",
  sep = "\t",
  col.names = FALSE
)

message("Unresolved DMRs written: ", nrow(unresolved))

EOF
  
  bedtools getfasta \
    -fi "${LOCAL_GENOME}" \
    -bed "${OUTDIR}/unresolved_dmrs.bed" \
    -fo "${OUTDIR}/unresolved_dmrs.fa" \
    -name
  
  cmscan \
    --cpu "${THREADS}" \
    --rfam \
    --nohmmonly \
    --cut_ga \
    --fmt 2 \
    --clanin "${RFAM_CLANIN}" \
    --oskip \
    --tblout "${OUTDIR}/unresolved_dmrs_rfam.tbl" \
    -o "${OUTDIR}/unresolved_dmrs_rfam.txt" \
    "${RFAM_CM}" \
    "${OUTDIR}/unresolved_dmrs.fa"
  
  grep -v "^#" "${OUTDIR}/unresolved_dmrs_rfam.tbl" \
    | grep -v "^$" \
    > "${OUTDIR}/unresolved_dmrs_rfam.clean.tbl"
  
  echo "Finished ${COMPARISON}"
  
done

echo "All Rfam searches complete."
