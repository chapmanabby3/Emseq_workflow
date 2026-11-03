###############################################################################
# Select polar cod proteins for eggNOG-mapper
#
# Selects proteins that:
#   1. Have no Swiss-Prot hit, AND
#   2. Have no Atlantic cod hit, OR
#      have a weak Atlantic cod hit, OR
#      match an uncharacterized/LOC Atlantic cod protein
#
# Output:
#   - One non-redundant FASTA of proteins for eggNOG-mapper
#   - A CSV key linking protein_id to DMR files
###############################################################################

library(data.table)
library(tidyverse)

###############################################################################
# EDIT THESE SETTINGS
###############################################################################
### Set working directory as same wd from atl cod blast ###

setwd("<PATH/TO/WORKING_DIRECTORY>")

input_files <- c(
  "high_vs_control/DMRs_with_gene_symbols_high_vs_control_pfg_with_atlantic_cod.csv",
  "low_vs_control/DMRs_with_gene_symbols_low_vs_control_pfg_with_atlantic_cod.csv",
  "high_vs_low/DMRs_with_gene_symbols_high_vs_low_pfg_with_atlantic_cod.csv"
)

polar_cod_proteins_faa <- "<PATH/TO/GCA_964260565.1_borsaiasm_protein.faa>"

eggnog_dir <- "eggnog_mapper"
dir.create(eggnog_dir, showWarnings = FALSE)

query_fasta <- file.path(eggnog_dir, "unresolved_polar_cod_proteins_for_eggnog.faa")
key_csv <- file.path(eggnog_dir, "unresolved_protein_key.csv")

###############################################################################

read_dmr_table <- function(file) {
  
  tbl <- fread(
    file,
    na.strings = c("", "NA", "NaN")
  )
  
  required_columns <- c(
    "dmr_id",
    "protein_id",
    "sseqid",
    "blast_atl_sseqid",
    "blast_atl_annotation_quality",
    "blast_atl_protein_name"
  )
  
  missing_columns <- setdiff(required_columns, colnames(tbl))
  
  if (length(missing_columns) > 0) {
    stop(
      "Missing required columns in ", file, ":\n",
      paste(missing_columns, collapse = ", ")
    )
  }
  
  tbl %>%
    mutate(
      source_file = basename(file)
    )
}

dmr_all <- rbindlist(
  lapply(input_files, read_dmr_table),
  use.names = TRUE,
  fill = TRUE
)

# Clean Atlantic cod protein descriptions
dmr_all <- dmr_all %>%
  mutate(
    protein_id = as.character(protein_id),
    sseqid = as.character(sseqid),
    blast_atl_sseqid = as.character(blast_atl_sseqid),
    blast_atl_annotation_quality = as.character(blast_atl_annotation_quality),
    blast_atl_protein_name = as.character(blast_atl_protein_name),
    
    no_swissprot = is.na(sseqid) | sseqid == "" | sseqid == "NA",
    
    no_atl_hit = is.na(blast_atl_sseqid) |
      blast_atl_sseqid == "" |
      blast_atl_sseqid == "NA" |
      blast_atl_annotation_quality == "no_atlantic_cod_hit",
    
    weak_atl_hit = blast_atl_annotation_quality == "weak_atlantic_cod_hit",
    
    uncharacterized_atl_hit = str_detect(
      tolower(coalesce(blast_atl_protein_name, "")),
      "uncharacterized|loc\\d+|hypothetical"
    )
  )

# Proteins worth sending to eggNOG
selected <- dmr_all %>%
  filter(
    no_swissprot,
    no_atl_hit | weak_atl_hit | uncharacterized_atl_hit
  ) %>%
  distinct(protein_id, source_file) %>%
  arrange(protein_id, source_file)

n_proteins <- n_distinct(selected$protein_id)

message("Unique polar cod proteins selected for eggNOG-mapper: ", n_proteins)

if (n_proteins == 0) {
  stop("No unresolved proteins were selected. Check the input files and filters.")
}

fwrite(selected, key_csv)

###############################################################################
# Extract selected protein sequences
###############################################################################

message("Reading polar cod protein FASTA...")

fasta_lines <- readLines(polar_cod_proteins_faa)

header_positions <- which(str_detect(fasta_lines, "^>"))

if (length(header_positions) == 0) {
  stop("No FASTA headers found in: ", polar_cod_proteins_faa)
}

sequence_start <- header_positions + 1
sequence_end <- c(
  header_positions[-1] - 1,
  length(fasta_lines)
)

fasta_ids <- str_remove(fasta_lines[header_positions], "^>") %>%
  str_extract("^\\S+")

selected_ids <- unique(selected$protein_id)

selected_positions <- which(fasta_ids %in% selected_ids)

if (length(selected_positions) == 0) {
  stop(
    "None of the selected protein IDs were found in the polar cod FASTA.\n",
    "Check that protein_id values match the FASTA headers."
  )
}

output_lines <- unlist(
  lapply(
    selected_positions,
    function(i) {
      c(
        fasta_lines[header_positions[i]],
        fasta_lines[sequence_start[i]:sequence_end[i]]
      )
    }
  )
)

writeLines(output_lines, query_fasta)

message("Wrote eggNOG query FASTA: ", query_fasta)
message("Proteins written: ", length(selected_positions))
message("Protein-to-DMR key written: ", key_csv)
