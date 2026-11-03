###############################################################################
# Merge Rfam/Infernal ncRNA results into DMR tables
#
# Input: *_with_atlantic_cod.csv files containing methylation_direction and feature_class
# Output: *_with_rfam.csv
#
# Adds:
#   rfam_family
#   rfam_accession
#   rfam_evalue
#   rfam_score
#   rfam_query_start
#   rfam_query_end
#   rfam_strand
#   rfam_description
#   has_rfam_hit
#
# Updates:
#   final_annotation_tier
#   final_gene_label
#   feature_class
###############################################################################

library(data.table)
library(tidyverse)

setwd("<PATH/TO/WORKING_DIRECTORY")

rfam_dir <- "rfam_unresolved"

input_files <- c(
  "high_vs_control/<GROUP1_VS_GROUP2_MOST_RECENT_ANNOTATED_DMR_TABLES.csv>",
  "low_vs_control/<GROUP1_VS_GROUP3_MOST_RECENT_ANNOTATED_DMR_TABLES.csv>",
  "high_vs_low/<GROUP2_VS_GROUP3_MOST_RECENT_ANNOTATED_DMR_TABLES.csv>"
)

overwrite_inputs <- FALSE

parse_rfam_tbl <- function(tbl_file) {
  
  if (!file.exists(tbl_file)) {
    warning("Rfam table not found: ", tbl_file)
    return(NULL)
  }
  
  raw <- readLines(tbl_file)
  raw <- raw[!str_detect(raw, "^#")]
  raw <- raw[raw != ""]
  
  if (length(raw) == 0) {
    message("No Rfam hits found in: ", tbl_file)
    return(NULL)
  }
  
  # cmscan --fmt 2 tblout is whitespace-delimited with a variable number of
  # trailing description words. Split each line manually:
  #   field 1  = row index
  #   field 2  = Rfam family name
  #   field 3  = Rfam accession
  #   field 4  = target name (bedtools -name output:
  #              "dmr_id::chr:flankedStart-flankedEnd")
  #   fields 5-26 = fixed columns; description = fields 29+
  parsed <- str_split(raw, "\\s+")

  hits_clean <- tibble(
    rfam_family = map_chr(parsed, 2),
    rfam_accession = map_chr(parsed, 3),
    dmr_id = map_chr(parsed, ~ str_split(.x[4], "::", simplify = TRUE)[1]),
    rfam_strand = map_chr(parsed, 12),
    rfam_score = suppressWarnings(as.numeric(map_chr(parsed, 17))),
    rfam_evalue = suppressWarnings(as.numeric(map_chr(parsed, 18))),
    rfam_query_start = suppressWarnings(as.numeric(map_chr(parsed, 27))),
    rfam_query_end = suppressWarnings(as.numeric(map_chr(parsed, 28))),
    rfam_description = map_chr(parsed, function(x) {
      if (length(x) >= 29) paste(x[29:length(x)], collapse = " ") else NA_character_
    })
  )

  hits_clean %>%
    arrange(dmr_id, rfam_evalue, desc(rfam_score)) %>%
    distinct(dmr_id, .keep_all = TRUE)
}

for (input_file in input_files) {
  
  message("\n=============================================================")
  message("Processing: ", input_file)
  message("=============================================================")
  
  if (!file.exists(input_file)) {
    warning("File not found, skipping: ", input_file)
    next
  }
  
  comparison <- dirname(input_file)
  
  tbl_file <- file.path(
    rfam_dir,
    comparison,
    "unresolved_dmrs_rfam.clean.tbl"
  )
  
  rfam_hits <- parse_rfam_tbl(tbl_file)
  
  dmr <- fread(input_file, na.strings = c("", "NA", "NaN"))
  
  if (is.null(rfam_hits)) {
    dmr <- dmr %>%
      mutate(
        rfam_family = NA_character_,
        rfam_accession = NA_character_,
        rfam_evalue = NA_real_,
        rfam_score = NA_real_,
        rfam_query_start = NA_real_,
        rfam_query_end = NA_real_,
        rfam_strand = NA_character_,
        rfam_description = NA_character_,
        has_rfam_hit = FALSE
      )
  } else {
    dmr <- dmr %>%
      left_join(
        rfam_hits,
        by = "dmr_id"
      ) %>%
      mutate(
        has_rfam_hit = !is.na(rfam_family) & rfam_family != ""
      )
  }
  
  dmr <- dmr %>%
    mutate(
      final_annotation_tier = case_when(
        has_swissprot ~ "Swiss-Prot_supported",
        
        has_usable_atl & has_named_atl ~
          "Atlantic_cod_ortholog_supported",
        
        has_eggnog_name ~ "eggNOG_ortholog_supported",
        
        has_usable_atl ~ "Uncharacterized_Atlantic_cod_locus",
        
        has_eggnog_function ~ "eggNOG_function_supported",
        
        !is.na(eggnog_PFAMs) & eggnog_PFAMs != "-" ~
          "eggNOG_PFAM_domain_supported",
        
        has_rfam_hit ~ "Rfam_ncRNA_supported",
        
        TRUE ~ "Unresolved"
      ),
      
      final_gene_label = case_when(
        has_swissprot & !is.na(gene_symbol) &
          gene_symbol != "" & gene_symbol != geneId ~
          gene_symbol,
        
        has_usable_atl & has_named_atl ~
          blast_atl_gene_label,
        
        has_eggnog_name ~ eggnog_preferred_name,
        
        has_rfam_hit & !is.na(rfam_family) ~
          paste0(rfam_family, "_ncRNA"),
        
        TRUE ~ NA_character_
      ),
      
      feature_class = case_when(
        has_rfam_hit & final_annotation_tier == "Rfam_ncRNA_supported" ~
          "ncRNA",
        
        TRUE ~ feature_class
      )
    )
  
  if (overwrite_inputs) {
    output_file <- input_file
  } else {
    output_file <- str_replace(
      input_file,
      "\\.csv$",
      "_with_rfam.csv"
    )
  }
  
  fwrite(dmr, output_file)
  
  message("Wrote: ", output_file)
  
  message("\nFinal annotation tier summary:")
  print(dmr %>% count(final_annotation_tier))
  
  message("\nFeature class summary:")
  print(dmr %>% count(feature_class))
  
  message("\nRfam families detected:")
  print(
    dmr %>%
      filter(has_rfam_hit) %>%
      count(rfam_family, sort = TRUE)
  )
}

message("\nAll Rfam results merged.")
