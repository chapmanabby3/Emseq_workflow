###############################################################################
# Merge DMR methylation direction, gene/BLAST annotation, and eggNOG results
#
# For each comparison, this script:
#   1. Reads DMRs_filtered_*.txt, which contains methylation values
#   2. Reads DMRs_with_gene_symbols_..._with_atlantic_cod.csv
#   3. Joins methylation information by dmr_id, or chromosome/start/end
#   4. Adds methylation_direction and feature_class
#   5. Merges eggNOG annotations by protein_id
#   6. Writes/overwrites *_with_eggnog.csv
#
# Direction definition:
#   diff.Methy = meanMethy1 - meanMethy2
#   Positive = hypermethylated in group 1
#   Negative = hypomethylated in group 1
###############################################################################

library(data.table)
library(tidyverse)


###############################################################################
# EDIT THESE SETTINGS
###############################################################################
# Simplest to choose same working directory from 03_blastp_at.R #
setwd("<PATH/TO/WORKING_DIRECTORY>")

# Each element contains:
#   filtered_file: original DMR-calling output with methylation values
#   annotated_file: DMR table with polar cod genes and Atlantic cod BLASTP
#   group1: biological label for meanMethy1
#   group2: biological label for meanMethy2

comparisons <- list(
  
  list(
    name = "high_vs_control",
    filtered_file = "<PATH/TO/ORIGINAL DMR-CALLING OUTPUT WITH METHYLATION VALUES.txt>",
    annotated_file = "<PATH/TO/DRM_TABLE_WITH_PC_GNES_AND_ATLANTIC_COD_BLASTP.csv>",
    group1 = "high",
    group2 = "control"
  ),
  
  list(
    name = "low_vs_control",
    filtered_file = "<PATH/TO/ORIGINAL DMR-CALLING OUTPUT WITH METHYLATION VALUES.txt>",
    annotated_file = "<PATH/TO/DRM_TABLE_WITH_PC_GNES_AND_ATLANTIC_COD_BLASTP.csv>",
    group1 = "low",
    group2 = "control"
  ),
  
  list(
    name = "high_vs_low",
    filtered_file = "<PATH/TO/ORIGINAL DMR-CALLING OUTPUT WITH METHYLATION VALUES.txt>",
    annotated_file = "<PATH/TO/DRM_TABLE_WITH_PC_GNES_AND_ATLANTIC_COD_BLASTP.csv>",
    group1 = "high",
    group2 = "low"
  )
)

eggnog_annotations <- "eggnog_mapper/polar_cod_eggnog.emapper.annotations"

# This script overwrites the existing *_with_eggnog.csv files.
overwrite_inputs <- TRUE


###############################################################################


###############################################################################
# Read eggNOG output
###############################################################################

message("Reading eggNOG annotations: ", eggnog_annotations)

eggnog_raw <- readLines(eggnog_annotations)

header_line <- which(str_detect(eggnog_raw, "^#query"))
header_line <- header_line[1]

if (is.na(header_line)) {
  stop("Could not find eggNOG header line beginning with '#query'.")
}

eggnog <- read_tsv(
  eggnog_annotations,
  comment = "##",
  show_col_types = FALSE
)

eggnog <- eggnog %>%
  rename_with(
    ~ str_remove(.x, "^#"),
    everything()
  )

if ("max_annot_lvl" %in% colnames(eggnog) &&
    !"max_og_lvl" %in% colnames(eggnog)) {
  eggnog <- eggnog %>%
    rename(max_og_lvl = max_annot_lvl)
}

expected_cols <- c(
  "query",
  "seed_ortholog",
  "evalue",
  "score",
  "eggNOG_OGs",
  "max_og_lvl",
  "COG_category",
  "Description",
  "Preferred_name",
  "GOs",
  "EC",
  "KEGG_ko",
  "KEGG_Pathway",
  "KEGG_Module",
  "KEGG_Reaction",
  "PFAMs"
)

missing_cols <- setdiff(expected_cols, colnames(eggnog))

if (length(missing_cols) > 0) {
  message(
    "Warning: these expected eggNOG columns were not found:\n",
    paste(missing_cols, collapse = ", ")
  )
}

eggnog_clean <- eggnog %>%
  transmute(
    protein_id = as.character(query),
    
    eggnog_seed_ortholog = seed_ortholog,
    eggnog_evalue = as.numeric(evalue),
    eggnog_score = as.numeric(score),
    eggnog_eggNOG_OGs = eggNOG_OGs,
    eggnog_max_og_lvl = max_og_lvl,
    eggnog_COG_category = COG_category,
    eggnog_description = Description,
    eggnog_preferred_name = Preferred_name,
    eggnog_GOs = GOs,
    eggnog_EC = EC,
    eggnog_KEGG_ko = KEGG_ko,
    eggnog_KEGG_Pathway = KEGG_Pathway,
    eggnog_KEGG_Module = KEGG_Module,
    eggnog_KEGG_Reaction = KEGG_Reaction,
    eggnog_PFAMs = PFAMs
  ) %>%
  distinct(protein_id, .keep_all = TRUE)

message("eggNOG proteins annotated: ", nrow(eggnog_clean))


###############################################################################
# Function: merge filtered DMR methylation data with annotated DMR table
###############################################################################

merge_methylation_direction <- function(filtered_file,
                                        annotated_file,
                                        group1,
                                        group2) {
  
  message("\nReading filtered DMR file: ", filtered_file)
  filtered <- fread(filtered_file, na.strings = c("", "NA", "NaN"))
  
  message("Reading annotated DMR file: ", annotated_file)
  annotated <- fread(annotated_file, na.strings = c("", "NA", "NaN"))
  
  required_filtered <- c(
    "chr",
    "start",
    "end",
    "meanMethy1",
    "meanMethy2",
    "diff.Methy"
  )
  
  required_annotated <- c(
    "seqnames",
    "start",
    "end"
  )
  
  missing_filtered <- setdiff(required_filtered, names(filtered))
  missing_annotated <- setdiff(required_annotated, names(annotated))
  
  if (length(missing_filtered) > 0) {
    stop(
      "Missing columns in filtered DMR file:\n",
      paste(missing_filtered, collapse = ", ")
    )
  }
  
  if (length(missing_annotated) > 0) {
    stop(
      "Missing columns in annotated DMR file:\n",
      paste(missing_annotated, collapse = ", ")
    )
  }
  
  # Standardize filtered DMR table
  filtered <- filtered %>%
    mutate(
      chr = as.character(chr),
      start = as.integer(start),
      end = as.integer(end),
      
      mean_methylation_group1 = as.numeric(meanMethy1),
      mean_methylation_group2 = as.numeric(meanMethy2),
      methylation_difference = as.numeric(diff.Methy),
      
      methylation_direction = case_when(
        is.na(methylation_difference) ~ NA_character_,
        methylation_difference > 0 ~
          paste0("hypermethylated_in_", group1),
        methylation_difference < 0 ~
          paste0("hypomethylated_in_", group1),
        methylation_difference == 0 ~ "no_difference"
      )
    )
  
  # Standardize annotated DMR table
  annotated <- annotated %>%
    mutate(
      seqnames = as.character(seqnames),
      start = as.integer(start),
      end = as.integer(end),
      protein_id = as.character(protein_id)
    )
  
  # Prefer joining by dmr_id if both files contain it.
  if ("dmr_id" %in% names(filtered) && "dmr_id" %in% names(annotated)) {
    
    message("Joining filtered and annotated DMR tables by dmr_id.")
    
    merged <- annotated %>%
      left_join(
        filtered %>%
          select(
            dmr_id,
            mean_methylation_group1,
            mean_methylation_group2,
            methylation_difference,
            methylation_direction
          ) %>%
          distinct(dmr_id, .keep_all = TRUE),
        by = "dmr_id"
      )
    
  } else {
    
    message("dmr_id not available in both files; joining by genomic coordinates.")
    
    merged <- annotated %>%
      left_join(
        filtered %>%
          select(
            chr,
            start,
            end,
            mean_methylation_group1,
            mean_methylation_group2,
            methylation_difference,
            methylation_direction
          ) %>%
          distinct(chr, start, end, .keep_all = TRUE),
        by = c(
          "seqnames" = "chr",
          "start" = "start",
          "end" = "end"
        )
      )
  }
  
  unmatched <- sum(is.na(merged$methylation_direction))
  
  message("DMR rows without matched methylation direction: ", unmatched)
  
  if (unmatched > 0) {
    warning(
      unmatched,
      " annotated DMR rows did not match a filtered DMR. ",
      "Check chromosome/start/end formatting or dmr_id consistency."
    )
  }
  
  merged %>%
    mutate(
      feature_class = "protein_coding"
    )
}


###############################################################################
# Process each comparison
###############################################################################

for (comp in comparisons) {
  
  message("\n=============================================================")
  message("Processing comparison: ", comp$name)
  message("Group 1: ", comp$group1)
  message("Group 2: ", comp$group2)
  message("=============================================================")
  
  if (!file.exists(comp$filtered_file)) {
    warning("Filtered DMR file not found, skipping: ", comp$filtered_file)
    next
  }
  
  if (!file.exists(comp$annotated_file)) {
    warning("Annotated DMR file not found, skipping: ", comp$annotated_file)
    next
  }
  
  dmr <- merge_methylation_direction(
    filtered_file = comp$filtered_file,
    annotated_file = comp$annotated_file,
    group1 = comp$group1,
    group2 = comp$group2
  )
  
  dmr_eggnog <- dmr %>%
    left_join(
      eggnog_clean,
      by = "protein_id"
    )
  
  dmr_eggnog <- dmr_eggnog %>%
    mutate(
      has_swissprot = !is.na(sseqid) &
        sseqid != "" &
        sseqid != "NA",
      
      has_usable_atl = !is.na(blast_atl_annotation_quality) &
        blast_atl_annotation_quality == "usable_atlantic_cod_hit",
      
      has_named_atl = !is.na(blast_atl_gene_label) &
        blast_atl_gene_label != "" &
        !str_detect(blast_atl_gene_label, "uncharacterized|LOC\\d+"),
      
      has_eggnog_name = !is.na(eggnog_preferred_name) &
        eggnog_preferred_name != "" &
        eggnog_preferred_name != "-",
      
      has_eggnog_function = !is.na(eggnog_description) &
        eggnog_description != "" &
        eggnog_description != "-" &
        !str_detect(
          tolower(eggnog_description),
          "uncharacterized|hypothetical"
        ),
      
      final_annotation_tier = case_when(
        has_swissprot ~ "Swiss-Prot_supported",
        
        has_usable_atl & has_named_atl ~
          "Atlantic_cod_ortholog_supported",
        
        has_eggnog_name ~ "eggNOG_ortholog_supported",
        
        has_usable_atl ~ "Uncharacterized_Atlantic_cod_locus",
        
        has_eggnog_function ~ "eggNOG_function_supported",
        
        !is.na(eggnog_PFAMs) & eggnog_PFAMs != "-" ~
          "eggNOG_PFAM_domain_supported",
        
        TRUE ~ "Unresolved"
      ),
      
      final_gene_label = case_when(
        has_swissprot & !is.na(gene_symbol) &
          gene_symbol != "" & gene_symbol != geneId ~
          gene_symbol,
        
        has_usable_atl & has_named_atl ~
          blast_atl_gene_label,
        
        has_eggnog_name ~ eggnog_preferred_name,
        
        TRUE ~ NA_character_
      )
    )
  
  if (overwrite_inputs) {
    output_file <- comp$annotated_file
  } else {
    output_file <- str_replace(
      comp$annotated_file,
      "\\.csv$",
      "_with_eggnog.csv"
    )
  }
  
  fwrite(dmr_eggnog, output_file)
  
  message("Wrote: ", output_file)
  
  message("\nMethylation direction summary:")
  print(dmr_eggnog %>% count(methylation_direction))
  
  message("\nFinal annotation tier summary:")
  print(dmr_eggnog %>% count(final_annotation_tier))
  
  message("\nGenes with a final label: ",
          sum(!is.na(dmr_eggnog$final_gene_label)))
}

message("\nAll methylation directions, gene annotations, and eggNOG results merged.")
