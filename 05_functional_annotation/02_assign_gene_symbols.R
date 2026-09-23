###############################################################################
# Files required for downstream analysis
###############################################################################

# Transfer the following files to the directory containing the R analysis:
#
#   blast_results.txt
#   <FEATURE_TABLE>.txt
#   <ANNOTATED_DMC_TABLE>.txt
#   <ANNOTATED_DMR_TABLE>.txt
#
# I create 3 directories inside working directory: group1_vs_group2, group3_vs_group2, group1_vs_group3
# I place the respective DMC & DMR .txt files in each sub-directory
###############################################################################

###############################################################################
# Functional Annotation
# Step 2: Assign gene symbols to annotated DMCs
#
# Description:
#   Maps predicted proteins to UniProtKB/Swiss-Prot annotations and
#   assigns standardized gene symbols to DMC-associated loci.
#
# Inputs:
#   - blast_results.txt
#   - feature_table.txt
#   - annotated_DMCs.txt
#   - annotated_DMRs.txt
#
# Output:
#   - DMCs_with_gene_symbols.csv
#   - DMRs_with_gene_symbols.csv
###############################################################################

library(data.table)
library(tidyverse)

# ============================================================
# READ THE RAW FILE AND RE-PARSEAR
# ============================================================

# Read line by line and divide by the real number.#

#input blast_results.txt as whatever you named your blast results .txt file#

setwd("<WORKING_DIRECTORY>")

blast_file   <- "blast_results.txt"
feature_file <- "<FEATURE_TABLE>.txt"

#####################################################

# One row = one pairwise comparison to annotate
comparisons <- tibble(
  comparison = c(
    "high_vs_control",
    "low_vs_control",
    "high_vs_low"
  ),
dmc_file = c(
  "group1_vs_group3/DSS_sign_DMCs_annotated_sex_group1_vs_group3.txt",
  "group2_vs_group3/DSS_sign_DMCs_annotated_sex_group2_vs_group3.txt",
  "group1_vs_group2/DSS_sign_DMCs_annotated_sex_group1_vs_group2.txt"
),
dmc_output_file = c(
  "group1_vs_group3/DMCs_with_gene_symbols_group1_vs_group3_pfg.csv",
  "group2_vs_group3/DMCs_with_gene_symbols_group2_vs_group3_pfg.csv",
  "group1_vs_group2/DMCs_with_gene_symbols_group1_vs_group2_pfg.csv"
),
dmr_file = c(
  "group1_vs_group3/DSS_DMRs_annotated_sex_group1_vs_group3.txt",
  "group2_vs_group3/DSS_DMRs_annotated_sex_group2_vs_group3.txt",
  "group1_vs_group2/DSS_DMRs_annotated_sex_group1_vs_group2.txt"
),
dmr_output_file = c(
  "group1_vs_group3/DMRs_with_gene_symbols_group1_vs_group3_pfg.csv",
  "group2_vs_group3/DMRs_with_gene_symbols_group2_vs_group3_pfg.csv",
  "group1_vs_group2/DMRs_with_gene_symbols_group1_vs_group2_pfg.csv"
 )
)

blast_raw <- readLines("blast_results.txt")
# Separate by actual tabulations #
blast_parsed <- str_split(blast_raw, "\t", simplify = TRUE)
colnames(blast_parsed) <- c("protein_id", "sseqid", "pident", "evalue", "stitle")

blast_res <- as.data.table(blast_parsed) %>%
  mutate(
    protein_id = str_trim(protein_id),
    pident     = as.numeric(pident),
    evalue     = as.numeric(evalue),
    gene_symbol = str_extract(stitle, "(?<=GN=)\\S+"),
    protein_name = str_trim(str_extract(stitle, "^[^O]+(?=OS=)"))
  ) %>%
  # Mejor hit por proteína
  group_by(protein_id) %>%
  slice_min(evalue, n = 1, with_ties = FALSE) %>%
  ungroup()

# Verify
head(blast_res %>% select(protein_id, gene_symbol, pident, evalue))
#output should look as below#
#       protein_id gene_symbol  pident    evalue
# CAL8233918.1       Vgll4  52.874  2.65e-78
# CAL8233919.1        Syn2  67.562       0.0

# ============================================================
# JOIN with feature table and methylation
# ============================================================
feat <- fread("<FEATURE_TABLE>.txt", 
              skip = 1, # Skips the line with the '#'
              header = FALSE)
colnames(feat) <- c(
  "feature", "class", "assembly", "assembly_unit", "seq_type", 
  "chromosome", "genomic_accession", "start", "end", "strand", 
  "product_accession", "non_redundant_refseq", "related_accession", 
  "name", "symbol", "GeneID", "locus_tag", 
  "feature_interval_length", "product_length", "attributes"
)

locus_to_protein <- feat %>%
  filter(feature == "CDS", product_accession != "") %>%
  select(locus_tag, protein_id = product_accession) %>%
  distinct()

locus_to_symbol <- locus_to_protein %>%
  left_join(blast_res, by = "protein_id")

######################
# Both DMC and DMR tables are plain tab-delimited files with clean headers,
# so the same join/annotate logic works for either one.

annotate_and_write <- function(in_file, out_file, label) {
  message("Processing ", label, ": ", in_file)

  tbl <- read_tsv(
    in_file,
    trim_ws = TRUE,
    show_col_types = FALSE
  )

  tbl_final <- tbl %>%
    left_join(locus_to_symbol, by = c("geneId" = "locus_tag")) %>%
    mutate(
      gene_label = case_when(
        !is.na(gene_symbol) & pident >= 40 & evalue <= 1e-5 ~ gene_symbol,
        TRUE ~ geneId
      )
    )

  write_csv(tbl_final, out_file)
  message("  -> wrote ", out_file)
  cat("Total", label, ":", nrow(tbl_final), "\n")
  cat("With gene symbol:", sum(!is.na(tbl_final$gene_symbol)), "\n")
  cat("With symbol (quality filter):",
      sum(tbl_final$gene_label != tbl_final$geneId, na.rm = TRUE), "\n")

  invisible(tbl_final)
}

for (i in seq_len(nrow(comparisons))) {
  annotate_and_write(comparisons$dmc_file[i], comparisons$dmc_output_file[i], "DMCs")
  annotate_and_write(comparisons$dmr_file[i], comparisons$dmr_output_file[i], "DMRs")
}
