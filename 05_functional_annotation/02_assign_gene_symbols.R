###############################################################################
# Files required for downstream analysis
###############################################################################

# Transfer the following files to the directory containing the R analysis:
#
#   blast_results.txt
#   <FEATURE_TABLE>.txt
#   <ANNOTATED_DMC_TABLE>.csv
#
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
#   - annotated_DMCs.csv
#
# Output:
#   - DMCs_with_gene_symbols.csv
###############################################################################

library(data.table)
library(tidyverse)

###############################################################################
# User input
###############################################################################

setwd("<WORKING_DIRECTORY>")

blast_file   <- "blast_results.txt"
feature_file <- "<FEATURE_TABLE>.txt"
dmc_file     <- "<ANNOTATED_DMC_TABLE>.csv"

###############################################################################
# Read and parse BLAST results
###############################################################################

blast_raw <- readLines(blast_file)

blast_parsed <- str_split(blast_raw, "\t", simplify = TRUE)

colnames(blast_parsed) <- c(
  "protein_id",
  "sseqid",
  "pident",
  "evalue",
  "stitle"
)

blast_res <- as.data.table(blast_parsed) %>%
  mutate(
    protein_id = str_trim(protein_id),
    pident = as.numeric(pident),
    evalue = as.numeric(evalue),
    gene_symbol = str_extract(stitle, "(?<=GN=)\\S+"),
    protein_name = str_trim(str_extract(stitle, "^[^O]+(?=OS=)"))
  ) %>%
  group_by(protein_id) %>%
  slice_min(evalue, n = 1, with_ties = FALSE) %>%
  ungroup()

###############################################################################
# Read feature table
###############################################################################

feature_table <- fread(
  feature_file,
  skip = 1,
  header = FALSE
)

colnames(feature_table) <- c(
  "feature",
  "class",
  "assembly",
  "assembly_unit",
  "seq_type",
  "chromosome",
  "genomic_accession",
  "start",
  "end",
  "strand",
  "product_accession",
  "non_redundant_refseq",
  "related_accession",
  "name",
  "symbol",
  "GeneID",
  "locus_tag",
  "feature_interval_length",
  "product_length",
  "attributes"
)

###############################################################################
# Build locus-to-protein mapping
###############################################################################

locus_to_protein <- feature_table %>%
  filter(
    feature == "CDS",
    product_accession != ""
  ) %>%
  select(
    locus_tag,
    protein_id = product_accession
  ) %>%
  distinct()

locus_to_symbol <- locus_to_protein %>%
  left_join(blast_res, by = "protein_id")

###############################################################################
# Read annotated DMC table
###############################################################################

methyl <- read_csv(
  dmc_file,
  quote = '""',
  trim_ws = TRUE,
  show_col_types = FALSE
)

header_names <- colnames(methyl)[1]

real_names <- unlist(
  strsplit(
    gsub('"', "", header_names),
    ","
  )
)

methyl_clean <- methyl %>%
  separate_wider_delim(
    cols = 1,
    delim = ",",
    names = real_names,
    too_many = "merge"
  ) %>%
  mutate(
    across(
      everything(),
      ~ str_remove_all(.x, '\\"')
    )
  )

###############################################################################
# Assign gene symbols
###############################################################################

methyl_final <- methyl_clean %>%
  left_join(
    locus_to_symbol,
    by = c("geneId" = "locus_tag")
  ) %>%
  mutate(
    gene_label = case_when(
      !is.na(gene_symbol) &
        pident >= 40 &
        evalue <= 1e-5 ~ gene_symbol,
      TRUE ~ geneId
    )
  )

###############################################################################
# Summary statistics
###############################################################################

cat("Total DMCs:", nrow(methyl_final), "\n")

cat(
  "Matched to UniProt:",
  sum(!is.na(methyl_final$gene_symbol)),
  "\n"
)

cat(
  "Gene symbols assigned:",
  sum(
    methyl_final$gene_label != methyl_final$geneId,
    na.rm = TRUE
  ),
  "\n"
)

###############################################################################
# Export annotated table
###############################################################################

write_csv(
  methyl_final,
  "DMCs_with_gene_symbols.csv"
)

cat("\nAnnotation complete.\n")