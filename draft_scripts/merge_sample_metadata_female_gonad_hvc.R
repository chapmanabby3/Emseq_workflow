# ==============================================================================
# merge_sample_metadata_female_gonad_hvc.R
#
# Description:
# - Loads BSseq HDF5 for female gonad.
# - Constructs full sample metadata (sample, rep, condition, sex, tissue).
# - Loads DMRs (03_dmr.rds) and annotated DMR table.
# - Re-extracts methylation matrix from BSseq for those DMRs.
# - Writes:
#   1) DMR table with contrast- and sample-level metadata.
#   2) DMR methylation matrix (rows = DMRs, cols = samples).
#   3) Sample metadata table.
#
# Comparison: female gonad, high vs control (pfg_hvc)
# ==============================================================================

suppressPackageStartupMessages({
  library(bsseq)
  library(HDF5Array)
  library(GenomicRanges)
  library(dplyr)
})

# ==============================================================================
# 1. USER CONFIGURATION
# ==============================================================================

## Paths
h5_dir <- "/cluster/work/users/abchap73/R_analysis/gonad/female/dss_results/emseq_hdf5"

dmr_rds_dir <- "/cluster/work/users/abchap73/R_analysis/gonad/female/dss_results/old_results/DSS_pairwise_complete/high_vs_control"
dmr_anno_csv <- "/cluster/work/users/abchap73/borsaiasm_blast_fem_gon_high_vs_control_DMR/DMRs_annotated_pfg_hvc_with_gene_symbols.csv"

output_dir <- "/cluster/work/users/abchap73/tables_merged/gonad/female"

## Contrast labels
contrast_label   <- "pfg_hvc"  # polar female gonad, high vs control
group1_condition <- "high"
group2_condition <- "control"

## Tissue/sex constants for this analysis
sex_type   <- "female"
tissue_type <- "gonad"

# ==============================================================================
# 2. LOAD BSSEQ OBJECT
# ==============================================================================

message("Loading BSseq object from HDF5...")
se <- loadHDF5SummarizedExperiment(h5_dir)
bsseq_obj <- as(se, "BSseq")

all_sample_ids <- colnames(bsseq_obj)
message("Total samples in BSseq: ", length(all_sample_ids))

# ==============================================================================
# 3. BUILD SAMPLE METADATA (sample, rep, condition, sex, tissue)
# ==============================================================================

## You already know the exact sample vector for this analysis:
sample_info <- data.frame(
  id = c(
    "2_female_gonad_control_rep1",
    "4_female_gonad_high_rep1",
    "6_female_gonad_low_rep1",
    "8_female_gonad_control_rep2",
    "12_female_gonad_high_rep2",
    "14_female_gonad_high_rep3",
    "16_female_gonad_control_rep3",
    "18_female_gonad_low_rep3",
    "19_female_gonad_control_rep4",
    "22_female_gonad_low_rep5",
    "26_female_gonad_high_rep4",
    "27_female_gonad_control_rep5",
    "30_female_gonad_high_rep5"
  ),
  condition = c(
    "control",
    "high",
    "low",
    "control",
    "high",
    "high",
    "control",
    "low",
    "control",
    "low",
    "high",
    "control",
    "high"
  ),
  sex    = rep(sex_type, 13),
  tissue = rep(tissue_type, 13),
  stringsAsFactors = FALSE
)

rownames(sample_info) <- sample_info$id

## Extract sample number and replicate from id, e.g.:
## "2_female_gonad_control_rep1" -> sample = 2, rep = 1
parse_sample_and_rep <- function(id_vec) {
  parts <- strsplit(id_vec, "_")
  sample_num <- sapply(parts, function(p) as.integer(p[1]))
  # assume last element is "repN"
  rep_num <- sapply(parts, function(p) {
    last <- p[length(p)]
    if (grepl("^rep", last)) {
      as.integer(sub("rep", "", last))
    } else {
      NA_integer_
    }
  })
  data.frame(sample = sample_num, rep = rep_num, stringsAsFactors = FALSE)
}

samp_rep_df <- parse_sample_and_rep(sample_info$id)
sample_info$sample <- samp_rep_df$sample
sample_info$rep    <- samp_rep_df$rep

## Ensure BSseq pData matches this sample_info
pData(bsseq_obj) <- sample_info

message("Sample metadata constructed: ", nrow(sample_info), " samples.")

# ==============================================================================
# 4. LOAD DMRs (RDS) AND ANNOTATED DMR TABLE (CSV)
# ==============================================================================

message("Loading DMRs (03_dmr.rds)...")
dmr_rds_path <- file.path(dmr_rds_dir, "03_dmr.rds")
dmr <- readRDS(dmr_rds_path)

if (nrow(dmr) == 0) {
  stop("No DMRs found in 03_dmr.rds. Check path and contents.")
}

message("DMRs loaded: ", nrow(dmr))

message("Loading annotated DMR table: ", dmr_anno_csv)
dmr_anno <- read.csv(dmr_anno_csv, stringsAsFactors = FALSE)

message("Annotated DMR table loaded: ", nrow(dmr_anno), " rows.")

## Ensure that the DMRs in the CSV correspond to the RDS by dmr_id or coordinates.
## We'll assume they match row-for-row after sorting by chr/start/end.
## If you know they are already aligned, you can skip reordering.

order_dmr <- function(df) {
  if ("chr" %in% colnames(df)) {
    chr_col <- "chr"
  } else if ("seqnames" %in% colnames(df)) {
    chr_col <- "seqnames"
  } else {
    stop("No chr or seqnames column found in DMR table.")
  }

  start_col <- if ("start" %in% colnames(df)) "start" else "start.x"
  end_col   <- if ("end" %in% colnames(df)) "end" else "end.x"

  df[order(df[[chr_col]], df[[start_col]], df[[end_col]]), , drop = FALSE]
}

dmr_sorted <- order_dmr(dmr)
dmr_anno_sorted <- order_dmr(dmr_anno)

if (nrow(dmr_sorted) != nrow(dmr_anno_sorted)) {
  warning(
    "Row counts differ between DMR RDS and annotated CSV (",
    nrow(dmr_sorted), " vs ", nrow(dmr_anno_sorted),
    "). Merging will be by dmr_id if available, otherwise by chr/start/end."
  )
}

# ==============================================================================
# 5. CREATE DMR GRanges AND EXTRACT METHYLATION MATRIX FROM BSSEQ
# ==============================================================================

## Ensure DMR has dmr_id; if not, create one.
if (!"dmr_id" %in% colnames(dmr_sorted)) {
  dmr_sorted$dmr_id <- paste0(
    "DMR_",
    seq_len(nrow(dmr_sorted))
  )
}

dmr_gr <- GRanges(
  seqnames = dmr_sorted$chr,
  ranges   = IRanges(start = dmr_sorted$start, end = dmr_sorted$end),
  dmr_id   = dmr_sorted$dmr_id
)

message("Extracting methylation matrix from BSseq for ", length(dmr_gr), " DMRs...")

dmr_meth_mat <- getMeth(
  bsseq_obj,
  type = "raw",
  regions = dmr_gr,
  what = "perRegion"
)

rownames(dmr_meth_mat) <- dmr_gr$dmr_id
colnames(dmr_meth_mat) <- colnames(bsseq_obj)  # sample IDs

message("Methylation matrix dimensions: ",
        nrow(dmr_meth_mat), " DMRs x ", ncol(dmr_meth_mat), " samples.")

# ==============================================================================
# 6. BUILD CONTRAST-LEVEL METADATA
# ==============================================================================

group1_samples <- sample_info$id[sample_info$condition == group1_condition]
group2_samples <- sample_info$id[sample_info$condition == group2_condition]

contrast_meta <- data.frame(
  contrast        = contrast_label,
  group1_condition = group1_condition,
  group2_condition = group2_condition,
  group1_samples   = paste(group1_samples, collapse = ";"),
  group2_samples   = paste(group2_samples, collapse = ";"),
  sex              = sex_type,
  tissue           = tissue_type,
  stringsAsFactors = FALSE
)

# ==============================================================================
# 7. MERGE ANNOTATED DMR TABLE WITH CONTRAST METADATA
# ==============================================================================

## Add contrast metadata columns to every row of dmr_anno
dmr_anno_with_meta <- cbind(
  contrast_meta,
  dmr_anno_sorted
)

## Optionally, if you want to ensure the RDS-based DMR info matches the CSV,
## you could merge by dmr_id or coordinates. For now, we assume the CSV is
## your final annotated table and we just prepend contrast metadata.

# ==============================================================================
# 8. SAVE OUTPUTS
# ==============================================================================

## 1) DMR table with sample/contrast metadata
dmr_out_csv <- file.path(
  output_dir,
  paste0("DMRs_annotated_", contrast_label, "_with_sample_metadata.csv")
)

write.csv(
  dmr_anno_with_meta,
  dmr_out_csv,
  row.names = FALSE,
  quote = FALSE
)
message("Saved DMR table with metadata: ", dmr_out_csv)

## 2) DMR methylation matrix
meth_mat_out_csv <- file.path(
  output_dir,
  paste0("DMR_methylation_matrix_", contrast_label, ".csv")
)

write.csv(
  dmr_meth_mat,
  meth_mat_out_csv,
  quote = FALSE
)
message("Saved DMR methylation matrix: ", meth_mat_out_csv)

## 3) Sample metadata table
sample_meta_out_csv <- file.path(
  output_dir,
  paste0("sample_metadata_", contrast_label, ".csv")
)

sample_meta_out <- sample_info[, c("id", "sample", "rep", "condition", "sex", "tissue")]
colnames(sample_meta_out)[1] <- "sample_id"

write.csv(
  sample_meta_out,
  sample_meta_out_csv,
  row.names = FALSE,
  quote = FALSE
)
message("Saved sample metadata: ", sample_meta_out_csv)

message("Done! Files written to: ", output_dir)
