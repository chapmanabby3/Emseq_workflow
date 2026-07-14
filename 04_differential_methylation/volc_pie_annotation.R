################################################################
# MASTER ANALYSIS: Annotation + Volcano + DMC/DMR Pie Charts
#
# This script:
#   1. Loads DSS results (DML, DMC, DMR)
#   2. Annotates CpGs and regions using ChIPseeker
#   3. Generates volcano plots
#   4. Generates DMC/DMR pie charts
#   5. Supports multiple pairwise comparisons
################################################################

library(bsseq)
library(HDF5Array)
library(rtracklayer)
library(GenomicRanges)
library(ChIPseeker)
library(txdbmaker)
library(DSS)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(scales)

# ==========================================================
# USER INPUT
# ==========================================================

# Base directory containing DSS outputs
base_dir <- "/path/to/dss_results/DSS_pairwise_complete"

# HDF5 BSseq object directory
h5_dir <- "/path/to/dss_results/emseq_hdf5"

# Output directory
output_dir <- file.path(base_dir, "master_plots")
dir.create(output_dir, recursive = TRUE)

plots_dir <- file.path(output_dir, "plots")
dir.create(plots_dir, showWarnings = FALSE)

# Species annotation file (GTF/GFF)
gtf_file <- "genome_annotation.gtf.gz"

# If missing, download (optional template behavior)
if (!file.exists(gtf_file)) {
  stop("Please provide a valid GTF file path.")
}

txdb <- makeTxDbFromGFF(gtf_file, format = "gtf")

# ==========================================================
# COLORS
# ==========================================================

colors <- list(
  volcano = c(
    "Hyper" = "#D2691E",
    "Hypo" = "#87CEEB",
    "NotSig" = "#D3D3D3"
  ),
  pie = c(
    "Promoter"   = "#E86A33",
    "Exon"       = "#F4A261",
    "Intron"     = "#2A9D8F",
    "UTR"        = "#E9C46A",
    "Downstream" = "#A8DADC",
    "TSS"        = "#457B9D",
    "Intergenic" = "#F1FAEE",
    "Other"      = "#95A3A6"
  )
)

# ==========================================================
# COMPARISONS (EDIT THIS ONLY)
# ==========================================================

comparisons <- list(
  list(folder = "GROUP1_vs_GROUP2", label = "Comparison 1"),
  list(folder = "GROUP1_vs_GROUP3", label = "Comparison 2"),
  list(folder = "GROUP2_vs_GROUP3", label = "Comparison 3")
)

# ==========================================================
# CRITERIA SUMMARY
# ==========================================================

criteria_summary <- data.frame(
  Metric = c("DMC Criteria", "DMR Criteria"),
  Thresholds = c(
    "FDR < 0.05 | ΔMeth > 20%",
    "p < 0.01 | ΔMeth > 20%"
  )
)

write.csv(criteria_summary,
          file.path(output_dir, "criteria_summary.csv"),
          row.names = FALSE)

# ==========================================================
# LOAD BSSEQ
# ==========================================================

bsseq_obj <- as(loadHDF5SummarizedExperiment(h5_dir), "BSseq")

sample_info <- data.frame(
  id = colnames(bsseq_obj),
  condition = pData(bsseq_obj)$condition,
  sex = pData(bsseq_obj)$sex
)

pData(bsseq_obj) <- sample_info

# ==========================================================
# LOOP THROUGH COMPARISONS
# ==========================================================

for (comp in comparisons) {

  comp_folder <- comp$folder
  comp_label <- comp$label
  comp_path <- file.path(base_dir, comp_folder)

  message(paste("Processing:", comp_label))

  # -------------------------
  # Load DSS results
  # -------------------------

  dmc <- readRDS(file.path(comp_path, "02_dmc.rds"))
  dmr <- readRDS(file.path(comp_path, "03_dmr.rds"))
  dml <- readRDS(file.path(comp_path, "01_dml_test.rds"))

  if (nrow(dmc) == 0) {
    message("No significant DMCs found.")
    next
  }

  # ==========================================================
  # DMC ANNOTATION
  # ==========================================================

  dmc_gr <- GRanges(dmc$chr, IRanges(dmc$pos, width = 1))

  dmc_anno <- annotatePeak(
    dmc_gr,
    TxDb = txdb,
    tssRegion = c(-2000, 2000)
  )

  dmc_anno_df <- as.data.frame(dmc_anno)

  dmc_anno_df$feature_category <- case_when(
    grepl("Promoter", dmc_anno_df$annotation) ~ "Promoter",
    grepl("Exon", dmc_anno_df$annotation) ~ "Exon",
    grepl("Intron", dmc_anno_df$annotation) ~ "Intron",
    grepl("UTR", dmc_anno_df$annotation) ~ "UTR",
    grepl("Downstream", dmc_anno_df$annotation) ~ "Downstream",
    grepl("TSS", dmc_anno_df$annotation) ~ "TSS",
    grepl("Intergenic", dmc_anno_df$annotation) ~ "Intergenic",
    TRUE ~ "Other"
  )

  dmc_anno_df$Gene <- ifelse(
    is.na(dmc_anno_df$geneId),
    paste0(dmc_anno_df$seqnames, ":", dmc_anno_df$start),
    dmc_anno_df$geneId
  )

  # ==========================================================
  # DMR ANNOTATION
  # ==========================================================

  if (!is.null(dmr) && nrow(dmr) > 0) {

    dmr_gr <- GRanges(dmr$chr, IRanges(dmr$start, dmr$end))

    dmr_anno <- annotatePeak(
      dmr_gr,
      TxDb = txdb,
      tssRegion = c(-2000, 2000)
    )

    dmr_anno_df <- as.data.frame(dmr_anno)

    dmr_anno_df$feature_category <- case_when(
      grepl("Promoter", dmr_anno_df$annotation) ~ "Promoter",
      grepl("Exon", dmr_anno_df$annotation) ~ "Exon",
      grepl("Intron", dmr_anno_df$annotation) ~ "Intron",
      grepl("UTR", dmr_anno_df$annotation) ~ "UTR",
      grepl("Downstream", dmr_anno_df$annotation) ~ "Downstream",
      grepl("TSS", dmr_anno_df$annotation) ~ "TSS",
      grepl("Intergenic", dmr_anno_df$annotation) ~ "Intergenic",
      TRUE ~ "Other"
    )

    write.csv(dmr_anno_df,
              file.path(output_dir,
                        paste0("DMRs_annotated_", comp_folder, ".csv")),
              row.names = FALSE)
  }

  # ==========================================================
  # VOLCANO DATA PREP
  # ==========================================================

  dml_df <- as.data.frame(dml)

  dml_df$Gene <- dmc_anno_df$Gene[
    match(
      paste(dml_df$chr, dml_df$pos),
      paste(dmc_anno_df$seqnames, dmc_anno_df$start)
    )
  ]

  dml_df$direction <- "NotSig"

  dml_df$direction[dml_df$fdr < 0.05 & dml_df$diff > 0.20] <- "Hyper"
  dml_df$direction[dml_df$fdr < 0.05 & dml_df$diff < -0.20] <- "Hypo"

  top10 <- dml_df %>%
    filter(fdr < 0.05, abs(diff) > 0.20, !is.na(Gene)) %>%
    arrange(fdr) %>%
    slice_head(n = 10)

  # ==========================================================
  # VOLCANO PLOT
  # ==========================================================

  volcano <- ggplot(dml_df, aes(x = diff, y = -log10(fdr))) +
    geom_point(aes(color = direction), alpha = 0.6, size = 1.2) +
    scale_color_manual(values = colors$volcano) +
    geom_vline(xintercept = c(-0.20, 0.20), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_text_repel(data = top10, aes(label = Gene), size = 3) +
    labs(
      title = comp_label,
      x = "Methylation Difference",
      y = "-log10(FDR)"
    ) +
    theme_minimal()

  ggsave(file.path(plots_dir,
                   paste0("volcano_", comp_folder, ".png")),
         volcano, width = 9, height = 7, dpi = 300)

  # ==========================================================
  # PIE CHARTS (DMC)
  # ==========================================================

  pie_dmc <- dmc_anno_df %>%
    count(feature_category) %>%
    mutate(p = n / sum(n))

  pie_plot <- ggplot(pie_dmc, aes(x = 1, y = n, fill = feature_category)) +
    geom_col() +
    coord_polar(theta = "y") +
    scale_fill_manual(values = colors$pie) +
    theme_void()

  ggsave(file.path(plots_dir,
                   paste0("pie_DMC_", comp_folder, ".png")),
         pie_plot, width = 8, height = 6)

  message(paste("Completed:", comp_label))
}

message("ALL ANALYSES COMPLETE")
