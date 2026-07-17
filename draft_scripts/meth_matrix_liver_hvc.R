# ==============================================================================
# DRAFT SCRIPT: meth_matrix_liver_hvc.R
# 
# Description and Functions:
# 1. Loads an out-of-memory HDF5-backed BSseq object containing raw methylation calls.
# 2. Generates consistent DMR IDs (e.g., "DMR_1") prior to any genomic coordinate filtering.
# 3. Validates and filters DMR chromosomes against the seqlevels of a provided Polar cod GTF.
# 4. Annotates valid DMRs with genomic features (Promoter, Exon, Intron, etc.) via ChIPseeker.
# 5. Joins the annotations back to the full set of DMRs. Unmappable scaffolds are flagged as "Unannotated".
# 6. Extracts a matrix of average region-level methylation scores for each sample using getMeth().
# 7. Saves the annotated DMR table, the regional methylation matrix, and a distribution pie chart.
# 
# Note: This is an exploratory/draft script and may not yet be modularized for the main pipeline.
# ==============================================================================

library(bsseq)
library(HDF5Array)
library(rtracklayer)
library(GenomicRanges)
library(ChIPseeker)
library(txdbmaker)
library(DSS)
library(dplyr)
library(ggplot2)

# ============================================================
# PATHS
# ============================================================

base_dir <- "/cluster/work/users/abchap73/R_analysis/liver/female/dss_results/DSS_pairwise_complete/high_vs_control_CG3"
h5_dir <- "/cluster/work/users/abchap73/R_analysis/liver/female/dss_results/emseq_hdf5"

output_dir <- file.path(base_dir, "master_plots_final")
plots_dir  <- file.path(output_dir, "plots")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir, showWarnings = FALSE)

comp_label <- "High vs Unexposed"

message("Processing: ", comp_label)

# ============================================================
# LOAD BSSEQ OBJECT
# ============================================================

bsseq_obj <- as(loadHDF5SummarizedExperiment(h5_dir), "BSseq")

sample_info <- data.frame(
    id = colnames(bsseq_obj),
    condition = pData(bsseq_obj)$condition,
    sex = pData(bsseq_obj)$sex
)

pData(bsseq_obj) <- sample_info

# ============================================================
# LOAD GTF
# ============================================================

gtf_file <- "/cluster/work/users/abchap73/GCA_964260565.1_borsaiasm_genomic.gtf.gz"

txdb <- makeTxDbFromGFF(gtf_file, format = "gtf")

# ============================================================
# COLORS
# ============================================================

pie_colors <- c(
    Promoter    = "#E86A33",
    Exon        = "#F4A261",
    Intron      = "#2A9D8F",
    UTR         = "#E9C46A",
    Downstream  = "#A8DADC",
    TSS         = "#457B9D",
    Intergenic  = "#F1FAEE",
    Other       = "#95A3A6",
    Unannotated = "#CCCCCC"
)

# ============================================================
# LOAD DMRs
# ============================================================

dmr <- readRDS(file.path(base_dir, "03_dmr.rds"))

if (nrow(dmr) == 0) {

    stop("No DMRs found.")

}

message("DMRs found: ", nrow(dmr))

# ============================================================
# ANNOTATE DMRs
# ============================================================

dmr$dmr_id <- paste0("DMR_", seq_len(nrow(dmr)))

dmr_gr_full <- GRanges(
    seqnames = dmr$chr,
    ranges = IRanges(
        start = dmr$start,
        end   = dmr$end
    ),
    dmr_id = dmr$dmr_id
)

seqlevels_in_anno <- seqlevels(txdb)
keep_for_anno <- as.character(seqnames(dmr_gr_full)) %in% seqlevels_in_anno
dmr_gr_anno <- dmr_gr_full[keep_for_anno]

message("DMRs total: ", length(dmr_gr_full), " | kept for annotation: ", length(dmr_gr_anno), " | dropped (no matching chr in GTF): ", sum(!keep_for_anno))

if (length(dmr_gr_anno) > 0) {
    dmr_anno <- annotatePeak(
        dmr_gr_anno,
        TxDb = txdb,
        tssRegion = c(-2000, 2000)
    )

    dmr_anno_df <- as.data.frame(dmr_anno)
    dmr_anno_df$dmr_id <- dmr_gr_anno$dmr_id

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
} else {
    dmr_anno_df <- data.frame(dmr_id = character(0), feature_category = character(0))
}

# ============================================================
# MERGE DSS + ANNOTATION
# ============================================================

dmr_sig <- left_join(
    dmr,
    dmr_anno_df,
    by = "dmr_id"
)

# For unannotated scaffolds, set category to "Unannotated"
dmr_sig$feature_category[is.na(dmr_sig$feature_category)] <- "Unannotated"

# ============================================================
# EXTRACT REGION-LEVEL METHYLATION
# ============================================================

dmr_methylation_matrix <- getMeth(
    bsseq_obj,
    type = "raw",
    regions = dmr_gr_full,
    what = "perRegion"
)

rownames(dmr_methylation_matrix) <- dmr_gr_full$dmr_id

# ============================================================
# SAVE TABLES
# ============================================================

write.csv(
    dmr_sig,
    file.path(output_dir, "DMRs_annotated.csv"),
    row.names = FALSE
)

write.csv(
    dmr_methylation_matrix,
    file.path(output_dir, "DMR_methylation_matrix.csv")
)

# ============================================================
# DMR PIE CHART
# ============================================================

pie_data <- dmr_sig %>%
    count(feature_category, sort = TRUE) %>%
    mutate(
        percentage = round(100 * n / sum(n), 1)
    )

pie_plot <- ggplot(
    pie_data,
    aes(x = "", y = n, fill = feature_category)
) +
    geom_bar(
        stat = "identity",
        width = 1,
        color = "white"
    ) +
    coord_polar("y") +
    geom_text(
        aes(label = paste0(n, "\n(", percentage, "%)")),
        position = position_stack(vjust = 0.5),
        color = "black",
        size = 4,
        fontface = "bold"
    ) +
    scale_fill_manual(values = pie_colors) +
    labs(
        title = paste("Female Liver:", comp_label),
        subtitle = "DMR Genomic Feature Distribution"
    ) +
    theme_void(base_size = 12) +
    theme(
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "right"
    )

ggsave(
    file.path(plots_dir, "DMR_pie_chart.png"),
    pie_plot,
    width = 10,
    height = 8,
    dpi = 300,
    bg = "white"
)

message("Done!")
message("Annotated DMR table saved.")
message("DMR methylation matrix saved.")
message("DMR pie chart saved.")
message("Output directory: ", output_dir)
