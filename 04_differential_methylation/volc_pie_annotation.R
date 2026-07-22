# ==========================================================
# USER INPUT
# ==========================================================

base_dir <- "/path/to/DSS_pairwise_complete"
h5_dir <- "/path/to/emseq_hdf5"

gtf_file <- "your_annotation.gtf.gz"

organism_label <- "Your species / tissue"

comparisons <- list(
  list(folder="comparison1", label="Group A vs Group B"),
  list(folder="comparison2", label="Group A vs Group C"),
  list(folder="comparison3", label="Group B vs Group C")
)
################################################################
# FINAL MASTER TEMPLATE: Annotation + Volcano (Dual) + DMC/DMR Pie Charts
#
# This script:
#   1. Loads DSS results (DML, DMC, DMR)
#   2. Annotates CpGs and regions using ChIPseeker
#   3. Generates volcano plots (full range + poster view)
#   4. Generates DMC and DMR genomic feature pie charts
#   5. Supports multiple pairwise DSS comparisons
#
# Replace only the USER INPUT section.
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

# Base directory containing DSS pairwise outputs
#
# Expected structure:
#
# base_dir/
# ├── comparison1/
# │   ├── 01_dml_test.rds
# │   ├── 02_dmc.rds
# │   └── 03_dmr.rds
# ├── comparison2/
# └── comparison3/
#

base_dir <- "/path/to/DSS_pairwise_complete"


# HDF5 BSseq object directory

h5_dir <- "/path/to/emseq_hdf5"


# Output directory

output_dir <- file.path(base_dir, "master_plots_final")

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


plots_dir <- file.path(output_dir, "plots")

dir.create(
  plots_dir,
  showWarnings = FALSE
)


# Genome annotation file

gtf_file <- "genome_annotation.gtf.gz"


# If annotation is missing, stop

if (!file.exists(gtf_file)) {
  stop(
    "Annotation file not found. Please provide a valid GTF/GFF file."
  )
}


# Create TxDb object

txdb <- makeTxDbFromGFF(
  gtf_file,
  format = "gtf"
)



# Species/tissue label used in figures

analysis_label <- "Species Tissue"



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
# DSS COMPARISONS
# ==========================================================


# Folder names must match your DSS output directories

comparisons <- list(

  list(
    folder = "comparison1",
    label = "Group A vs Group B"
  ),

  list(
    folder = "comparison2",
    label = "Group A vs Group C"
  ),

  list(
    folder = "comparison3",
    label = "Group B vs Group C"
  )

)



# ==========================================================
# CRITERIA SUMMARY
# ==========================================================


criteria_summary <- data.frame(

  Metric = c(
    "DMC Criteria",
    "DMR Criteria"
  ),


  Thresholds = c(
    "FDR<0.05, |Delta Meth|>20%",
    "p<0.01, |Delta Meth|>20%"
  )

)


write.csv(
  criteria_summary,
  file.path(
    output_dir,
    "criteria_summary.csv"
  ),
  row.names = FALSE
)



# ==========================================================
# LOAD BSSEQ OBJECT
# ==========================================================


bsseq_obj <- as(
  loadHDF5SummarizedExperiment(h5_dir),
  "BSseq"
)



sample_info <- data.frame(

  id = colnames(bsseq_obj),

  condition = pData(bsseq_obj)$condition,

  sex = pData(bsseq_obj)$sex

)


pData(bsseq_obj) <- sample_info



message(
  "MASTER ANALYSIS STARTED"
)



# ==========================================================
# PROCESS EACH COMPARISON
# ==========================================================


for (comp in comparisons) {


  comp_folder <- comp$folder

  comp_label <- comp$label


  comp_path <- file.path(
    base_dir,
    comp_folder
  )


  message(
    paste(
      "Processing:",
      comp_label
    )
  )



  # ========================================================
  # LOAD DSS FILES
  # ========================================================


  dmc <- readRDS(
    file.path(
      comp_path,
      "02_dmc.rds"
    )
  )


  dmr <- readRDS(
    file.path(
      comp_path,
      "03_dmr.rds"
    )
  )


  dml <- readRDS(
    file.path(
      comp_path,
      "01_dml_test.rds"
    )
  )



  if (nrow(dmc) == 0) {

    message(
      "No DMCs found"
    )

    next

  }



  message(
    paste(
      "DMCs:",
      nrow(dmc),
      "DMRs:",
      nrow(dmr)
    )
  )



  # ========================================================
  # ANNOTATE DMCs
  # ========================================================



  dmc_gr <- GRanges(

    dmc$chr,

    IRanges(
      dmc$pos,
      width = 1
    )

  )



  dmc_anno <- annotatePeak(

    dmc_gr,

    TxDb = txdb,

    tssRegion = c(
      -2000,
      2000
    )

  )



  dmc_anno_df <- as.data.frame(
    dmc_anno
  )



  dmc_anno_df$feature_category <- case_when(

    grepl(
      "Promoter",
      dmc_anno_df$annotation
    ) ~ "Promoter",

    grepl(
      "Exon",
      dmc_anno_df$annotation
    ) ~ "Exon",

    grepl(
      "Intron",
      dmc_anno_df$annotation
    ) ~ "Intron",

    grepl(
      "UTR",
      dmc_anno_df$annotation
    ) ~ "UTR",

    grepl(
      "Downstream",
      dmc_anno_df$annotation
    ) ~ "Downstream",

    grepl(
      "TSS",
      dmc_anno_df$annotation
    ) ~ "TSS",

    grepl(
      "Intergenic",
      dmc_anno_df$annotation
    ) ~ "Intergenic",

    TRUE ~ "Other"

  )



  dmc_anno_df$Gene <- ifelse(

    is.na(
      dmc_anno_df$geneId
    ),

    paste0(
      dmc_anno_df$seqnames,
      ":",
      dmc_anno_df$start
    ),

    dmc_anno_df$geneId

  )



  # ========================================================
  # ANNOTATE DMRs
  # ========================================================


  if (
    !is.null(dmr) &&
    nrow(dmr) > 0
  ) {


    dmr_gr <- GRanges(

      dmr$chr,

      IRanges(
        dmr$start,
        dmr$end
      )

    )


    dmr_anno <- annotatePeak(

      dmr_gr,

      TxDb = txdb,

      tssRegion = c(
        -2000,
        2000
      )

    )


    dmr_anno_df <- as.data.frame(
      dmr_anno
    )


    dmr_anno_df$feature_category <- case_when(

      grepl(
        "Promoter",
        dmr_anno_df$annotation
      ) ~ "Promoter",

      grepl(
        "Exon",
        dmr_anno_df$annotation
      ) ~ "Exon",

      grepl(
        "Intron",
        dmr_anno_df$annotation
      ) ~ "Intron",

      grepl(
        "UTR",
        dmr_anno_df$annotation
      ) ~ "UTR",

      grepl(
        "Downstream",
        dmr_anno_df$annotation
      ) ~ "Downstream",

      grepl(
        "TSS",
        dmr_anno_df$annotation
      ) ~ "TSS",

      grepl(
        "Intergenic",
        dmr_anno_df$annotation
      ) ~ "Intergenic",

      TRUE ~ "Other"

    )


    dmr_sig <- dmr_anno_df


    write.csv(

      dmr_sig,

      file.path(
        output_dir,
        paste0(
          "DMRs_annotated_",
          comp_folder,
          ".csv"
        )
      ),

      row.names = FALSE

    )

  }
    # ========================================================
  # TOP 10 LABELS
  # ========================================================

  top10 <- dml_df %>%
    filter(
      fdr < 0.05,
      abs(diff) > delta_threshold,
      !is.na(Gene),
      Gene != ""
    ) %>%
    arrange(fdr) %>%
    slice_head(n = 10)


  # ========================================================
  # VOLCANO FULL RANGE
  # ========================================================

  volcano_full <- ggplot(
    dml_df,
    aes(
      x = diff,
      y = -log10(fdr)
    )
  ) +
    geom_point(
      aes(color = direction),
      alpha = 0.6,
      size = 1.2
    ) +
    scale_color_manual(
      values = colors$volcano
    ) +
    geom_vline(
      xintercept = c(-delta_threshold, delta_threshold),
      linetype = "dashed",
      color = "gray40",
      linewidth = 0.8
    ) +
    geom_hline(
      yintercept = -log10(fdr_threshold),
      linetype = "dashed",
      color = "gray40",
      linewidth = 0.8
    ) +
    geom_text_repel(
      data = top10,
      aes(label = Gene),
      size = 3.2,
      max.overlaps = Inf,
      box.padding = 0.3,
      point.padding = 0.3,
      segment.size = 0.3
    ) +
    labs(
      title = comp_label,
      subtitle = paste0(
        "Full range + Top 10 genes (FDR<",
        fdr_threshold,
        ", |ΔMeth|>",
        delta_threshold * 100,
        "%)"
      ),
      x = "Methylation Difference",
      y = "-log10(FDR)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 14
      ),
      plot.subtitle = element_text(
        hjust = 0.5,
        size = 11,
        color = "gray50"
      ),
      legend.position = "bottom"
    )


  ggsave(
    file.path(
      plots_dir,
      paste0(
        "volcano_full_",
        gsub(" ", "_", comp_label),
        ".png"
      )
    ),
    volcano_full,
    width = 9,
    height = 7,
    dpi = 300,
    bg = "white"
  )


  # ========================================================
  # VOLCANO POSTER VERSION
  # ========================================================

  x_lim <- max(
    abs(dml_df$diff),
    na.rm = TRUE
  ) * 1.05


  y_lim <- max(
    -log10(
      pmin(
        dml_df$fdr[dml_df$fdr > 0],
        1e-10
      )
    ),
    na.rm = TRUE
  ) * 1.05


  volcano_poster <- volcano_full +
    xlim(-x_lim, x_lim) +
    ylim(0, y_lim) +
    labs(
      subtitle = paste0(
        "Poster view + Top 10 genes (FDR<",
        fdr_threshold,
        ", |ΔMeth|>",
        delta_threshold * 100,
        "%)"
      )
    )


  ggsave(
    file.path(
      plots_dir,
      paste0(
        "volcano_poster_",
        gsub(" ", "_", comp_label),
        ".png"
      )
    ),
    volcano_poster,
    width = 9,
    height = 7,
    dpi = 300,
    bg = "white"
  )


  # ========================================================
  # DMC PIE CHART
  # ========================================================

  pie_data_dmc <- dmc_anno_df %>%
    count(feature_category, sort = TRUE) %>%
    mutate(
      percentage = round(100 * n / sum(n), 1),
      cumulative = cumsum(n),
      midpoint = cumulative - n / 2
    )


  pie_dmc <- ggplot(
    pie_data_dmc,
    aes(
      x = 1,
      y = n,
      fill = feature_category
    )
  ) +
    geom_col(
      color = "white",
      linewidth = 0.8
    ) +
    coord_polar(theta = "y") +
    geom_label_repel(
      aes(
        y = midpoint,
        label = paste0(
          n,
          "\n(",
          percentage,
          "%)"
        )
      ),
      color = "black",
      fill = "white",
      fontface = "bold",
      size = 4,
      show.legend = FALSE,
      min.segment.length = 0,
      max.overlaps = Inf,
      seed = 123
    ) +
    scale_fill_manual(
      values = colors$pie
    ) +
    xlim(c(0.5, 2.2)) +
    labs(
      title = comp_label,
      subtitle = paste0(
        "DMC Genomic Feature Distribution (FDR<",
        fdr_threshold,
        ", |ΔMeth|>",
        delta_threshold * 100,
        "%)"
      )
    ) +
    theme_void(base_size = 12) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 14
      ),
      plot.subtitle = element_text(
        hjust = 0.5,
        size = 11,
        color = "gray50"
      ),
      legend.position = "right"
    )


  ggsave(
    file.path(
      plots_dir,
      paste0(
        "pie_DMCs_",
        gsub(" ", "_", comp_label),
        ".png"
      )
    ),
    pie_dmc,
    width = 10,
    height = 8,
    dpi = 300,
    bg = "white"
  )


  # ========================================================
  # DMR PIE CHART
  # ========================================================

  if (!is.null(dmr) &&
      nrow(dmr) > 0 &&
      nrow(dmr_sig) > 0) {


    pie_data_dmr <- dmr_sig %>%
      count(feature_category, sort = TRUE) %>%
      mutate(
        percentage = round(100 * n / sum(n), 1),
        cumulative = cumsum(n),
        midpoint = cumulative - n / 2
      )


    pie_dmr <- ggplot(
      pie_data_dmr,
      aes(
        x = 1,
        y = n,
        fill = feature_category
      )
    ) +
      geom_col(
        color = "white",
        linewidth = 0.8
      ) +
      coord_polar(theta = "y") +
      geom_label_repel(
        aes(
          y = midpoint,
          label = paste0(
            n,
            "\n(",
            percentage,
            "%)"
          )
        ),
        color = "black",
        fill = "white",
        fontface = "bold",
        size = 4,
        show.legend = FALSE,
        min.segment.length = 0,
        max.overlaps = Inf,
        seed = 123
      ) +
      scale_fill_manual(
        values = colors$pie
      ) +
      xlim(c(0.5, 2.2)) +
      labs(
        title = comp_label,
        subtitle = paste0(
          "DMR Genomic Feature Distribution (p<",
          dmr_p_threshold,
          ", |ΔMeth|>",
          delta_threshold * 100,
          "%)"
        )
      ) +
      theme_void(base_size = 12) +
      theme(
        plot.title = element_text(
          hjust = 0.5,
          face = "bold",
          size = 14
        ),
        plot.subtitle = element_text(
          hjust = 0.5,
          size = 11,
          color = "gray50"
        ),
        legend.position = "right"
      )


    ggsave(
      file.path(
        plots_dir,
        paste0(
          "pie_DMRs_",
          gsub(" ", "_", comp_label),
          ".png"
        )
      ),
      pie_dmr,
      width = 10,
      height = 8,
      dpi = 300,
      bg = "white"
    )

    message("DMR pie chart saved")
  }


  # ========================================================
  # SAVE ANNOTATED TABLES
  # ========================================================

  write.csv(
    dmc_anno_df,
    file.path(
      output_dir,
      paste0(
        "DMCs_annotated_",
        comp_folder,
        ".csv"
      )
    ),
    row.names = FALSE
  )


  message(
    paste(
      "Completed:",
      comp_label
    )
  )

}


# ==========================================================
# FINISHED
# ==========================================================

message("COMPLETE!")
message("Output tables:")
message(output_dir)

message("Output plots:")
message(plots_dir)
