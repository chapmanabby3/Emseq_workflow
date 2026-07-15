###############################################################
## DSS Pairwise Differential Methylation Analysis (TEMPLATE)
###############################################################

###############################################################
## USER INPUT
###############################################################

# Biological variables (for metadata only)
SEX <- "female"
TISSUE <- "liver"

# Pairwise comparison
GROUP1_COND <- "high"
GROUP2_COND <- "control"
COMP_DIR_NAME <- "GROUP1_vs_GROUP2"

# Output directory (EDIT THIS PATH)
output_dir <- "/path/to/R_analysis/dss_results"

# Sample information (REPLACE WITH YOUR OWN SAMPLES)
sample_info <- data.frame(
  id = c(
    "sample1",
    "sample2",
    "sample3",
    "sample4"
  ),
  condition = c(
    "high",
    "high",
    "control",
    "control"
  ),
  sex = rep(SEX, 4)
)

###############################################################
## END USER INPUT
###############################################################

library(BiocParallel)
library(bsseq)
library(HDF5Array)
library(DSS)
library(dplyr)
library(S4Vectors)
library(GenomicRanges)

#-------------------------------------------------------------
# Output setup
#-------------------------------------------------------------

sub_dir <- file.path(output_dir, "DSS_pairwise_complete")
comp_dir <- file.path(sub_dir, COMP_DIR_NAME)

ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
  }
}

ensure_dir(comp_dir)

#-------------------------------------------------------------
# Load BSseq object
#-------------------------------------------------------------

message("Importing BSseq object...")

h5_bsseq_save_dir <- file.path(output_dir, "emseq_hdf5")

tryCatch({
  se <- loadHDF5SummarizedExperiment(h5_bsseq_save_dir)
  bsseq_obj <- as(se, "BSseq")
  message("BSseq object imported successfully.")
}, error = function(e) {
  stop("Failed to load BSseq object: ", e$message)
})

#-------------------------------------------------------------
# Assign metadata
#-------------------------------------------------------------

rownames(sample_info) <- sample_info$id
pData(bsseq_obj) <- sample_info

#-------------------------------------------------------------
# Coverage filtering
#-------------------------------------------------------------

message("Applying coverage filter (>= 5)...")

bs_subset <- bsseq_obj[, sample_info$id]

bs_filtered <- bs_subset[
  rowMeans(getCoverage(bs_subset, type = "Cov"), na.rm = TRUE) >= 5,
]

rm(bs_subset, bsseq_obj)
gc(verbose = FALSE)

message(paste("Filtered CpGs:", nrow(bs_filtered)))

#-------------------------------------------------------------
# Define groups
#-------------------------------------------------------------

group1_samples <- sample_info$id[sample_info$condition == GROUP1_COND]
group2_samples <- sample_info$id[sample_info$condition == GROUP2_COND]

message(paste("Running:", GROUP1_COND, "vs", GROUP2_COND))
message(paste("Group1 samples:", length(group1_samples)))
message(paste("Group2 samples:", length(group2_samples)))

#-------------------------------------------------------------
# DSS analysis
#-------------------------------------------------------------

tryCatch({

  message("Running DMLtest by chromosome...")

  chrs <- seqlevels(bs_filtered)
  dml_test_list <- list()

  for (chr in chrs) {

    keep_idx <- which(as.character(seqnames(bs_filtered)) == chr)

    if (length(keep_idx) == 0) next

    bs_chr <- bs_filtered[keep_idx, ]
    cov_matrix <- getCoverage(bs_chr, type = "Cov")

    g1_has_data <- sum(rowSums(cov_matrix[, group1_samples, drop = FALSE] > 0)) > 0
    g2_has_data <- sum(rowSums(cov_matrix[, group2_samples, drop = FALSE] > 0)) > 0

    if (g1_has_data && g2_has_data) {

      message(paste("Processing chromosome:", chr))

      try({

        tmp_dml <- DMLtest(
          bs_chr,
          group1 = group1_samples,
          group2 = group2_samples
        )

        if (!is.null(tmp_dml)) {
          dml_test_list[[chr]] <- tmp_dml
        }

      }, silent = TRUE)
    }

    rm(bs_chr, keep_idx, cov_matrix)
    gc(verbose = FALSE)
  }

  #-----------------------------------------------------------
  # Combine results
  #-----------------------------------------------------------

  dml_test <- do.call(rbind, unname(dml_test_list))

  rm(dml_test_list)
  gc(verbose = FALSE)

  if (is.null(dml_test) || nrow(dml_test) == 0) {
    stop("DMLtest resulted in 0 valid comparisons.")
  }

  #-----------------------------------------------------------
  # Call DMCs and DMRs
  # 3 DMC minimum added for DMRs #
  #-----------------------------------------------------------

  dmc <- callDML(dml_test, delta = 0.20, p.threshold = 0.05)
  dmr <- callDMR(dml_test, delta = 0.20, p.threshold = 0.01, minlen = 100, minCG = 3)

  #-----------------------------------------------------------
  # Save outputs
  #-----------------------------------------------------------

  saveRDS(dml_test, file = file.path(comp_dir, "01_dml_test.rds"))

  if (!is.null(dmc)) {
    saveRDS(dmc, file = file.path(comp_dir, "02_dmc.rds"))
  }

  if (!is.null(dmr)) {
    saveRDS(dmr, file = file.path(comp_dir, "03_dmr.rds"))
  }

  write.csv(as.data.frame(dml_test),
            file.path(comp_dir, "01_dml_test.csv"),
            row.names = FALSE)

  if (!is.null(dmc)) {
    write.csv(as.data.frame(dmc),
              file.path(comp_dir, "02_dmc.csv"),
              row.names = FALSE)
  }

  if (!is.null(dmr)) {
    write.csv(as.data.frame(dmr),
              file.path(comp_dir, "03_dmr.csv"),
              row.names = FALSE)
  }

}, error = function(e) {
  message(paste("CRITICAL ERROR:", e$message))
})
