# ==============================================================================
# DRAFT SCRIPT: dss_female_liver_high_vs_control_CG3.R
# 
# Description and Functions:
# 1. Bypasses the computationally expensive DMLtest function by directly loading 
#    a previously calculated DML test object (01_dml_test.rds).
# 2. Calls Differentially Methylated Loci (DMLs) with standard parameters (delta=0.2, p<0.05).
# 3. Calls Differentially Methylated Regions (DMRs) using strict filtering parameters:
#    - minlen = 100 (Minimum length of 100 bp)
#    - minCG = 3 (Minimum of 3 CpG sites per region)
# 4. Exports both the DML and DMR results as .rds objects and .csv tables for 
#    downstream plotting and annotation.
# 
# Note: This is an exploratory/draft script and may not yet be modularized for the main pipeline.
# ==============================================================================

library(DSS)

input_dir <- "/cluster/work/users/abchap73/R_analysis/liver/female/dss_results/DSS_pairwise_complete/high_vs_control"
output_dir <- "/cluster/work/users/abchap73/R_analysis/liver/female/dss_results/DSS_pairwise_complete/high_vs_control_CG3"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  message(paste("Created directory:", output_dir))
}

message("Loading existing DML test results...")
dml_test <- readRDS(file.path(input_dir, "01_dml_test.rds"))

message("Calling DMLs...")
dmc <- callDML(dml_test, delta = 0.20, p.threshold = 0.05)

message("Calling DMRs with minlen = 100, minCG = 3...")
dmr <- callDMR(dml_test, delta = 0.20, p.threshold = 0.01, minlen = 100, minCG = 3)

message("Saving results...")
if (!is.null(dmc)) {
  saveRDS(dmc, file = file.path(output_dir, "02_dmc.rds"))
  write.csv(as.data.frame(dmc), file = file.path(output_dir, "02_dmc.csv"), row.names = FALSE)
}

if (!is.null(dmr)) {
  saveRDS(dmr, file = file.path(output_dir, "03_dmr.rds"))
  write.csv(as.data.frame(dmr), file = file.path(output_dir, "03_dmr.csv"), row.names = FALSE)
}

message("Done!")
