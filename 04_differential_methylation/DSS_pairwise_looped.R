########################################################################################
### DSS Pairwise Template: Looped Comparisons, 3 comparisons ###
########################################################################################


library(BiocParallel)
library(bsseq)
library(HDF5Array)
library(DSS)
library(dplyr)
library(S4Vectors)
library(GenomicRanges)


# --- Main Analysis Function ---
run_dss_analysis <- function(bsseq_obj, sample_info, output_dir, multi_factor = FALSE) {
  
  # Helper function to ensure directory exists
  ensure_dir <- function(path) {
    if (!dir.exists(path)) {
      dir.create(path, recursive = TRUE)
      message(paste("Created directory:", path))
    }
  }
  
  ensure_dir(output_dir)
  saved_files <- list()
  
  if (multi_factor) {
    # --- Multi-factor analysis: condition + sex ---
    message("Starting multi-factor DSS analysis (condition + sex)...")
    sub_dir <- file.path(output_dir, "DSSmulticonditionsex")
    ensure_dir(sub_dir)
    
    design <- sample_info
    BSobj.fit <- DSS::DMLfit.multiFactor(bsseq_obj, design = design, formula = ~condition + sex)
    dml_results_multi <- DSS::DMLtest.multiFactor(BSobj.fit, coef = "conditioncontinuous")
    
    output_path <- file.path(sub_dir, "dml_results_multifactor.rds")
    saveRDS(dml_results_multi, file = output_path)
    saved_files$multifactor_results <- output_path
    
  } else {
    # --- Two-group comparisons (split by sex) ---
    message("Starting explicitly defined two-group DSS analyses...")
    sub_dir <- file.path(output_dir, "DSS_Pairwise_Comparisons")
    ensure_dir(sub_dir)
    
    # -----------------------------------------------------------------
    # INTERNAL FUNCTION: Performs a specific comparison based on arguments.
    # -----------------------------------------------------------------
    run_specific_comparison <- function(bs_obj, metadata, sex_group, g1_cond, g2_cond) {
      comp_name <- paste0(g1_cond, "_vs_", g2_cond)
      message(paste("\n---> Starting comparison:", comp_name, "for", sex_group, "samples..."))
      
      # Identify samples for the two groups
      group1_samples <- metadata$id[metadata$condition == g1_cond]
      group2_samples <- metadata$id[metadata$condition == g2_cond]
      
      # Safety check: Do we have samples for both groups?
      if(length(group1_samples) == 0 | length(group2_samples) == 0) {
        message(paste("Skipping", comp_name, "- missing samples."))
        return(NULL)
      }
      
      # IMPORTANT: Subset BSseq object ONLY for the samples in these two groups
      # This saves RAM and ensures accurate calculations
      samples_to_keep <- c(group1_samples, group2_samples)
      bs_subset <- bs_obj[, samples_to_keep]
      
      message(paste("Splitting by chromosome to save memory..."))
      
      # Extract chromosome names from the BSseq object
      chrs <- as.character(unique(seqnames(bs_subset)))
      dml_list <- list()
      
      for (chr in chrs) {
        message(paste("  Processing chromosome:", chr))
        
        # Subset BSseq object for just this chromosome
        # We index using row logicals because it is a SummarizedExperiment
        bs_subset_chr <- subset(bs_subset, seqnames(bs_subset) == chr)
        
        # If a chromosome ends up empty (e.g. some obscure scaffold), skip it
        if (nrow(bs_subset_chr) == 0) {
            next
        }
        
        # Run DMLtest for just this chromosome
        dml_chr <- tryCatch({
          DMLtest(
            bs_subset_chr,
            group1 = group1_samples,
            group2 = group2_samples,
            smoothing = TRUE
          )
        }, error = function(e) {
          message(paste("    Smoothing failed for", chr, "- retrying without smoothing. Error:", e$message))
          tryCatch({
            DMLtest(
              bs_subset_chr,
              group1 = group1_samples,
              group2 = group2_samples,
              smoothing = FALSE
            )
          }, error = function(e2) {
            message(paste("    Non-smoothed DMLtest also failed for", chr, "- skipping chromosome. Error:", e2$message))
            return(NULL)
          })
        })
        
        dml_list[[chr]] <- dml_chr
        rm(bs_subset_chr, dml_chr)
        gc(verbose = FALSE)
      }
      
      # Combine all the chromosome results into a single data frame
      message("Combining chromosome results...")
      dml_result <- do.call(rbind, dml_list)
      
      # Save results
      output_file <- paste0("dml_result_", comp_name, "_", sex_group, ".rds")
      output_path <- file.path(sub_dir, output_file)
      saveRDS(dml_result, file = output_path)
      message(paste("Saved:", output_file))
      
      # Clean memory
      rm(bs_subset, dml_result)
      gc(verbose = FALSE)
      
      return(output_path)
    }
    # -----------------------------------------------------------------

    ### Below, GROUP 1 & GROUP 2 = sex 1 (female) & sex 2 (male) ###
    ### We run this script for each tissue & separate by sex with Groups 1 & " ###
    ### For example for Group 1, change all "group1" to "female" and "group1_label" to "female" ###
    
    # --- GROUP 1 EXPLICIT COMPARISONS ---
    group1_samples <- sample_info[sample_info$sex == "group1_label", ]
    if(nrow(group1_samples) > 0) {
      saved_files$group1_high_vs_control <- run_specific_comparison(bsseq_obj, group1_samples, "group1_label", "high", "control")
      saved_files$group1_low_vs_control  <- run_specific_comparison(bsseq_obj, group1_samples, "group1_label", "low", "control")
      saved_files$group1_high_vs_low     <- run_specific_comparison(bsseq_obj, group1_samples, "group1_label", "high", "low")
    }
    
    # --- GROUP 2 EXPLICIT COMPARISONS (Ready for the future) ---
    group2_samples <- sample_info[sample_info$sex == "group2_label", ]
    if(nrow(group2_samples) > 0) {
      saved_files$group2_high_vs_control <- run_specific_comparison(bsseq_obj, group2_samples, "group2_label", "high", "control")
      saved_files$group2_low_vs_control  <- run_specific_comparison(bsseq_obj, group2_samples, "group2_label", "low", "control")
      saved_files$group2_high_vs_low     <- run_specific_comparison(bsseq_obj, group2_samples, "group2_label", "high", "low")
    }
  }
  
  message("\nAll explicit DSS analyses complete.")
  return(saved_files)
}


# --- Main Script Execution ---
### output directory should be "directory hdf5 object is stored/output directory name" use complete path ###

output_dir <- "<PATH_TO_OUTPUT_DIRECTORY>"
h5_bsseq_save_dir <- file.path(output_dir, "emseq_hdf5")


message("Importing BS object...")
tryCatch({
  se <- loadHDF5SummarizedExperiment(h5_bsseq_save_dir)
  bsseq_obj <- as(se, "BSseq")
  message("Raw BS object imported.")
}, error = function(e) {
  stop("Failed to load BSseq object: ", e$message)
})


# Sample metadata
# For metadata, use root sample IDs, such as "01_male_liver_control_rep1" #
# condition should be treatment group in the order the groups appear in metadata, i.e. "control", "low", "control", "high", ec. #
# for sex =c..."group1_label" = "male" or "female"; n_group = number males or females in your metadata #

sample_info <- data.frame(
    id = c(
    "<SAMPLE_ID_1>", "<SAMPLE_ID_2>", "<SAMPLE_ID_3>",
    "<SAMPLE_ID_4>", "<SAMPLE_ID_5>", "<SAMPLE_ID_6>",
    "<SAMPLE_ID_7>", "<SAMPLE_ID_8>", "<SAMPLE_ID_9>",
    "<SAMPLE_ID_10>", "<SAMPLE_ID_11>", "<SAMPLE_ID_12>",
    "<SAMPLE_ID_13>", "<SAMPLE_ID_14>", "<SAMPLE_ID_15>",
    "<SAMPLE_ID_16>", "<SAMPLE_ID_17>", "<SAMPLE_ID_18>",
    "<SAMPLE_ID_19>", "<SAMPLE_ID_20>", "<SAMPLE_ID_21>",
    "<SAMPLE_ID_22>", "<SAMPLE_ID_23>", "<SAMPLE_ID_24>",
    "<SAMPLE_ID_25>", "<SAMPLE_ID_26>", "<SAMPLE_ID_27>"
  ),
  condition = c(
    "<CONDITION_1>", "<CONDITION_2>", "<CONDITION_3>", "<CONDITION_4>",
    "<CONDITION_5>", "<CONDITION_6>", "<CONDITION_7>", "<CONDITION_8>",
    "<CONDITION_9>", "<CONDITION_10>", "<CONDITION_11>", "<CONDITION_12>",
    "<CONDITION_13>", "<CONDITION_14>", "<CONDITION_15>", "<CONDITION_16>",
    "<CONDITION_17>", "<CONDITION_18>", "<CONDITION_19>", "<CONDITION_20>",
    "<CONDITION_21>", "<CONDITION_22>", "<CONDITION_23>", "<CONDITION_24>",
    "<CONDITION_25>", "<CONDITION_26>", "<CONDITION_27>"
  ),
  sex = c(
    rep("group1_label", <N_GROUP1>),
    rep("group2_label", <N_GROUP2>)
  ),
  stringsAsFactors = FALSE
)


rownames(sample_info) <- sample_info$id
rownames(sample_info) <- sample_info$id
pData(bsseq_obj) <- sample_info
message("Sample metadata loaded and assigned.")
gc()


# --- Run Analysis ---
multi_factor_analysis <- FALSE 


dml_results_paths <- run_dss_analysis(bsseq_obj, sample_info, output_dir, multi_factor = multi_factor_analysis)


print("Saved file path(s):")
print(dml_results_paths)
