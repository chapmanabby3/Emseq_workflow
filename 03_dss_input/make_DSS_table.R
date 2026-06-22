#=========================================================
# Create BSseq object for DSS
#
# This script:
#   1. Reads Bismark CpG reports
#   2. Merges CpG coordinates across samples
#   3. Creates HDF5-backed methylation and coverage matrices
#   4. Constructs a BSseq object
#   5. Filters low-coverage CpGs
#   6. Saves the object for downstream DSS analysis
#=========================================================

#-----------------------------
# Working directory
#-----------------------------
setwd("PATH/TO/CpG_FILES")

#-----------------------------
# Required packages
#-----------------------------
library(DSS)
library(HDF5Array)
library(bsseq)
library(S4Vectors)
library(GenomicRanges)
library(data.table)
library(rhdf5)
library(SummarizedExperiment)

#-----------------------------
# Output directory
#-----------------------------
output_dir <- "PATH/TO/DSS_RESULTS"

if (!dir.exists(output_dir))
    dir.create(output_dir, recursive = TRUE)

h5_emseq_file <- file.path(output_dir, "emseq.h5")

#-----------------------------
# Sample metadata
#-----------------------------

file.list_c <- c(
    "sample1.CpG.txt.gz",
    "sample2.CpG.txt.gz",
    "sample3.CpG.txt.gz"
)

sample_info <- data.frame(

    id = c(
        "sample1",
        "sample2",
        "sample3"
    ),

    condition = c(
        "control",
        "low",
        "high"
    ),

    sex = c(
        "female",
        "female",
        "female"
    )
)

rownames(sample_info) <- sample_info$id

#-------------------------------------------------------
# Check input files
#-------------------------------------------------------

if (!all(file.exists(file.list_c))) {

    missing_files <- file.list_c[!file.exists(file.list_c)]

    stop(
        "Missing input files:\n",
        paste(missing_files, collapse = "\n")
    )

}

#=========================================================
# Create HDF5 methylation matrices
#=========================================================

message("Collecting CpG coordinates...")

num_samples <- length(file.list_c)

all_cpg_coords <- vector("list", num_samples)

for (i in seq_len(num_samples)) {

    temp <- fread(file.list_c[i], select = c(1,2), header = FALSE)

    all_cpg_coords[[i]] <- temp

}

unique_cpg_df <- unique(rbindlist(all_cpg_coords))

setorder(unique_cpg_df, V1, V2)

setnames(unique_cpg_df, c("chr","pos"))

unique_cpg_lookup <- copy(unique_cpg_df)

unique_cpg_lookup[, global_idx := .I]

num_sites <- nrow(unique_cpg_lookup)

if (file.exists(h5_emseq_file))
    file.remove(h5_emseq_file)

h5createFile(h5_emseq_file)

chunk_rows <- min(100000, num_sites)

h5createDataset(
    h5_emseq_file,
    "M",
    dims = c(num_sites, num_samples),
    storage.mode = "integer",
    chunk = c(chunk_rows, num_samples)
)

h5createDataset(
    h5_emseq_file,
    "N",
    dims = c(num_sites, num_samples),
    storage.mode = "integer",
    chunk = c(chunk_rows, num_samples)
)

h5_file <- H5Fopen(h5_emseq_file)

#=========================================================
# Populate HDF5 matrices
#=========================================================

for (i in seq_len(num_samples)) {

    message("Processing ", file.list_c[i])

    sample_data <- fread(
        file.list_c[i],
        select = c(1,2,4,5),
        header = FALSE
    )

    setnames(sample_data,
             c("chr","pos","meth","unmeth"))

    sample_data[, N := meth + unmeth]

    merged_data <- unique_cpg_lookup[
        sample_data,
        on = .(chr,pos)
    ]

    merged_data[is.na(meth), `:=`(meth = 0, N = 0)]

    setorder(merged_data, global_idx)

    h5write(
        merged_data$meth,
        h5_file,
        "M",
        index = list(NULL,i)
    )

    h5write(
        merged_data$N,
        h5_file,
        "N",
        index = list(NULL,i)
    )

}

H5Fclose(h5_file)

#=========================================================
# Create BSseq object
#=========================================================

M <- HDF5Array(h5_emseq_file, "M")
N <- HDF5Array(h5_emseq_file, "N")

gr_sites <- GRanges(

    seqnames = unique_cpg_df$chr,

    ranges = IRanges(
        start = unique_cpg_df$pos,
        width = 1
    )

)

BSobj <- BSseq(

    M = M,

    Cov = N,

    pData = DataFrame(sample_info),

    gr = gr_sites,

    sampleNames = sample_info$id

)

# Filter low coverage CpGs

BSobj <- BSobj[
    rowMeans(getCoverage(BSobj), na.rm = TRUE) >= 5,
]

#=========================================================
# Save object
#=========================================================

save_dir <- file.path(output_dir, "emseq_hdf5")

saveHDF5SummarizedExperiment(

    BSobj,

    dir = save_dir,

    replace = TRUE

)

message("Finished.")
