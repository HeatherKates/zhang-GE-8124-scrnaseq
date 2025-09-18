# First, create the file copying script (run before knitting):
# 07_prepare_report.R

library(here)

# Set dataset
DATASET <- "YL"  # Change for each analysis

# Define paths
ANALYSIS_DIR <- paste0("/blue/zhangw/BCBSR_BIOINFORMATICS/hkates/TEMI/GE-8124_scrnaseq/02_R/",DATASET, "_analysis")
OUTPUT_DIR <- file.path(ANALYSIS_DIR, "analysis_outputs")
WEB_DIR <- paste0("/orange/cancercenter-dept/web/public/BCB-SR/Zhangw/TEMI/GE-8124/", DATASET)
BASE_PATH <- paste0("/blue/zhangw/BCBSR_BIOINFORMATICS/hkates/TEMI/GE-8124_scrnaseq/01_preprocessing/02_cellranger/", DATASET, "_run/", DATASET, "_hashing/outs")

# Create web directory
dir.create(WEB_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(WEB_DIR, "plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(WEB_DIR, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(WEB_DIR, "web_summaries"), recursive = TRUE, showWarnings = FALSE)

# Copy analysis outputs
file.copy(file.path(OUTPUT_DIR, "plots"), WEB_DIR, recursive = TRUE, overwrite = TRUE)
file.copy(file.path(OUTPUT_DIR, "tables"), WEB_DIR, recursive = TRUE, overwrite = TRUE)

# Copy Cell Ranger web summaries
sample_dirs <- list.dirs(file.path(BASE_PATH, "per_sample_outs"), full.names = TRUE, recursive = FALSE)
for (sample_dir in sample_dirs) {
  sample_name <- basename(sample_dir)
  web_summary_src <- file.path(sample_dir, "web_summary.html")
  if (file.exists(web_summary_src)) {
    file.copy(web_summary_src, file.path(WEB_DIR, "web_summaries", paste0(sample_name, "_web_summary.html")))
  }
}

message("Files copied to: ", WEB_DIR)
