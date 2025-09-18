# make-cloupe-integrated-dataset.R
# =============================================================================
# PROTAC scRNA-seq Analysis: Convert Integrated Dataset to Loupe Browser Format
# =============================================================================

library(Seurat)
library(loupeR)
library(dplyr)

# Set parameters
DATASET <- "INTEGRATED"
BASE_PATH <- "/blue/zhangw/BCBSR_BIOINFORMATICS/hkates/TEMI/GE-8124_scrnaseq/02_R/"
ANALYSIS_DIR <- file.path(BASE_PATH, "integrated_analysis")
OUTPUT_DIR <- file.path(ANALYSIS_DIR, "analysis_outputs")
WEB_DIR = "/orange/cancercenter-dept/web/secure/zhangw/GE-8124/integrated"

# Load integrated Seurat object
message("Loading integrated Seurat object...")
obj <- readRDS(file.path(OUTPUT_DIR, "objects", "integrated_with_vdj.rds"))

# Verify what assays and reductions are available
message("Available assays: ", paste(names(obj@assays), collapse = ", "))
message("Available reductions: ", paste(names(obj@reductions), collapse = ", "))

# Set default assay to RNA
DefaultAssay(obj) <- "RNA"

# FIX 1: Join layers in Seurat v5 (this fixes the warning and potential issues)
message("Joining layers in RNA assay...")
obj[["RNA"]] <- JoinLayers(obj[["RNA"]])

# Prepare metadata for Loupe Browser
message("Preparing metadata...")
# Convert factor columns to character (Loupe Browser requirement)
metadata_cols <- colnames(obj@meta.data)
for(col in metadata_cols) {
  if(is.factor(obj@meta.data[[col]])) {
    obj@meta.data[[col]] <- as.character(obj@meta.data[[col]])
  }
}

# Create enhanced categorical metadata for filtering
obj$cluster_annotation <- paste0("Cluster_", obj$seurat_clusters)
obj$bcr_status <- ifelse(obj$has_bcr %in% TRUE, "BCR_positive", "BCR_negative")
obj$tcr_status <- ifelse(obj$has_tcr %in% TRUE, "TCR_positive", "TCR_negative")
obj$immune_type <- case_when(
  obj$has_bcr %in% TRUE ~ "B_cells",
  obj$has_tcr %in% TRUE ~ "T_cells",
  TRUE ~ "Other"
)

# Create analysis group labels for better visualization
obj$analysis_group_clean <- case_when(
  obj$analysis_group == "AL_analysis" ~ "Aged_Lymph_node",
  obj$analysis_group == "AT_analysis" ~ "Aged_Tumor",
  obj$analysis_group == "YL_analysis" ~ "Young_Lymph_node",
  obj$analysis_group == "YT_analysis" ~ "Young_Tumor",
  TRUE ~ obj$analysis_group
)
obj$group_treatment <- paste(obj$analysis_group_clean, obj$treatment, sep = "_")

# FIX 2: Use create_loupe_from_seurat() without manual cluster/projection setup
# This should handle the formatting automatically
message("Creating Loupe Browser file...")
output_file <- file.path(WEB_DIR, paste0(DATASET, "_analysis.cloupe"))

# Ensure output directory exists
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

# Convert to Loupe format using simplified approach
create_loupe_from_seurat(
  obj,
  output_name = output_file,
  force = TRUE  # Overwrite if exists
)

message("Loupe Browser file created: ", output_file)
message("File size: ", round(file.size(output_file) / 1024^2, 2), " MB")

# Create a detailed summary for the researchers
summary_info <- data.frame(
  Dataset = DATASET,
  Total_Cells = ncol(obj),
  Total_Genes = nrow(obj),
  Total_Clusters = length(unique(obj$seurat_clusters)),
  Analysis_Groups = length(unique(obj$analysis_group)),
  BCR_Positive = sum(obj$has_bcr %in% TRUE, na.rm = TRUE),
  TCR_Positive = sum(obj$has_tcr %in% TRUE, na.rm = TRUE),
  File_Path = output_file,
  File_Size_MB = round(file.size(output_file) / 1024^2, 2)
)

# Create group-specific summary
group_summary <- obj@meta.data %>%
  group_by(analysis_group_clean) %>%
  summarise(
    cells = n(),
    bcr_pos = sum(has_bcr %in% TRUE, na.rm = TRUE),
    tcr_pos = sum(has_tcr %in% TRUE, na.rm = TRUE),
    .groups = 'drop'
  )

# Save summaries
write.csv(summary_info, file.path(WEB_DIR, paste0(DATASET, "_loupe_summary.csv")),
          row.names = FALSE)
write.csv(group_summary, file.path(WEB_DIR, paste0(DATASET, "_group_breakdown.csv")),
          row.names = FALSE)

message("=== Integrated Loupe Browser Conversion Complete ===")
message("Total cells: ", ncol(obj))
message("BCR+ cells: ", sum(obj$has_bcr %in% TRUE, na.rm = TRUE))
message("TCR+ cells: ", sum(obj$has_tcr %in% TRUE, na.rm = TRUE))
message("\nResearchers can:")
message("1. Download the .cloupe file and open in Loupe Browser")
message("2. Filter by analysis group (Aged/Young Lymph_node/Tumor)")
message("3. Compare treatments across integrated dataset")
message("4. Explore BCR/TCR patterns across all samples")
message("5. Analyze integrated clusters spanning multiple conditions")

print(group_summary)
