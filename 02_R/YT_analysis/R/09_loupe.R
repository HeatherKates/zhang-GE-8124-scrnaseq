# =============================================================================
# PROTAC scRNA-seq Analysis: Convert to Loupe Browser Format
# =============================================================================

library(Seurat)
library(loupeR)
library(dplyr)
select <- dplyr::select
# Set parameters
DATASET <- "YT"  # Change for each analysis
ANALYSIS_DIR <- paste0("/blue/zhangw/BCBSR_BIOINFORMATICS/hkates/TEMI/GE-8124_scrnaseq/02_R/",DATASET, "_analysis")
OUTPUT_DIR <- file.path(ANALYSIS_DIR, "analysis_outputs")
WEB_DIR <- paste0("/orange/cancercenter-dept/web/public/BCB-SR/Zhangw/TEMI/GE-8124/", DATASET)

# Load Seurat object
message("Loading Seurat object...")
obj <- readRDS(file.path(OUTPUT_DIR, "objects", "06_final_analysis.rds"))

# Prepare metadata for Loupe Browser
message("Preparing metadata...")
# Convert factor columns to character (Loupe Browser requirement)
metadata_cols <- colnames(obj@meta.data)
for(col in metadata_cols) {
  if(is.factor(obj@meta.data[[col]])) {
    obj@meta.data[[col]] <- as.character(obj@meta.data[[col]])
  }
}

# Create categorical metadata for filtering
obj$cluster_annotation <- paste0("Cluster_", obj$seurat_clusters)
obj$bcr_status <- ifelse(obj$has_bcr, "BCR_positive", "BCR_negative")
obj$tcr_status <- ifelse(obj$has_tcr, "TCR_positive", "TCR_negative")
obj$immune_type <- case_when(
  obj$has_bcr ~ "B_cells",
  obj$has_tcr ~ "T_cells",
  TRUE ~ "Other"
)

# Select key metadata columns for Loupe Browser
loupe_metadata <- obj@meta.data %>%
  select(
    sample_id,
    treatment,
    cluster_annotation,
    bcr_status,
    tcr_status,
    immune_type,
    nCount_RNA,
    nFeature_RNA,
    percent.mt
  )

# Create projections list (UMAP and PCA)
message("Preparing projections...")
projections <- list()

# Add UMAP if available
if("umap" %in% names(obj@reductions)) {
  projections[["UMAP"]] <- obj@reductions$umap@cell.embeddings
}

# Add PCA if available
if("pca" %in% names(obj@reductions)) {
  projections[["PCA"]] <- obj@reductions$pca@cell.embeddings[, 1:10]  # First 10 PCs
}

# Create clusters list
message("Preparing clusters...")
clusters <- list()

# Add Seurat clusters
clusters[["Seurat_Clusters"]] <- setNames(
  as.character(obj$seurat_clusters), 
  colnames(obj)
)

# Add treatment groups
clusters[["Treatment"]] <- setNames(
  as.character(obj$treatment), 
  colnames(obj)
)

# Add immune cell types
clusters[["Immune_Type"]] <- setNames(
  as.character(obj$immune_type), 
  colnames(obj)
)

# Create the Loupe file
message("Creating Loupe Browser file...")
output_file <- file.path(WEB_DIR, paste0(DATASET, "_analysis.cloupe"))

# Ensure output directory exists
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

# Convert to Loupe format
create_loupe_from_seurat(
  obj,
  output_name = output_file,
  executable_path = NULL     # Use system PATH
)

message("Loupe Browser file created: ", output_file)
message("File size: ", round(file.size(output_file) / 1024^2, 2), " MB")

# Create a summary for the researchers
summary_info <- data.frame(
  Dataset = DATASET,
  Total_Cells = ncol(obj),
  Total_Genes = nrow(obj),
  Clusters = length(unique(obj$seurat_clusters)),
  BCR_Positive = sum(obj$has_bcr),
  TCR_Positive = sum(obj$has_tcr),
  File_Path = output_file,
  File_Size_MB = round(file.size(output_file) / 1024^2, 2)
)

write.csv(summary_info, file.path(WEB_DIR, paste0(DATASET, "_loupe_summary.csv")), 
          row.names = FALSE)

message("=== Loupe Browser Conversion Complete ===")
message("Researchers can:")
message("1. Download the .cloupe file and open in Loupe Browser")
message("2. Filter cells by any metadata column (treatment, clusters, etc.)")
message("3. Explore gene expression interactively")
message("4. Create custom cell selections and export lists")