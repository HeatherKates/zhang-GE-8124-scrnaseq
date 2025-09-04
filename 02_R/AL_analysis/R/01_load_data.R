# =============================================================================
# PROTAC scRNA-seq Analysis: Data Loading and Initial Setup (AL Group)
# =============================================================================

library(Seurat)
library(dplyr)
library(readr)
library(stringr)

# Configuration
DATASET <- "AL"  # "YT","AT","YL"
BASE_PATH <- paste0("/blue/zhangw/BCBSR_BIOINFORMATICS/hkates/TEMI/GE-8124_scrnaseq/processing/02_cellranger/", DATASET, "_run/", DATASET, "_hashing/outs")
samples <- paste0(DATASET, c("_Vehicle_1", "_Vehicle_2", "_Vehicle_3", "_PZ15227_1", "_PZ15227_2", "_PZ15227_3"))
age = "Aged" # "Young"
tissue = "Lymph node" # "Tumor"

#BASE_PATH <- "/blue/zhangw/BCBSR_BIOINFORMATICS/hkates/TEMI/GE-8124_scrnaseq/processing/02_cellranger/AL_run/AL_hashing/outs"
OUTPUT_DIR <- "../analysis_outputs"
WEB_DIR <- paste0("/orange/cancercenter-dept/web/public/BCB-SR/Zhangw/TEMI/GE-8124/",DATASET)

# Create directories
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "objects"), recursive = TRUE, showWarnings = FALSE)

message("=== PROTAC scRNA-seq Analysis: AL Group (Aged Lymph Node) ===")
message("Loading Cell Ranger Multi outputs...")

# Sample information
#samples <- c("AL_Vehicle_1", "AL_Vehicle_2", "AL_Vehicle_3", 
#             "AL_PZ15227_1", "AL_PZ15227_2", "AL_PZ15227_3")

sample_metadata <- data.frame(
  sample_id = samples,
  treatment = rep(c("Vehicle", "PZ15227"), each = 3),
  replicate = rep(paste0("Rep", 1:3), 2),
  age = age,
  tissue = tissue,
  condition = paste0(DATASET,"_", rep(c("Vehicle", "PZ15227"), each = 3)),
  stringsAsFactors = FALSE
)

load_sample_data <- function(sample_id, base_path) {
  message("Loading sample: ", sample_id)
  # Paths
  count_path <- file.path(base_path, "per_sample_outs", sample_id, "count",
                          "sample_filtered_feature_bc_matrix")
  vdj_b_path <- file.path(base_path, "per_sample_outs", sample_id, "vdj_b",
                          "filtered_contig_annotations.csv")
  vdj_t_path <- file.path(base_path, "per_sample_outs", sample_id, "vdj_t",
                          "filtered_contig_annotations.csv")
  metrics_path <- file.path(base_path, "per_sample_outs", sample_id, "metrics_summary.csv")
  
  # Load count data
  counts <- Read10X(count_path)
  
  # Get joint barcodes between Gene Expression and Antibody Capture
  gex_bcs <- colnames(counts[["Gene Expression"]])
  if ("Antibody Capture" %in% names(counts)) {
    hto_bcs <- colnames(counts[["Antibody Capture"]])
    joint_bcs <- intersect(gex_bcs, hto_bcs)
    message("  - Gene Expression: ", length(gex_bcs), " barcodes")
    message("  - Antibody Capture: ", length(hto_bcs), " barcodes")
    message("  - Joint barcodes: ", length(joint_bcs), " barcodes")
    
    # Subset both matrices to joint barcodes
    gex_counts <- counts[["Gene Expression"]][, joint_bcs]
    hto_counts <- counts[["Antibody Capture"]][, joint_bcs]
  } else {
    joint_bcs <- gex_bcs
    gex_counts <- counts[["Gene Expression"]]
    hto_counts <- NULL
  }
  
  # Create Seurat object with subsetted Gene Expression data
  obj <- CreateSeuratObject(
    counts = gex_counts,
    project = sample_id,
    min.cells = 3,
    min.features = 200
  )
  
  # Add HTO data if present - subset to cells that survived filtering
  if (!is.null(hto_counts)) {
    # Get the cells that survived Seurat filtering
    surviving_cells <- colnames(obj)
    # Subset HTO data to only include surviving cells
    hto_counts_filtered <- hto_counts[, surviving_cells]
    obj[["HTO"]] <- CreateAssayObject(counts = hto_counts_filtered)
    message("  - Added HTO data: ", nrow(hto_counts_filtered), " features, ", 
            ncol(hto_counts_filtered), " cells")
  }
  
  # Add sample metadata
  sample_info <- sample_metadata[sample_metadata$sample_id == sample_id, ]
  for (col in names(sample_info)) {
    if (col != "sample_id") {
      obj[[col]] <- sample_info[[col]]
    }
  }
  
  # Load VDJ data (unchanged)
  vdj_data <- list()
  if (file.exists(vdj_b_path)) {
    vdj_b <- read_csv(vdj_b_path, show_col_types = FALSE)
    vdj_data$bcr <- vdj_b
    message("  - BCR contigs: ", nrow(vdj_b))
  }
  if (file.exists(vdj_t_path)) {
    vdj_t <- read_csv(vdj_t_path, show_col_types = FALSE)
    vdj_data$tcr <- vdj_t
    message("  - TCR contigs: ", nrow(vdj_t))
  }
  
  # Load metrics
  if (file.exists(metrics_path)) {
    metrics <- read_csv(metrics_path, show_col_types = FALSE)
    vdj_data$metrics <- metrics
  }
  
  return(list(seurat = obj, vdj = vdj_data))
}
# Load all samples
message("Loading individual samples...")
sample_data <- list()

for (sample_id in samples) {
  sample_data[[sample_id]] <- load_sample_data(sample_id, BASE_PATH)
}

# Extract Seurat objects for merging
seurat_objects <- lapply(sample_data, function(x) x$seurat)

# Merge all samples
message("Merging samples...")
merged_obj <- merge(
  seurat_objects[[1]], 
  y = seurat_objects[2:length(seurat_objects)], 
  add.cell.ids = names(seurat_objects),
  project = "AL_PROTAC_Study"
)

message("Merged object: ", ncol(merged_obj), " cells, ", nrow(merged_obj), " genes")

# Add comprehensive metadata
for (col in colnames(sample_metadata)) {
  if (col != "sample_id") {
    # Extract sample ID from cell names (format: SAMPLE_BARCODE)
    cell_sample_ids <- str_extract(colnames(merged_obj), "^[^_]+_[^_]+_[^_]+")
    merged_obj[[col]] <- sample_metadata[[col]][match(cell_sample_ids, sample_metadata$sample_id)]
  }
}

# Compile VDJ data
message("Compiling VDJ data...")
all_vdj_data <- list(
  bcr = do.call(rbind, lapply(sample_data, function(x) x$vdj$bcr)),
  tcr = do.call(rbind, lapply(sample_data, function(x) x$vdj$tcr)),
  metrics = do.call(rbind, lapply(sample_data, function(x) x$vdj$metrics))
)

# Save initial objects
saveRDS(merged_obj, file.path(OUTPUT_DIR, "objects", "01_merged_raw.rds"))
saveRDS(all_vdj_data, file.path(OUTPUT_DIR, "objects", "01_vdj_data.rds"))
saveRDS(sample_metadata, file.path(OUTPUT_DIR, "objects", "01_sample_metadata.rds"))

# Create summary table
initial_summary <- data.frame(
  sample_id = samples,
  n_cells = sapply(seurat_objects, ncol),
  n_genes = sapply(seurat_objects, nrow),
  treatment = sample_metadata$treatment,
  replicate = sample_metadata$replicate
)

write_csv(initial_summary, file.path(OUTPUT_DIR, "tables", "01_sample_summary.csv"))

message("=== Data Loading Complete ===")
message("Files saved:")
message("- Merged object: ", file.path(OUTPUT_DIR, "objects", "01_merged_raw.rds"))
message("- VDJ data: ", file.path(OUTPUT_DIR, "objects", "01_vdj_data.rds"))
message("- Sample summary: ", file.path(OUTPUT_DIR, "tables", "01_sample_summary.csv"))

