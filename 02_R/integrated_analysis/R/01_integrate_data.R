# integrate_harmony.R - Complete version with VDJ
library(Seurat)
library(harmony)
library(dplyr)
library(ggplot2)
library(readr)

# Define paths and groups (same as original)
base_path <- "/blue/zhangw/BCBSR_BIOINFORMATICS/hkates/TEMI/GE-8124_scrnaseq/02_R/"
groups <- c("AL_analysis", "AT_analysis", "YL_analysis", "YT_analysis")
output_dir <- file.path(base_path, "integrated_analysis/analysis_outputs")

# Create output directories
dir.create(file.path(output_dir, "objects"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "tables"), recursive = TRUE, showWarnings = FALSE)

message("=== Loading Seurat Objects with VDJ Data ===")
# Read objects from each group (SAME AS ORIGINAL)
seurat_objects <- list()
metadata_list <- list()
vdj_summaries <- list()

for(group in groups) {
  message("Loading ", group, "...")
  # Load main object with VDJ data
  seurat_objects[[group]] <- readRDS(file.path(base_path, group, "analysis_outputs/objects/03_with_vdj.rds"))
  # Add group identifier to metadata
  seurat_objects[[group]]$analysis_group <- group
  # Load additional metadata and summaries
  metadata_list[[group]] <- readRDS(file.path(base_path, group, "analysis_outputs/objects/01_sample_metadata.rds"))
  # Load VDJ summary
  vdj_summaries[[group]] <- read_csv(file.path(base_path, group, "analysis_outputs/tables/03_vdj_summary.csv"))
  
  message("  Cells: ", ncol(seurat_objects[[group]]))
  message("  BCR+ cells: ", sum(seurat_objects[[group]]$has_bcr, na.rm = TRUE))
  message("  TCR+ cells: ", sum(seurat_objects[[group]]$has_tcr, na.rm = TRUE))
}

message("=== Starting Harmony Integration ===")
# Merge objects (instead of complex integration)
integrated_obj <- merge(
  seurat_objects[[1]], 
  y = seurat_objects[-1],
  add.cell.ids = names(seurat_objects)
)

# Standard preprocessing
message("Preprocessing for Harmony...")
integrated_obj <- NormalizeData(integrated_obj, normalization.method = "LogNormalize", scale.factor = 10000)
integrated_obj <- FindVariableFeatures(integrated_obj, selection.method = "vst", nfeatures = 2000)
integrated_obj <- ScaleData(integrated_obj, verbose = FALSE)
integrated_obj <- RunPCA(integrated_obj, npcs = 50, verbose = FALSE)

# Find optimal number of PCs (elbow plot) - SAME AS ORIGINAL
pdf(file.path(output_dir, "plots/integration_elbow_plot.pdf"))
ElbowPlot(integrated_obj, ndims = 50)
dev.off()

# Run Harmony integration
message("Running Harmony...")
integrated_obj <- RunHarmony(
  integrated_obj,
  group.by.vars = "analysis_group",
  assay.use = "RNA",
  verbose = TRUE
)

# Run UMAP and clustering using Harmony reduction
integrated_obj <- RunUMAP(integrated_obj, reduction = "harmony", dims = 1:30)
integrated_obj <- FindNeighbors(integrated_obj, reduction = "harmony", dims = 1:30)
integrated_obj <- FindClusters(integrated_obj, resolution = 0.5)

message("=== Creating Summary Plots ===")
# Integration overview plots (SAME AS ORIGINAL)
p1 <- DimPlot(integrated_obj, reduction = "umap", group.by = "analysis_group",
              cols = c("#E31A1C", "#1F78B4", "#33A02C", "#FF7F00"))
ggsave(file.path(output_dir, "plots/integration_by_group.pdf"), p1, width = 10, height = 8)

p2 <- DimPlot(integrated_obj, reduction = "umap", group.by = "seurat_clusters", label = TRUE)
ggsave(file.path(output_dir, "plots/integration_clusters.pdf"), p2, width = 10, height = 8)

# VDJ visualization (SAME AS ORIGINAL)
p3 <- DimPlot(integrated_obj, reduction = "umap", group.by = "has_bcr",
              cols = c("grey90", "red")) + ggtitle("BCR+ Cells")
p4 <- DimPlot(integrated_obj, reduction = "umap", group.by = "has_tcr",
              cols = c("grey90", "blue")) + ggtitle("TCR+ Cells")
ggsave(file.path(output_dir, "plots/integration_vdj_overview.pdf"),
       (p3 | p4), width = 16, height = 6)

# Treatment visualization (SAME AS ORIGINAL)
if("treatment" %in% colnames(integrated_obj@meta.data)) {
  p5 <- DimPlot(integrated_obj, reduction = "umap", group.by = "treatment")
  ggsave(file.path(output_dir, "plots/integration_by_treatment.pdf"), p5, width = 10, height = 8)
}

message("=== Creating Summary Tables ===")
# Integration summary (SAME AS ORIGINAL)
integration_summary <- integrated_obj@meta.data %>%
  group_by(analysis_group) %>%
  summarise(
    total_cells = n(),
    bcr_positive = sum(has_bcr, na.rm = TRUE),
    tcr_positive = sum(has_tcr, na.rm = TRUE),
    bcr_percent = round(100 * bcr_positive / total_cells, 1),
    tcr_percent = round(100 * tcr_positive / total_cells, 1),
    .groups = 'drop'
  )
write_csv(integration_summary, file.path(output_dir, "tables/integration_summary.csv"))

# Combined VDJ summary (SAME AS ORIGINAL)
combined_vdj_summary <- bind_rows(vdj_summaries, .id = "analysis_group")
write_csv(combined_vdj_summary, file.path(output_dir, "tables/combined_vdj_summary.csv"))

# Save integrated object
saveRDS(integrated_obj, file.path(output_dir, "objects/integrated_with_vdj.rds"))

message("=== Integration Complete ===")
message("Total integrated cells: ", ncol(integrated_obj))
message("BCR+ cells: ", sum(integrated_obj$has_bcr, na.rm = TRUE))
message("TCR+ cells: ", sum(integrated_obj$has_tcr, na.rm = TRUE))
print(integration_summary)