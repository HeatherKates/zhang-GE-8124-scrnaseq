# =============================================================================
# PROTAC scRNA-seq Analysis: Cell Type Annotation and Differential Expression
# =============================================================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(readr)

# Load data
OUTPUT_DIR <- "../analysis_outputs"
obj <- readRDS(file.path(OUTPUT_DIR, "objects", "04_clustered.rds"))

message("=== Cell Type Annotation ===")

# Find markers for all clusters
cluster_markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25, 
                                  logfc.threshold = 0.25, test.use = "wilcox")

write_csv(cluster_markers, file.path(OUTPUT_DIR, "tables", "05_cluster_markers.csv"))

# Top markers per cluster for visualization
top_markers <- cluster_markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

# Known immune cell markers for mouse
immune_markers <- list(
  "T_cells" = c("Cd3e", "Cd3d", "Cd3g"),
  "CD4_T" = c("Cd4"),
  "CD8_T" = c("Cd8a", "Cd8b1"),  
  "Regulatory_T" = c("Foxp3", "Il2ra", "Ctla4", "Ikzf2"),
  "NK_cells" = c("Nkg7", "Gzma", "Klrb1c", "Ncr1"),
  "B_cells" = c("Cd19", "Ms4a1", "Cd79a", "Pax5"),
  "Plasma_cells" = c("Jchain", "Mzb1", "Xbp1", "Sdc1"),
  "Macrophages" = c("Adgre1", "Cd68", "Itgam", "Fcgr1"),
  "Dendritic_cells" = c("Itgax", "Flt3", "Zbtb46"),
  "Monocytes" = c("Ly6c2", "Ccr2", "Nr4a1"),
  "Neutrophils" = c("S100a8", "S100a9", "Ly6g")
)
library(patchwork)

# Plot immune markers
for (cell_type in names(immune_markers)) {
  markers <- immune_markers[[cell_type]]
  available_markers <- markers[markers %in% rownames(obj)]

  if (length(available_markers) > 0) {
    p <- FeaturePlot(obj, features = available_markers, ncol = 2) &
      theme(plot.title = element_text(size = 10))  # Make individual titles smaller
    
    # Add overall title
    p <- p + plot_annotation(title = cell_type, 
                             theme = theme(plot.title = element_text(size = 16, hjust = 0.5)))
    
    ggsave(file.path(OUTPUT_DIR, "plots", paste0("05_markers_", cell_type, ".pdf")),
           p, width = 12, height = 8)
  }
}

# Save objects and data for next step
saveRDS(obj, file.path(OUTPUT_DIR, "objects", "05_annotated.rds"))
saveRDS(cluster_markers, file.path(OUTPUT_DIR, "objects", "05_cluster_markers.rds"))
saveRDS(immune_markers, file.path(OUTPUT_DIR, "objects", "05_immune_markers.rds"))

message("=== Cell Annotation Complete ===")
message("Saved annotated object and marker data for differential expression analysis")
