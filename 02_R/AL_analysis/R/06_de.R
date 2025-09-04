# =============================================================================
# PROTAC scRNA-seq Analysis: Differential Expression Analysis
# =============================================================================
library(Seurat)
library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)

# Load data
OUTPUT_DIR <- "../analysis_outputs"
obj <- readRDS(file.path(OUTPUT_DIR, "objects", "05_annotated.rds"))
cluster_markers <- readRDS(file.path(OUTPUT_DIR, "objects", "05_cluster_markers.rds"))

message("=== Differential Expression Analysis ===")
message("Starting with ", ncol(obj), " cells across ", length(levels(obj$seurat_clusters)), " clusters")

# Treatment differential expression
message("Performing treatment differential expression...")
Idents(obj) <- "treatment"
treatment_de <- FindMarkers(obj, ident.1 = "PZ15227", ident.2 = "Vehicle",
                            test.use = "MAST",
                            latent.vars = "sample_id",
                            logfc.threshold = 0.25)
treatment_de$gene <- rownames(treatment_de)
write_csv(treatment_de, file.path(OUTPUT_DIR, "tables", "06_treatment_de_MAST.csv"))

# Per-cluster treatment effects using MAST
cluster_treatment_de <- list()

for (cluster in unique(obj$seurat_clusters)) {
  cells_in_cluster <- subset(obj, seurat_clusters == cluster)
  
  if (length(unique(cells_in_cluster$treatment)) == 2 && 
      min(table(cells_in_cluster$treatment)) >= 10) {
    
    Idents(cells_in_cluster) <- "treatment"
    
    # MAST with sample as latent variable
    de_result <- FindMarkers(cells_in_cluster, 
                             ident.1 = "PZ15227", ident.2 = "Vehicle",
                             test.use = "MAST",
                             latent.vars = "sample_id",  # This accounts for sample effects!
                             logfc.threshold = 0.25)
    
    de_result$gene <- rownames(de_result)
    de_result$cluster <- cluster
    
    cluster_treatment_de[[paste0("cluster_", cluster)]] <- de_result
    write_csv(de_result, file.path(OUTPUT_DIR, "tables", paste0("06_treatment_de_cluster_", cluster, "_MAST.csv")))
  }
}

# Summary statistics
summary_stats <- obj@meta.data %>%
  group_by(treatment, seurat_clusters) %>%
  summarise(
    n_cells = n(),
    bcr_positive = sum(has_bcr, na.rm = TRUE),
    tcr_positive = sum(has_tcr, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  pivot_wider(names_from = treatment, values_from = c(n_cells, bcr_positive, tcr_positive),
              names_sep = "_", values_fill = 0)

write_csv(summary_stats, file.path(OUTPUT_DIR, "tables", "06_cluster_summary.csv"))

# Save final object
saveRDS(obj, file.path(OUTPUT_DIR, "objects", "06_final_analysis.rds"))

message("=== Analysis Complete ===")
message("Treatment DE genes (adj. p. val. < 0.05): ", sum(treatment_de$p_val_adj < 0.05))
message("Files saved in: ", OUTPUT_DIR)
