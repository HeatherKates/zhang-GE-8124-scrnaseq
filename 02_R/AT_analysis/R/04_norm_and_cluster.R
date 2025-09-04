# =============================================================================
# PROTAC scRNA-seq Analysis: Normalization and Clustering  
# =============================================================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

# Load data
OUTPUT_DIR <- "../analysis_outputs"
obj <- readRDS(file.path(OUTPUT_DIR, "objects", "03_with_vdj.rds"))

message("=== Normalization and Clustering ===")

# Normalization and variable feature selection
message("Performing normalization...")
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)

# Identify highly variable features
top10 <- head(VariableFeatures(obj), 10)
plot1 <- VariableFeaturePlot(obj)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
ggsave(file.path(OUTPUT_DIR, "plots", "04_variable_features.pdf"), plot2, width = 12, height = 8)

# Scale data and PCA
message("Scaling data and running PCA...")
all.genes <- rownames(obj)
obj <- ScaleData(obj, features = all.genes)
obj <- RunPCA(obj, features = VariableFeatures(object = obj))

# PCA visualization
p1 <- DimPlot(obj, reduction = "pca", group.by = "treatment") + ggtitle("Treatment")
p2 <- DimPlot(obj, reduction = "pca", group.by = "sample_id") + ggtitle("Sample")
p3 <- DimPlot(obj, reduction = "pca", group.by = "has_bcr") + ggtitle("BCR Status")
p4 <- DimPlot(obj, reduction = "pca", group.by = "has_tcr") + ggtitle("TCR Status")

pca_plots <- (p1 + p2) / (p3 + p4)
ggsave(file.path(OUTPUT_DIR, "plots", "04_pca_overview.pdf"), pca_plots, width = 15, height = 12)

# Elbow plot
elbow_plot <- ElbowPlot(obj, ndims = 50)
ggsave(file.path(OUTPUT_DIR, "plots", "04_elbow_plot.pdf"), elbow_plot, width = 10, height = 6)

# Clustering and UMAP
message("Running clustering and UMAP...")
obj <- FindNeighbors(obj, dims = 1:25)
obj <- FindClusters(obj, resolution = 0.5)
obj <- RunUMAP(obj, dims = 1:25)

# UMAP visualization
p1 <- DimPlot(obj, reduction = "umap", label = TRUE, label.size = 6) + 
  ggtitle("Clusters") + NoLegend()
p2 <- DimPlot(obj, reduction = "umap", group.by = "treatment") + 
  ggtitle("Treatment")
p3 <- DimPlot(obj, reduction = "umap", group.by = "sample_id") + 
  ggtitle("Sample")
p4 <- DimPlot(obj, reduction = "umap", group.by = "has_bcr") + 
  ggtitle("BCR+ Cells")

umap_overview <- (p1 + p2) / (p3 + p4)
ggsave(file.path(OUTPUT_DIR, "plots", "04_umap_overview.pdf"), umap_overview, width = 15, height = 12)

# join layers, weird timing
obj <- JoinLayers(obj)

# VDJ-specific UMAP
if (sum(obj$has_tcr, na.rm = TRUE) > 0) {
  p5 <- DimPlot(obj, reduction = "umap", group.by = "has_tcr") + 
    ggtitle("TCR+ Cells")
  ggsave(file.path(OUTPUT_DIR, "plots", "04_umap_tcr.pdf"), p5, width = 8, height = 6)
}

# Sample distribution across clusters
sample_dist <- DimPlot(obj, reduction = "umap", split.by = "treatment") +
  ggtitle("Treatment Distribution")
ggsave(file.path(OUTPUT_DIR, "plots", "04_umap_split_treatment.pdf"), sample_dist, width = 12, height = 6)

# Save processed object
saveRDS(obj, file.path(OUTPUT_DIR, "objects", "04_clustered.rds"))

message("=== Clustering Complete ===")
message("Final object: ", ncol(obj), " cells, ", length(levels(obj$seurat_clusters)), " clusters")
