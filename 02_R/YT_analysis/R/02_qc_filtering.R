# =============================================================================
# PROTAC scRNA-seq Analysis: Quality Control and Filtering
# =============================================================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

# Load data
OUTPUT_DIR <- "../analysis_outputs"
obj <- readRDS(file.path(OUTPUT_DIR, "objects", "01_merged_raw.rds"))
sample_metadata <- readRDS(file.path(OUTPUT_DIR, "objects", "01_sample_metadata.rds"))

message("=== Quality Control Analysis ===")
message("Starting with ", ncol(obj), " cells")

# Calculate QC metrics
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")
obj[["percent.ribo"]] <- PercentageFeatureSet(obj, pattern = "^Rpl|^Rps")
obj[["percent.hb"]] <- PercentageFeatureSet(obj, pattern = "^Hb[^(p)]")
# Add sample_id to existing object
cell_sample_ids <- str_extract(colnames(obj), "^[^_]+_[^_]+_[^_]+")
obj[["sample_id"]] <- cell_sample_ids

# Verify it worked
table(obj$sample_id)
# QC metrics summary
qc_summary <- obj@meta.data %>%
  group_by(sample_id, treatment) %>%
  summarise(
    n_cells = n(),
    median_genes = median(nFeature_RNA),
    median_umis = median(nCount_RNA),
    median_mt = median(percent.mt),
    median_ribo = median(percent.ribo),
    .groups = 'drop'
  )

write_csv(qc_summary, file.path(OUTPUT_DIR, "tables", "02_qc_summary_before.csv"))

# QC Visualizations
library(ggplot2)
p1 <- VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        group.by = "sample_id", ncol = 3, pt.size = 0) +
  ggtitle("QC Metrics by Sample")

p2 <- VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
              group.by = "treatment", ncol = 3, pt.size = 0.1) +
  ggtitle("QC Metrics by Treatment")

p3 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "percent.mt", 
                     group.by = "treatment") +
  ggtitle("UMI vs Mitochondrial %")

p4 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", 
                     group.by = "treatment") +
  ggtitle("UMI vs Gene Count")

# Save QC plots
ggsave(file.path(OUTPUT_DIR, "plots", "02_qc_by_sample.pdf"), p1, width = 15, height = 10)
ggsave(file.path(OUTPUT_DIR, "plots", "02_qc_by_treatment.pdf"), p2, width = 12, height = 8)
ggsave(file.path(OUTPUT_DIR, "plots", "02_qc_scatter.pdf"), (p3 | p4), width = 15, height = 6)

# Define QC thresholds (adjust based on your data)
min_genes <- 200
max_genes <- 5000
max_mt <- 20
min_umis <- 1000

message("Applying QC filters:")
message("- Min genes: ", min_genes)
message("- Max genes: ", max_genes) 
message("- Max MT%: ", max_mt)
message("- Min UMIs: ", min_umis)

# Apply filters
obj_filtered <- subset(obj, subset = 
                         nFeature_RNA > min_genes & 
                         nFeature_RNA < max_genes & 
                         percent.mt < max_mt &
                         nCount_RNA > min_umis &
                         percent.hb < 5
)

message("After filtering: ", ncol(obj_filtered), " cells retained")

# Filter summary by sample
filter_summary <- data.frame(
  sample_id = sample_metadata$sample_id,
  treatment = sample_metadata$treatment,
  cells_before = sapply(sample_metadata$sample_id, function(x) sum(obj$sample_id == x, na.rm = TRUE)),
  cells_after = sapply(sample_metadata$sample_id, function(x) sum(obj_filtered$sample_id == x, na.rm = TRUE))
) %>%
  mutate(
    cells_removed = cells_before - cells_after,
    percent_retained = round(100 * cells_after / cells_before, 2)
  )

write_csv(filter_summary, file.path(OUTPUT_DIR, "tables", "02_filter_summary.csv"))

# Post-filtering QC
p5 <- VlnPlot(obj_filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
              group.by = "treatment", ncol = 3, pt.size = 0.1) +
  ggtitle("QC Metrics After Filtering")

ggsave(file.path(OUTPUT_DIR, "plots", "02_qc_after_filtering.pdf"), p5, width = 12, height = 8)

# Save filtered object
saveRDS(obj_filtered, file.path(OUTPUT_DIR, "objects", "02_filtered.rds"))

message("=== QC and Filtering Complete ===")
message("Filtered object saved: ", file.path(OUTPUT_DIR, "objects", "02_filtered.rds"))
print(filter_summary)
