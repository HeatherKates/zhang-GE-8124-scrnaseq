# =============================================================================
# PROTAC scRNA-seq Analysis: VDJ Integration
# =============================================================================

library(Seurat)
library(dplyr)
library(readr)
library(tidyr)

# Load data
OUTPUT_DIR <- "../analysis_outputs"
obj <- readRDS(file.path(OUTPUT_DIR, "objects", "02_filtered.rds"))
vdj_data <- readRDS(file.path(OUTPUT_DIR, "objects", "01_vdj_data.rds"))

message("=== VDJ Data Integration ===")

# Function to process VDJ data for integration
process_vdj_data <- function(vdj_df, chain_type = "BCR") {
  if (is.null(vdj_df) || nrow(vdj_df) == 0) {
    message("No ", chain_type, " data found")
    return(NULL)
  }
  
  message("Processing ", chain_type, " data: ", nrow(vdj_df), " contigs")
  
  # Basic filtering - only productive, high confidence contigs
  vdj_filtered <- vdj_df %>%
    filter(
      high_confidence == TRUE,
      productive == TRUE,
      is_cell == TRUE
    )
  
  message("After filtering: ", nrow(vdj_filtered), " contigs")
  
  if (chain_type == "BCR") {
    # Process BCR data (IGH, IGK, IGL)
    vdj_summary <- vdj_filtered %>%
      group_by(barcode) %>%
      summarise(
        bcr_chains = paste(unique(chain), collapse = ";"),
        bcr_heavy = any(chain == "IGH"),
        bcr_light_kappa = any(chain == "IGK"),
        bcr_light_lambda = any(chain == "IGL"),
        bcr_v_genes = paste(unique(v_gene[!is.na(v_gene)]), collapse = ";"),
        bcr_j_genes = paste(unique(j_gene[!is.na(j_gene)]), collapse = ";"),
        bcr_cdr3s = paste(unique(cdr3[!is.na(cdr3)]), collapse = ";"),
        bcr_productive_contigs = n(),
        .groups = 'drop'
      ) %>%
      mutate(
        bcr_complete = bcr_heavy & (bcr_light_kappa | bcr_light_lambda),
        bcr_chain_type = case_when(
          bcr_light_kappa & !bcr_light_lambda ~ "IGH_IGK",
          bcr_light_lambda & !bcr_light_kappa ~ "IGH_IGL",
          bcr_light_kappa & bcr_light_lambda ~ "IGH_IGK_IGL",
          TRUE ~ "IGH_only"
        )
      )
    
  } else {
    # Process TCR data (TRA, TRB)
    vdj_summary <- vdj_filtered %>%
      group_by(barcode) %>%
      summarise(
        tcr_chains = paste(unique(chain), collapse = ";"),
        tcr_alpha = any(chain == "TRA"),
        tcr_beta = any(chain == "TRB"),
        tcr_v_genes = paste(unique(v_gene[!is.na(v_gene)]), collapse = ";"),
        tcr_j_genes = paste(unique(j_gene[!is.na(j_gene)]), collapse = ";"),
        tcr_cdr3s = paste(unique(cdr3[!is.na(cdr3)]), collapse = ";"),
        tcr_productive_contigs = n(),
        .groups = 'drop'
      ) %>%
      mutate(
        tcr_complete = tcr_alpha & tcr_beta,
        tcr_chain_type = case_when(
          tcr_alpha & tcr_beta ~ "TRA_TRB",
          tcr_alpha & !tcr_beta ~ "TRA_only", 
          !tcr_alpha & tcr_beta ~ "TRB_only",
          TRUE ~ "unknown"
        )
      )
  }
  
  return(vdj_summary)
}

# Process BCR and TCR data
bcr_summary <- process_vdj_data(vdj_data$bcr, "BCR")
tcr_summary <- process_vdj_data(vdj_data$tcr, "TCR")

# Add VDJ information to Seurat object
# Initialize all cells as negative
obj$has_bcr <- FALSE
obj$has_tcr <- FALSE
obj$bcr_complete <- FALSE
obj$tcr_complete <- FALSE

# Add BCR data
if (!is.null(bcr_summary)) {
  # Match barcodes (need to account for sample prefixes)
  obj_barcodes <- colnames(obj)
  # Extract just the barcode part (after the sample prefix)
  obj_barcode_only <- sub(".*_", "", obj_barcodes)
  
  bcr_match <- match(obj_barcode_only, bcr_summary$barcode)
  
  for (col in names(bcr_summary)) {
    if (col != "barcode") {
      obj[[paste0("bcr_", col)]] <- bcr_summary[[col]][bcr_match]
    }
  }
  
  obj$has_bcr <- !is.na(bcr_match)
  obj$bcr_complete <- obj$has_bcr & !is.na(obj$bcr_bcr_complete) & obj$bcr_bcr_complete
  
  message("BCR+ cells: ", sum(obj$has_bcr, na.rm = TRUE))
}

# Add TCR data  
if (!is.null(tcr_summary)) {
  tcr_match <- match(obj_barcode_only, tcr_summary$barcode)
  
  for (col in names(tcr_summary)) {
    if (col != "barcode") {
      obj[[paste0("tcr_", col)]] <- tcr_summary[[col]][tcr_match]
    }
  }
  
  obj$has_tcr <- !is.na(tcr_match)
  obj$tcr_complete <- obj$has_tcr & !is.na(obj$tcr_tcr_complete) & obj$tcr_tcr_complete
  
  message("TCR+ cells: ", sum(obj$has_tcr, na.rm = TRUE))
}

# VDJ Summary statistics
vdj_summary_stats <- obj@meta.data %>%
  group_by(treatment, sample_id) %>%
  summarise(
    total_cells = n(),
    bcr_positive = sum(has_bcr, na.rm = TRUE),
    tcr_positive = sum(has_tcr, na.rm = TRUE),
    bcr_complete = sum(bcr_complete, na.rm = TRUE),
    tcr_complete = sum(tcr_complete, na.rm = TRUE),
    bcr_percent = round(100 * bcr_positive / total_cells, 1),
    tcr_percent = round(100 * tcr_positive / total_cells, 1),
    .groups = 'drop'
  )

write_csv(vdj_summary_stats, file.path(OUTPUT_DIR, "tables", "03_vdj_summary.csv"))

# Save object with VDJ data
saveRDS(obj, file.path(OUTPUT_DIR, "objects", "03_with_vdj.rds"))

message("=== VDJ Integration Complete ===")
message("Total cells: ", ncol(obj))
message("BCR+ cells: ", sum(obj$has_bcr, na.rm = TRUE))
message("TCR+ cells: ", sum(obj$has_tcr, na.rm = TRUE))
print(vdj_summary_stats)