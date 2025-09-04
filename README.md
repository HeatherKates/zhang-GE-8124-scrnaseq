# PROTAC scRNA-seq Analysis: Effect on Immune Landscape in Young vs Aged Mice

This repository contains the analysis pipeline for Dr. Weizhou Zhang's study examining the effect of PROTAC treatment (PZ15227) on the immune landscape in young and old mice using single-cell RNA sequencing with cell hashing and VDJ sequencing.

## Study Design

**Experimental Groups:**
- **AL**: Aged mice, Lymph node samples
- **AT**: Aged mice, Tumor samples  
- **YL**: Young mice, Lymph node samples
- **YT**: Young mice, Tumor samples

**Treatment Conditions:**
- Vehicle control (3 biological replicates per group)
- PZ15227 treatment (3 biological replicates per group)

**Technologies:**
- 10X Genomics single-cell gene expression
- Cell hashing for sample demultiplexing
- VDJ-B (BCR) sequencing
- VDJ-T (TCR) sequencing

## Repository Structure

├── 01_preprocessing/          # Data preprocessing pipeline
│   ├── 01_fastqc/            # Quality control of raw reads
│   └── 02_cellranger/        # Cell Ranger Multi processing
│       ├── *_run/            # Group-specific runs (AL, AT, YL, YT)
│       ├── refs/             # Reference genomes and feature files
│       └── symlinks/         # Symlinked FASTQ files
├── 02_R/                     # R analysis pipeline
│   ├── AL_analysis/          # Aged lymph node analysis
│   ├── AT_analysis/          # Aged tumor analysis  
│   ├── YL_analysis/          # Young lymph node analysis
│   └── YT_analysis/          # Young tumor analysis
└── 03_reports/               # Summary reports and presentations

## Analysis Pipeline

### 1. Preprocessing (`01_preprocessing/`)

**FastQC (`01_fastqc/`):**
- Quality assessment of raw sequencing reads
- Run: `sbatch 01_fastqc.sbatch`

**Cell Ranger Multi (`02_cellranger/`):**
- Processes gene expression, cell hashing, and VDJ data simultaneously
- Each group processed separately to enable cell hashing demultiplexing
- Configuration files (`multi_config_*.csv`) specify:
  - Gene expression reference: `refdata-gex-GRCm39-2024-A`
  - VDJ reference: `refdata-cellranger-vdj-GRCm38-alts-ensembl-7.0.0`
  - Cell hashing antibodies: `feature_ref.csv`
  - Sample-to-hashtag assignments

**Example run:**
```bash
cd 01_preprocessing/02_cellranger/AL_run/
sbatch cellranger_AL.sbatch

### 2. R Analysis Pipeline (`02_R/`)

Each group follows an identical 9-step Seurat-based analysis:

1. **`01_load_data.R`** - Load Cell Ranger outputs, merge samples, compile VDJ data
2. **`02_qc_filtering.R`** - Quality control filtering and cell hashing demultiplexing  
3. **`03_vdj_integration.R`** - Integrate BCR/TCR data with gene expression
4. **`04_norm_and_cluster.R`** - Normalization, dimensionality reduction, clustering
5. **`05_cell_annotation.R`** - Cell type annotation
6. **`06_de.R`** - Differential expression analysis
7. **`07_prepare_report.R`** - Prepare data for reporting
8. **`08_report.Rmd`** - Generate HTML analysis report
9. **`09_loupe.R`** - Create Loupe files for 10X visualization

## Reproducing the Analysis

### Prerequisites
- Access to UF HiPerGator with Cell Ranger 9.0.0
- R 4.x with Seurat, dplyr, readr packages
- Raw FASTQ files (not included in repository)

### Step 1: Preprocessing
```bash
# Navigate to each group and run Cell Ranger
cd 01_preprocessing/02_cellranger/AL_run/
sbatch cellranger_AL.sbatch

# Repeat for AT, YL, YT groups

# For each group (AL, AT, YL, YT):
cd 02_R/AL_analysis/R/
Rscript 01_load_data.R
Rscript 02_qc_filtering.R
Rscript 03_vdj_integration.R
Rscript 04_norm_and_cluster.R
Rscript 05_cell_annotation.R
Rscript 06_de.R
Rscript 07_prepare_report.R

# Render analysis report
rmarkdown::render("08_report.Rmd")

# Create Loupe files
source("09_loupe.R")

## Results and Outputs

**Analysis outputs are automatically published to:**
- Main results: https://data.rc.ufl.edu/pub/cancercenter-dept/BCB-SR/Zhangw/TEMI/GE-8124/
- Group-specific results:
  - AL: https://data.rc.ufl.edu/pub/cancercenter-dept/BCB-SR/Zhangw/TEMI/GE-8124/AL/
  - AT: https://data.rc.ufl.edu/pub/cancercenter-dept/BCB-SR/Zhangw/TEMI/GE-8124/AT/
  - YL: https://data.rc.ufl.edu/pub/cancercenter-dept/BCB-SR/Zhangw/TEMI/GE-8124/YL/
  - YT: https://data.rc.ufl.edu/pub/cancercenter-dept/BCB-SR/Zhangw/TEMI/GE-8124/YT/

**Output files include:**
- Interactive HTML reports (`08_report.html`)
- 10X Loupe files for visualization
- Processed Seurat objects
- QC metrics and plots
- Differential expression results

## Notes

- Large reference files (`refdata-*`) and analysis outputs are excluded from the repository
- Cell Ranger outputs require significant computational resources and storage
- All paths in scripts assume UF HiPerGator file system structure
