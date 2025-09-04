#!/bin/bash

# Create refs directory if it doesn't exist
mkdir -p refs
cd refs

#########################
### 1. GEX REFERENCE  ###
#########################

echo "Downloading mm10 GEX reference..."

wget "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCm39-2024-A.tar.gz"
tar -xzf refdata-gex-GRCm39-2024-A.tar.gz
rm refdata-gex-GRCm39-2024-A.tar.gz

#########################
### 2. VDJ REFERENCE  ###
#########################

echo "Downloading mm10 VDJ reference..."

wget "https://cf.10xgenomics.com/supp/cell-vdj/refdata-cellranger-vdj-GRCm38-alts-ensembl-7.0.0.tar.gz"
tar -xzf refdata-cellranger-vdj-GRCm38-alts-ensembl-7.0.0.tar.gz
rm refdata-cellranger-vdj-GRCm38-alts-ensembl-7.0.0.tar.gz

echo "References built:"
echo " - GEX: refs/refdata-gex-GRCm39-2024-A"
echo " - VDJ: refs/refdata-cellranger-vdj-GRCm38-alts-ensembl-7.0.0"
cd ../
