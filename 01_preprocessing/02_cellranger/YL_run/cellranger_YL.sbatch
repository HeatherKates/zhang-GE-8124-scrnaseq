#!/bin/bash
#SBATCH --job-name=YLpri_cellranger
#SBATCH --output=logs/cellranger_YL.%j.out
#SBATCH --error=logs/cellranger_YL.%j.err
#SBATCH --qos=cancercenter-dept-b
#SBATCH --account=cancercenter-dept
#SBATCH --mem=350G
#SBATCH --cpus-per-task=24
#SBATCH --time=12:00:00

cd ${SLURM_SUBMIT_DIR}
module load cellranger/9.0.0
export TENX_IGNORE_DEPRECATED_OS=1

cellranger multi \
    --id=YL_hashingpri \
    --csv=multi_config_YL.csv \
    --localcores=24 \
    --localmem=14
