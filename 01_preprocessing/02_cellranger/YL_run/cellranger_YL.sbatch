#!/bin/bash
#SBATCH --job-name=YL_cellranger
#SBATCH --output=logs/cellranger_YL.%j.out
#SBATCH --error=logs/cellranger_YL.%j.err
#SBATCH --qos=cancercenter-dept-b
#SBATCH --account=cancercenter-dept
#SBATCH --mem=500G
#SBATCH --cpus-per-task=64
#SBATCH --time=24:00:00

cd ${SLURM_SUBMIT_DIR}
module load cellranger/9.0.0
export TENX_IGNORE_DEPRECATED_OS=1

cellranger multi \
    --id=YL_hashing \
    --csv=multi_config_YL.csv \
    --localcores=${SLURM_CPUS_PER_TASK} \
    --localmem=$(( ${SLURM_MEM_PER_NODE} / 1024 ))
