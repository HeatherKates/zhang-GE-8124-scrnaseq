#!/bin/bash
#SBATCH --job-name=AL_cellranger
#SBATCH --output=logs/cellranger_AL.%j.out
#SBATCH --error=logs/cellranger_AL.%j.err
#SBATCH --qos=cancercenter-dept-b
#SBATCH --account=cancercenter-dept
#SBATCH --time=24:00:00
#SBATCH --mem=250G
#SBATCH --cpus-per-task=32

cd ${SLURM_SUBMIT_DIR}
module load cellranger/9.0.0
export TENX_IGNORE_DEPRECATED_OS=1

cellranger multi \
    --id=AL_hashing \
    --csv=multi_config_AL.csv \
    --localcores=${SLURM_CPUS_PER_TASK} \
    --localmem=$(( ${SLURM_MEM_PER_NODE} / 1024 ))
