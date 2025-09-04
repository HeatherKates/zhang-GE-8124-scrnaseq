#!/bin/bash
#SBATCH --job-name=YT_cellranger
#SBATCH --output=logs/cellranger_YT.%j.out
#SBATCH --error=logs/cellranger_YT.%j.err
#SBATCH --qos=timgarrett-b
#SBATCH --account=timgarrett
#SBATCH --time=8:00:00
#SBATCH --mem=250G
#SBATCH --cpus-per-task=32

cd ${SLURM_SUBMIT_DIR}
module load cellranger/9.0.0
export TENX_IGNORE_DEPRECATED_OS=1

cellranger multi \
    --id=YT_hashing \
    --csv=multi_config_YT.csv \
    --localcores=${SLURM_CPUS_PER_TASK} \
    --localmem=$(( ${SLURM_MEM_PER_NODE} / 1024 ))
