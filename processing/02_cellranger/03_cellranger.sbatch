#!/bin/bash
#SBATCH --job-name=cellranger_multi
#SBATCH --output=logs/cellrangermulti.out
#SBATCH --error=logs/cellrangermulti.err
#SBATCH --qos=timgarrett-b
#SBATCH --account=timgarrett
#SBATCH --time=24:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=16

module load cellranger/9.0.0

cellranger multi \
  --id=mouse_hashing_run \
  --csv=multi_config.csv
