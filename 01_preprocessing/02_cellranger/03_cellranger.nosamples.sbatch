#!/bin/bash
#SBATCH --job-name=cellranger_multi
#SBATCH --output=logs/cellrangermulti.%j.out
#SBATCH --error=logs/cellrangermulti%j.err
#SBATCH --qos=timgarrett-b
#SBATCH --account=timgarrett
#SBATCH --time=96:00:00
#SBATCH --mem=900G
#SBATCH --cpus-per-task=32

module load cellranger/9.0.0

cellranger multi \
  --id=mouse_hashing_run_nosamples \
  --csv=multi_config.nosamples.csv
