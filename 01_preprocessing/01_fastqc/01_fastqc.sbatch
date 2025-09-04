#!/bin/bash
#SBATCH --job-name=fastqc
#SBATCH --output=01_fastqc/logs/fastqc_%A_%a.out
#SBATCH --error=01_fastqc/logs/fastqc_%A_%a.err
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=2:00:00
#SBATCH --qos=timgarrett-b
#SBATCH --account=timgarrett
#SBATCH --array=1-44

module load fastqc

LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" 01_fastqc/sample_fastqs.tsv)
SAMPLE=$(echo "$LINE" | cut -f1)
FQ1=$(echo "$LINE" | cut -f2)
FQ2=$(echo "$LINE" | cut -f3)
OUTDIR=01_fastqc/reports

echo "Running FastQC for $SAMPLE:"
echo "  R1: $FQ1"
echo "  R2: $FQ2"

fastqc -t 2 -o "$OUTDIR" "$FQ1" "$FQ2"
