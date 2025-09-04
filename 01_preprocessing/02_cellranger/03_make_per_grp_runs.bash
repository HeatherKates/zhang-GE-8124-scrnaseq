# Generate all configs and scripts
for GROUP in AL AT YL YT; do
    mkdir -p ${GROUP}_run/logs
    
    # Create config file
    cat > ${GROUP}_run/multi_config_${GROUP}.csv << EOF
[gene-expression]
reference,/blue/zhangw/BCBSR_BIOINFORMATICS/hkates/TEMI/GE-8124_scrnaseq/processing/02_cellranger/refs/refdata-gex-GRCm39-2024-A
create-bam,true

[vdj]
reference,/blue/zhangw/BCBSR_BIOINFORMATICS/hkates/TEMI/GE-8124_scrnaseq/processing/02_cellranger/refs/refdata-cellranger-vdj-GRCm38-alts-ensembl-7.0.0

[feature]
reference,/blue/zhangw/BCBSR_BIOINFORMATICS/hkates/TEMI/GE-8124_scrnaseq/processing/02_cellranger/refs/feature_ref.csv

[libraries]
fastq_id,fastqs,feature_types,physical_library_id
BCR-${GROUP},/blue/zhangw/BCBSR_BIOINFORMATICS/hkates/TEMI/GE-8124_scrnaseq/processing/02_cellranger/symlinks/vdjb_fastqs,VDJ-B,BCR
TCR-${GROUP},/blue/zhangw/BCBSR_BIOINFORMATICS/hkates/TEMI/GE-8124_scrnaseq/processing/02_cellranger/symlinks/vdjt_fastqs,VDJ-T,TCR
CSP-${GROUP},/blue/zhangw/BCBSR_BIOINFORMATICS/hkates/TEMI/GE-8124_scrnaseq/processing/02_cellranger/symlinks/ab_fastqs,Antibody Capture,CSP
GE-${GROUP},/blue/zhangw/BCBSR_BIOINFORMATICS/hkates/TEMI/GE-8124_scrnaseq/processing/02_cellranger/symlinks/gex_fastqs,Gene Expression,GE

[samples]
sample_id,cmo_ids
${GROUP}_Vehicle_1,HTO1
${GROUP}_Vehicle_2,HTO2
${GROUP}_Vehicle_3,HTO3
${GROUP}_PZ15227_1,HTO4
${GROUP}_PZ15227_2,HTO5
${GROUP}_PZ15227_3,HTO6
EOF

    # Create SLURM script
    cat > ${GROUP}_run/cellranger_${GROUP}.sbatch << EOF
#!/bin/bash
#SBATCH --job-name=cellranger_${GROUP}
#SBATCH --output=logs/cellranger_${GROUP}.%j.out
#SBATCH --error=logs/cellranger_${GROUP}.%j.err
#SBATCH --qos=timgarrett-b
#SBATCH --account=timgarrett
#SBATCH --time=8:00:00
#SBATCH --mem=250G
#SBATCH --cpus-per-task=48

cd \${SLURM_SUBMIT_DIR}
module load cellranger/9.0.0
export TENX_IGNORE_DEPRECATED_OS=1

cellranger multi \\
    --id=${GROUP}_hashing \\
    --csv=multi_config_${GROUP}.csv \\
    --localcores=\${SLURM_CPUS_PER_TASK} \\
    --localmem=\$(( \${SLURM_MEM_PER_NODE} / 1024 ))
EOF

done
