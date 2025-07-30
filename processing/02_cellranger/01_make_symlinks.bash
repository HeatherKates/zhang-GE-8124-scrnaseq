#!/bin/bash

csv="multi_config.csv"
outdir="symlinks"

mkdir -p "$outdir"

# Skip header
tail -n +2 "$csv" | while IFS=',' read -r fastq_id fastqs feature_types physical_library_id; do
  # Normalize feature type name for directory
  ft_dir=$(echo "$feature_types" | tr '[:upper:]' '[:lower:]' | tr ' ' '_')
  target_dir="${outdir}/${ft_dir}_fastqs"
  mkdir -p "$target_dir"

  # Remove quotes and split comma-separated fastq dirs
  clean_fastqs=$(echo "$fastqs" | sed 's/"//g')
  IFS=',' read -ra paths <<< "$clean_fastqs"

  for fq_path in "${paths[@]}"; do
    for fq_file in "$fq_path"/*fastq.gz; do
      ln -sf "$fq_file" "$target_dir/"
    done
  done

done
