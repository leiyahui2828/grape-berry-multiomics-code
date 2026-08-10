#!/usr/bin/env bash
set -euo pipefail
source ./00_config.sh
mapfile -t SAMPLES < <(get_samples)

mkdir -p raw_bam logs
for sample in "${SAMPLES[@]}"; do
    run_meth_tools bismark \
      -p "$BISMARK_THREADS" \
      --genome ref/ \
      -o raw_bam/ \
      --basename "$sample" \
      -1 "clean_data/${sample}_1.fastq.gz" \
      -2 "clean_data/${sample}_2.fastq.gz" \
      > "logs/${sample}_bismark.log" 2>&1
done
