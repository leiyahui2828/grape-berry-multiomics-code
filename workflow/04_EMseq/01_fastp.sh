#!/usr/bin/env bash
set -euo pipefail
source ./00_config.sh
mapfile -t SAMPLES < <(get_samples)

mkdir -p clean_data logs
for sample in "${SAMPLES[@]}"; do
    fastp \
      --thread "$FASTP_THREADS" \
      -i "raw_data/${sample}_1.fastq.gz" \
      -I "raw_data/${sample}_2.fastq.gz" \
      -o "clean_data/${sample}_1.fastq.gz" \
      -O "clean_data/${sample}_2.fastq.gz" \
      -h "clean_data/${sample}.html" \
      -j "clean_data/${sample}.json" \
      > "logs/${sample}.fastp.log" 2>&1
done
