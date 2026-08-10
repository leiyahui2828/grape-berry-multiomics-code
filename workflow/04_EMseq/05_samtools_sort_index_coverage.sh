#!/usr/bin/env bash
set -euo pipefail
source ./00_config.sh
mapfile -t SAMPLES < <(get_samples)

mkdir -p raw_bam logs
for sample in "${SAMPLES[@]}"; do
    samtools sort -@ "$SAMTOOLS_THREADS" \
      -o "raw_bam/${sample}.sorted.bam" \
      "raw_bam/${sample}_pe.bam" \
      > "logs/${sample}_bam_sorted.log" 2>&1

    samtools index "raw_bam/${sample}.sorted.bam" \
      > "logs/${sample}_bam_index.log" 2>&1

    samtools coverage "raw_bam/${sample}.sorted.bam" \
      > "raw_bam/${sample}.sorted.coverage.txt"
done
