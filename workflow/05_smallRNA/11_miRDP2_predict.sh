#!/usr/bin/env bash
set -euo pipefail
source ./00_config.sh
mapfile -t SAMPLES < <(get_samples)

[[ -f "$MIRDP2_PIPELINE" ]] || { echo "[ERROR] Missing $MIRDP2_PIPELINE" >&2; exit 1; }

mkdir -p novel_miRNA_plant logs

for sample in "${SAMPLES[@]}"; do
    bash "$MIRDP2_PIPELINE" \
      --genome ref/genome.fa \
      --index ref/genome \
      --fasta \
      --input "novel_miRNA_plant/${sample}.reads.fa" \
      --output novel_miRNA_plant/ \
      > "logs/${sample}_novel_mir.log" 2>&1
done

bash "$MIRDP2_PIPELINE" \
  --genome ref/genome.fa \
  --index ref/genome \
  --fasta \
  --input novel_miRNA_plant/all.reads.fasta \
  --output novel_miRNA_plant/ \
  > logs/all_novel_mir.log 2>&1
