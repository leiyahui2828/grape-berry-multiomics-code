#!/usr/bin/env bash
set -euo pipefail
source ./00_config.sh
mapfile -t SAMPLES < <(get_samples)

mkdir -p novel_miRNA_plant

for sample in "${SAMPLES[@]}"; do
    awk '$NF=="unknown" || $NF=="intron_sense" || $NF=="intron_antisense"' \
      "reads_mapping/${sample}.out.read_anno.txt" \
      > "novel_miRNA_plant/${sample}.reads.anno.txt"

    awk '{print ">"$1"\n"$2}' \
      "novel_miRNA_plant/${sample}.reads.anno.txt" \
      > "novel_miRNA_plant/${sample}.reads.fa"
done

cat novel_miRNA_plant/SR_PN*_rep*.reads.fa \
  > novel_miRNA_plant/all.reads.fasta
