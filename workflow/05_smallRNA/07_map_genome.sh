#!/usr/bin/env bash
set -euo pipefail
source ./00_config.sh
mapfile -t SAMPLES < <(get_samples)

[[ -f scripts/stat_bwt.pl ]] || { echo "[ERROR] Missing scripts/stat_bwt.pl" >&2; exit 1; }

mkdir -p reads_mapping logs

for sample in "${SAMPLES[@]}"; do
    bowtie \
      -f \
      -v 0 \
      -k 1 \
      --al "reads_mapping/${sample}.mapped.fa" \
      ref/genome \
      "clean_data/uniq_data/${sample}.fa" \
      > "reads_mapping/${sample}.genome.bwt" \
      2> "logs/${sample}.bowtie_genome.log"

    perl scripts/stat_bwt.pl \
      "clean_data/uniq_data/${sample}.fa" \
      "reads_mapping/${sample}.genome.bwt" \
      > "reads_mapping/${sample}.genome.stat.tsv"
done
