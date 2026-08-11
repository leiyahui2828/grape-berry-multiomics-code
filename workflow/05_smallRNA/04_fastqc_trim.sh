#!/usr/bin/env bash
set -euo pipefail
source ./00_config.sh
mapfile -t SAMPLES < <(get_samples)

mkdir -p raw_data clean_data clean_data/uniq_data logs quality_control

for sample in "${SAMPLES[@]}"; do
    fq="raw_data/${sample}_R1.fastq.gz"
    [[ -f "$fq" ]] || { echo "[ERROR] Missing $fq" >&2; exit 1; }

    fastqc \
      -f fastq \
      -d raw_data/ \
      "$fq" \
      > "logs/${sample}_fastqc.log" 2>&1

    cutadapt \
      -g ACGATC \
      -a AGATCG \
      -o "raw_data/${sample}_R1.fastq" \
      "$fq" \
      > "logs/${sample}_cutadapt.log" 2>&1

    # Final corrected block in the historical code: retain 18–30 nt reads.
    trim_galore \
      --small_rna \
      --length 18 \
      --max_length 30 \
      --quality 30 \
      --stringency 3 \
      --phred33 \
      --cores 6 \
      --dont_gzip \
      -o clean_data \
      "raw_data/${sample}_R1.fastq" \
      > "logs/${sample}_trim_galore.log" 2>&1
done
