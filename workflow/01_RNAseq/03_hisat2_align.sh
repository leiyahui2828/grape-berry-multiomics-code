#!/usr/bin/env bash
set -euo pipefail

SAMPLES_FILE="${SAMPLES_FILE:-00_samples.txt}"

if [[ ! -f "$SAMPLES_FILE" ]]; then
    echo "[ERROR] Sample list not found: $SAMPLES_FILE" >&2
    exit 1
fi

mapfile -t SAMPLES < <(grep -v '^\s*$' "$SAMPLES_FILE")

mkdir -p aligned

for sample in "${SAMPLES[@]}"; do
    echo "[HISAT2] ${sample}"
    hisat2 \
        -x genome.hisat2 \
        -1 "clean_data/${sample}_1.fastp.fastq.gz" \
        -2 "clean_data/${sample}_2.fastp.fastq.gz" \
        --new-summary \
        --rna-strandness RF \
        -S "aligned/${sample}.sam" \
        > "aligned/${sample}.hisat2.log" 2>&1
done
