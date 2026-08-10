#!/usr/bin/env bash
set -euo pipefail

SAMPLES_FILE="${SAMPLES_FILE:-00_samples.txt}"

if [[ ! -f "$SAMPLES_FILE" ]]; then
    echo "[ERROR] Sample list not found: $SAMPLES_FILE" >&2
    exit 1
fi

mapfile -t SAMPLES < <(grep -v '^[[:space:]]*$' "$SAMPLES_FILE")

mkdir -p clean_data logs

ADAPTER_R1="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
ADAPTER_R2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

for sample in "${SAMPLES[@]}"; do
    echo "[fastp] ${sample}"

    fastp \
        -i "raw_data/${sample}_R1.fastq.gz" \
        -o "clean_data/${sample}_R1.fastq.gz" \
        -I "raw_data/${sample}_R2.fastq.gz" \
        -O "clean_data/${sample}_R2.fastq.gz" \
        --html "clean_data/${sample}_fastp.html" \
        --json "clean_data/${sample}_fastp.json" \
        --thread 3 \
        --adapter_sequence "$ADAPTER_R1" \
        --adapter_sequence_r2 "$ADAPTER_R2" \
        --qualified_quality_phred 20 \
        --length_required 50 \
        > "logs/${sample}_fastp_quality_control.log" 2>&1
done
