#!/usr/bin/env bash
set -euo pipefail

SAMPLES_FILE="${SAMPLES_FILE:-00_samples.txt}"

if [[ ! -f "$SAMPLES_FILE" ]]; then
    echo "[ERROR] Sample list not found: $SAMPLES_FILE" >&2
    exit 1
fi

mapfile -t SAMPLES < <(grep -v '^[[:space:]]*$' "$SAMPLES_FILE")

mkdir -p Matrix/Raw logs

for sample in "${SAMPLES[@]}"; do
    echo "[hicBuildMatrix] ${sample}"

    mkdir -p "Matrix/Raw/${sample}"

    hicBuildMatrix \
        --samFiles \
            "aligned/${sample}_R1.bam" \
            "aligned/${sample}_R2.bam" \
        --genomeAssembly ref/genome.fa \
        --outBam "Matrix/Raw/${sample}.bam" \
        --outFileName "Matrix/Raw/${sample}.cool" \
        --QCfolder "Matrix/Raw/${sample}" \
        --restrictionCutFile ref/rest_site.bed \
        --binSize 1000 \
        --threads 6 \
        --minDistance 500 \
        --maxLibraryInsertSize 1500 \
        --inputBufferSize 1000000 \
        --restrictionSequence GATC \
        --danglingSequence GATC \
        > "logs/${sample}_hicBuildMatrix.log" 2>&1
done
