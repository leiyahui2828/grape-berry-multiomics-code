#!/usr/bin/env bash
set -euo pipefail

SAMPLES_FILE="${SAMPLES_FILE:-00_samples.txt}"

if [[ ! -f "$SAMPLES_FILE" ]]; then
    echo "[ERROR] Sample list not found: $SAMPLES_FILE" >&2
    exit 1
fi

mapfile -t SAMPLES < <(grep -v '^\s*$' "$SAMPLES_FILE")

mkdir -p quanti

HELPER="software/RunFeatureCounts/run-featurecounts.R"

if [[ ! -f "$HELPER" ]]; then
    echo "[ERROR] Missing original helper script: $HELPER" >&2
    exit 1
fi

for sample in "${SAMPLES[@]}"; do
    echo "[featureCounts] ${sample}"
    Rscript "$HELPER" \
        -b "aligned/${sample}.bam" \
        -g ref/genes.gtf \
        -s 2 \
        -o "quanti/${sample}"
done
