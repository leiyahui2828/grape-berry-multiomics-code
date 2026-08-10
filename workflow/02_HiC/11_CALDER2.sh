#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# 11_CALDER2.sh
# Hi-C compartment / subcompartment analysis using CALDER2
#
# Preserves the parameters recorded in the supplied hic.code:
#   Rscript scripts/calder
#   -t mcool
#   -b 40000
#   -p 4
#   -g others
#   -f ref/gene_density.bed
#   --chr_prefix Chr
###############################################################################

SAMPLES_FILE="${SAMPLES_FILE:-00_samples.txt}"
CALDER_SCRIPT="${CALDER_SCRIPT:-scripts/calder}"
RESOLUTION=40000
THREADS=4

if [[ ! -f "$SAMPLES_FILE" ]]; then
    echo "[ERROR] Sample list not found: $SAMPLES_FILE" >&2
    exit 1
fi

if [[ ! -f "$CALDER_SCRIPT" ]]; then
    echo "[ERROR] CALDER2 wrapper not found: $CALDER_SCRIPT" >&2
    echo "[INFO] Copy the exact historical scripts/calder used in the original analysis." >&2
    exit 1
fi

if [[ ! -f "ref/gene_density.bed" ]]; then
    echo "[ERROR] Missing reference file: ref/gene_density.bed" >&2
    exit 1
fi

mapfile -t SAMPLES < <(grep -v '^[[:space:]]*$' "$SAMPLES_FILE")

mkdir -p "Compartments/CALDER2/${RESOLUTION}" logs

for sample in "${SAMPLES[@]}"; do
    stage="${sample%%_Rep*}"
    matrix="Matrix/Final/${stage}.mcool"
    outdir="Compartments/CALDER2/${RESOLUTION}/${stage}"
    logfile="logs/${stage}_${RESOLUTION}_calder2.log"

    if [[ ! -f "$matrix" ]]; then
        echo "[ERROR] Missing matrix: $matrix" >&2
        exit 1
    fi

    mkdir -p "$outdir"

    echo "[CALDER2] ${stage}"

    Rscript "$CALDER_SCRIPT" \
        -i "$matrix" \
        -t mcool \
        -b "$RESOLUTION" \
        -p "$THREADS" \
        -o "$outdir" \
        -g others \
        -f ref/gene_density.bed \
        --chr_prefix Chr \
        > "$logfile" 2>&1
done

echo "[DONE] CALDER2 compartment/subcompartment analysis completed."
