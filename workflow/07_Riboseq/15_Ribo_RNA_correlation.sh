#!/usr/bin/env bash
set -euo pipefail
source ./00_config.sh
DIR="$OUTDIR/RNA_Ribo_integration"

[[ -f scripts/ribo_rna_correlation.R ]] || {
    echo "[ERROR] Missing scripts/ribo_rna_correlation.R" >&2
    exit 1
}

Rscript scripts/ribo_rna_correlation.R \
  --ribo "$DIR/Ribo-seq.TPM.matrix" \
  --rna "$DIR/RNA-seq.TPM.matrix" \
  --sample_info "$DIR/sample_info_all.txt" \
  --outdir "$DIR/01.ribo_rna_correlation_result"
