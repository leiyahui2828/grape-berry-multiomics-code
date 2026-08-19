#!/usr/bin/env bash
set -euo pipefail
source ./00_config.sh
DIR="$OUTDIR/RNA_Ribo_integration"

[[ -f scripts/ribo_rna_DE_compare_plot.R ]] || {
    echo "[ERROR] Missing scripts/ribo_rna_DE_compare_plot.R" >&2
    exit 1
}

Rscript scripts/ribo_rna_DE_compare_plot.R \
  --ribo_de "$DIR/Ribo-seq.DE_results" \
  --rna_de "$DIR/RNA-seq.DE_results" \
  --padj_cutoff 0.05 \
  --log2fc_cutoff 1 \
  --label_top 10 \
  --outdir "$DIR/02.DE_compare_result"
