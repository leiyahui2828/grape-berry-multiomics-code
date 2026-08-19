#!/usr/bin/env bash
set -euo pipefail
source ./00_config.sh
DIR="$OUTDIR/RNA_Ribo_integration"

[[ -f scripts/ribo_rna_TE.R ]] || {
    echo "[ERROR] Missing scripts/ribo_rna_TE.R" >&2
    exit 1
}

Rscript scripts/ribo_rna_TE.R \
  --ribo_expr "$DIR/Ribo-seq.TPM.matrix" \
  --rna_expr "$DIR/RNA-seq.TPM.matrix" \
  --sample_info "$DIR/sample_info_all.txt" \
  --mean_cutoff 1 \
  --outdir "$DIR/03.TE_out"

if [[ -f scripts/log2TE_feature_compare.R && \
      -f ref/genome.five_prime_utr_geneID.GC_NMFE.tab ]]; then
  Rscript scripts/log2TE_feature_compare.R \
    --logte "$DIR/03.TE_out/log2TE_matrix.tsv" \
    --feature ref/genome.five_prime_utr_geneID.GC_NMFE.tab \
    -o "$DIR/03.TE_feature_compare" \
    --min_group_n 10
fi
