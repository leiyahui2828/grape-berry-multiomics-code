#!/usr/bin/env bash
set -euo pipefail
source ./00_config.sh

DE_DIR="$OUTDIR/RPF_quant/DESeq2_dir"
OUT="$OUTDIR/RPF_quant/all.RPF.DE_results"

files=("$DE_DIR"/*.DESeq2.DE_results)
[[ -e "${files[0]}" ]] || { echo "[ERROR] No DESeq2 result files found" >&2; exit 1; }

sed -n '1s/^/GeneID\t/p' "${files[0]}" > "$OUT"
for file in "${files[@]}"; do
    sed 1d "$file"
done >> "$OUT"

if [[ -f scripts/DESeq2_DEG_plot.R ]]; then
  Rscript scripts/DESeq2_DEG_plot.R \
    -i "$OUT" \
    -e "$OUTDIR/RPF_quant/matrix.TPM.not_cross_norm" \
    -o "$OUTDIR/RPF_quant/DESeq2_dir_plot" \
    --top 50 \
    --label_top 10
fi
