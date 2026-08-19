#!/usr/bin/env bash
set -euo pipefail
source ./00_config.sh

HELPER="scripts/RunFeatureCounts/abundance_estimates_to_matrix.pl"
[[ -f "$HELPER" ]] || { echo "[ERROR] Missing $HELPER" >&2; exit 1; }

mkdir -p "$OUTDIR/RPF_quant"

ls "$OUTDIR"/featureCounts/*.count \
  > "$OUTDIR/RPF_quant/genes.quant_files.txt"

(
  cd "$OUTDIR/RPF_quant"
  perl "../../$HELPER" \
    --est_method featureCounts \
    --quant_files genes.quant_files.txt \
    --cross_sample_norm none
)

if [[ -f scripts/expression_QC.R ]]; then
  Rscript scripts/expression_QC.R \
    --expr "$OUTDIR/RPF_quant/matrix.TPM.not_cross_norm" \
    --group "$SAMPLE_INFO" \
    --outdir "$OUTDIR/RPF_quant/expression_QC"
fi
