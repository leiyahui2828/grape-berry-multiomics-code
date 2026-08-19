#!/usr/bin/env bash
set -euo pipefail
source ./00_config.sh

RUN_DE="scripts/DifferentialExpression/run_DE_analysis.pl"
[[ -f "$RUN_DE" ]] || { echo "[ERROR] Missing $RUN_DE" >&2; exit 1; }

mkdir -p "$OUTDIR/RPF_quant/DESeq2_dir"

perl "$RUN_DE" \
  --matrix "$OUTDIR/RPF_quant/matrix.counts.matrix" \
  --samples_file "$SAMPLE_INFO" \
  --output "$OUTDIR/RPF_quant/DESeq2_dir" \
  --contrasts "$CONTRASTS_FILE" \
  --method DESeq2
