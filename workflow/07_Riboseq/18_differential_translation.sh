#!/usr/bin/env bash
set -euo pipefail
source ./00_config.sh
DIR="$OUTDIR/RNA_Ribo_integration"

DTEG="scripts/translational_regulation/DTEG_modified.R"
PLOT="scripts/DESeq2_plot_summary.R"

[[ -f "$DTEG" ]] || { echo "[ERROR] Missing $DTEG" >&2; exit 1; }

mkdir -p "$DIR/04.DTEG"

while read -r treat control; do
    [[ -z "${treat:-}" ]] && continue
    [[ "$treat" =~ ^# ]] && continue

    comp="${treat}_vs_${control}"
    outdir="$DIR/04.DTEG/${comp}"

    Rscript "$DTEG" \
      --ribo_count "$DIR/Ribo-seq.counts.matrix" \
      --rna_count "$DIR/RNA-seq.counts.matrix" \
      --sample_info "$DIR/sample_info_all.txt" \
      --outdir "$outdir" \
      --batch 0 \
      --control "$control" \
      --treat "$treat" \
      --log2fc_cutoff 1

    if [[ -f "$PLOT" && -f "$outdir/fold_changes/deltaRNA_addsample.txt" ]]; then
      Rscript "$PLOT" \
        -i "$outdir/fold_changes/deltaRNA_addsample.txt" \
        -o "${outdir}_plot" \
        --padj_cutoff 0.05 \
        --log2fc_cutoff 1
    fi
done < "$CONTRASTS_FILE"
