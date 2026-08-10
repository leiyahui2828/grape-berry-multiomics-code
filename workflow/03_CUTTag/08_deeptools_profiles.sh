#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 ]]; then
    echo "Usage: bash $0 <H3K27ac|H3K4me3|H3K9ac|H3K27me3>" >&2
    exit 1
fi

MARK="$1"
source ./00_config.sh
PEAK_TYPE="$(get_peak_type "$MARK")"
OUT="$(result_dir "$MARK")"
mapfile -t SAMPLES < <(get_samples "$MARK")

if [[ "${#SAMPLES[@]}" -eq 0 ]]; then
    echo "[ERROR] No samples found for mark: $MARK" >&2
    exit 1
fi

mkdir -p "$OUT/deeptools" "$OUT/logs"

BWS=()
LABELS=()
for sample in "${SAMPLES[@]}"; do
    BWS+=("$OUT/clean_bams/${sample}.bw")
    LABELS+=("$sample")
done

multiBigwigSummary BED-file \
    --bwfiles "${BWS[@]}" \
    --labels "${LABELS[@]}" \
    --BED "$TSS_TES_BED" \
    -o "$OUT/clean_bams/bw_summary.gz" \
    -p 10 \
    > "$OUT/logs/summarize_multiple_bigwigs.log" 2>&1

plotCorrelation \
    -in "$OUT/clean_bams/bw_summary.gz" \
    -o "$OUT/deeptools/sample_correlation.pdf" \
    --corMethod pearson \
    --whatToPlot heatmap \
    -min 0.5 \
    --plotTitle "Pearson Correlation of Samples" \
    --outFileCorMatrix "$OUT/deeptools/sample_correlation.tab" \
    > "$OUT/logs/generate_correlation_plots.log" 2>&1

computeMatrix reference-point \
    -S "${BWS[@]}" \
    --samplesLabel "${LABELS[@]}" \
    -R "$TSS_TES_BED" "$TSS_TES_SHUFFLE_BED" \
    --referencePoint TSS \
    -b 5000 -a 5000 \
    --binSize 50 \
    -o "$OUT/deeptools/tss_matrix.gz" \
    -p 10 \
    > "$OUT/logs/plot_heatmap_reference_point.log" 2>&1

plotHeatmap \
    -m "$OUT/deeptools/tss_matrix.gz" \
    -o "$OUT/deeptools/tss_heatmap.pdf" \
    --missingDataColor 0.5 \
    >> "$OUT/logs/plot_heatmap_reference_point.log" 2>&1

computeMatrix scale-regions \
    -S "${BWS[@]}" \
    --samplesLabel "${LABELS[@]}" \
    -R "$TSS_TES_BED" "$TSS_TES_SHUFFLE_BED" \
    --regionBodyLength 4000 \
    -b 2000 -a 2000 \
    --binSize 50 \
    -o "$OUT/deeptools/tss_tes_matrix.gz" \
    -p 10 \
    > "$OUT/logs/plot_heatmap_scale_regions.log" 2>&1

plotHeatmap \
    -m "$OUT/deeptools/tss_tes_matrix.gz" \
    -o "$OUT/deeptools/tss_tes_heatmap.pdf" \
    --missingDataColor 0.5 \
    >> "$OUT/logs/plot_heatmap_scale_regions.log" 2>&1
