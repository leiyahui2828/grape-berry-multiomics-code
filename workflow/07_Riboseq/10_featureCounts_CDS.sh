#!/usr/bin/env bash
set -euo pipefail
source ./00_config.sh
mapfile -t SAMPLES < <(get_samples)

RUN_FEATURECOUNTS="scripts/RunFeatureCounts/run-featurecounts.R"
[[ -f "$RUN_FEATURECOUNTS" ]] || {
    echo "[ERROR] Missing $RUN_FEATURECOUNTS" >&2
    exit 1
}

mkdir -p "$OUTDIR/featureCounts" "$OUTDIR/logs"

for sample in "${SAMPLES[@]}"; do
    bam="$OUTDIR/STAR_mapping/${sample}/${sample}.Aligned.sortedByCoord.out.bam"

    Rscript "$RUN_FEATURECOUNTS" \
      -b "$bam" \
      -s "$RPF_STRAND" \
      -g "$GTF" \
      -i "$RPF_IS_PAIRED" \
      -f "$RPF_FEATURE" \
      -a "$RPF_ATTRIBUTE" \
      -o "$OUTDIR/featureCounts/${sample}" \
      > "$OUTDIR/logs/${sample}.featureCounts.log" 2>&1
done
