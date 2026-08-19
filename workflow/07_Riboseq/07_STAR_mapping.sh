#!/usr/bin/env bash
set -euo pipefail
source ./00_config.sh
mapfile -t SAMPLES < <(get_samples)

mkdir -p "$OUTDIR/STAR_mapping" "$OUTDIR/logs"

for sample in "${SAMPLES[@]}"; do
    fq="$OUTDIR/ncrna_filter/${sample}.ncrna.unmapped.fq.gz"
    sample_dir="$OUTDIR/STAR_mapping/${sample}"
    prefix="${sample_dir}/${sample}."
    mkdir -p "$sample_dir"

    STAR \
      --runThreadN "$STAR_THREADS" \
      --genomeDir "$STAR_INDEX" \
      --readFilesIn "$fq" \
      --readFilesCommand zcat \
      --outFileNamePrefix "$prefix" \
      --outSAMtype BAM SortedByCoordinate \
      --quantMode TranscriptomeSAM GeneCounts \
      --outFilterMultimapNmax "$STAR_MULTIMAP_NMAX" \
      --outFilterMismatchNmax "$STAR_MISMATCH_NMAX" \
      --outSAMattributes Standard \
      > "$OUTDIR/logs/${sample}.STAR.log" 2>&1
done
