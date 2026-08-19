#!/usr/bin/env bash
set -euo pipefail
source ./00_config.sh
mapfile -t SAMPLES < <(get_samples)

[[ -f scripts/run_riboWaltz.R ]] || {
    echo "[ERROR] Missing scripts/run_riboWaltz.R" >&2
    exit 1
}

mkdir -p "$OUTDIR/riboWaltz/bam_dir"

for sample in "${SAMPLES[@]}"; do
    bam="$OUTDIR/STAR_mapping/${sample}/${sample}.Aligned.toTranscriptome.out.bam"
    [[ -f "$bam" ]] || { echo "[ERROR] Missing $bam" >&2; exit 1; }

    ln -sf "$(realpath "$bam")" \
      "$OUTDIR/riboWaltz/bam_dir/${sample}.bam"
done

Rscript scripts/run_riboWaltz.R \
  --bam_dir "$OUTDIR/riboWaltz/bam_dir/" \
  --gtf "$GTF" \
  --fasta "$TRANSCRIPTS" \
  --sample_info "$SAMPLE_INFO" \
  --length_filter_mode custom \
  --min_len "$RPF_MIN_LEN" \
  --max_len "$RPF_MAX_LEN"
