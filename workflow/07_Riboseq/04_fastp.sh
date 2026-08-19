#!/usr/bin/env bash
set -euo pipefail
source ./00_config.sh
mapfile -t SAMPLES < <(get_samples)

mkdir -p "$OUTDIR/clean_data" "$OUTDIR/qc/fastp" "$OUTDIR/logs"

for sample in "${SAMPLES[@]}"; do
    fq="$RAW_DIR/$(get_fastq "$sample")"
    [[ -f "$fq" ]] || { echo "[ERROR] Missing $fq" >&2; exit 1; }

    fastp \
      -i "$fq" \
      -o "$OUTDIR/clean_data/${sample}.clean.fq.gz" \
      --trim_front1 0 \
      --adapter_sequence "$RIBO_ADAPTER" \
      --length_required "$FASTP_MIN_LEN" \
      --length_limit "$FASTP_MAX_LEN" \
      --html "$OUTDIR/qc/fastp/${sample}.fastp.html" \
      --json "$OUTDIR/qc/fastp/${sample}.fastp.json" \
      --thread "$FASTP_THREADS" \
      > "$OUTDIR/logs/${sample}.fastp.log" 2>&1
done

if [[ -f scripts/fastp_json2tab.R ]]; then
    ls "$OUTDIR"/qc/fastp/*.json \
      | awk -F "/" '{n=$NF; sub(/\.fastp\.json$/,"",n); print n"\t"$0}' \
      > "$OUTDIR/qc/fastp.json.list"

    Rscript scripts/fastp_json2tab.R \
      "$OUTDIR/qc/fastp.json.list" \
      "$OUTDIR/qc/fastp_stat.tsv"
fi
