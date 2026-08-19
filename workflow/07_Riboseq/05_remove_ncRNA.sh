#!/usr/bin/env bash
set -euo pipefail
source ./00_config.sh
mapfile -t SAMPLES < <(get_samples)

mkdir -p "$OUTDIR/ncrna_filter" "$OUTDIR/logs"

echo -e "Sample\tTotal_reads\tMapped_ncRNA_reads\tMapped_ncRNA_rate_percent\tUnmapped_reads\tUnmapped_rate_percent" \
  > "$OUTDIR/ncrna_filter/ncrna_mapping.stat"

for sample in "${SAMPLES[@]}"; do
    clean="$OUTDIR/clean_data/${sample}.clean.fq.gz"
    unmap="$OUTDIR/ncrna_filter/${sample}.ncrna.unmapped.fq"
    log="$OUTDIR/logs/${sample}.bowtie.ncrna.log"
    stat="$OUTDIR/ncrna_filter/${sample}.bowtie.ncrna.stat"

    bowtie \
      -q \
      -v "$NCRNA_MISMATCH" \
      -p "$BOWTIE_THREADS" \
      --un "$unmap" \
      -x "$NCRNA_INDEX" \
      "$clean" \
      > /dev/null \
      2> "$log"

    gzip -f "$unmap"

    awk -v sample="$sample" '
      BEGIN{OFS="\t";total="NA";mapped="NA";maprate="NA";unmapped="NA";unmaprate="NA"}
      /reads processed:/ {total=$NF}
      /at least one alignment:/ {mapped=$(NF-1);maprate=$NF;gsub(/[()%]/,"",maprate)}
      /failed to align:/ {unmapped=$(NF-1);unmaprate=$NF;gsub(/[()%]/,"",unmaprate)}
      END{print sample,total,mapped,maprate,unmapped,unmaprate}
    ' "$log" > "$stat"

    cat "$stat" >> "$OUTDIR/ncrna_filter/ncrna_mapping.stat"
done
