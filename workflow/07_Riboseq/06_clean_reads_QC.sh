#!/usr/bin/env bash
set -euo pipefail
source ./00_config.sh

mkdir -p "$OUTDIR/qc/length_distribution"

seqkit stats -j 10 \
  "$OUTDIR"/ncrna_filter/*.fq.gz \
  > "$OUTDIR/qc/clean_reads.stat.tsv"

if [[ -f scripts/fastq_length_distribution.py ]]; then
  python scripts/fastq_length_distribution.py \
    -i "$OUTDIR"/ncrna_filter/*.fq.gz \
    -o "$OUTDIR/qc/length_distribution" \
    -t 6
fi
