#!/usr/bin/env bash
set -euo pipefail
source ./00_config.sh
mkdir -p "$STAR_INDEX" logs

bowtie-build "$NCRNA_FASTA" "$NCRNA_INDEX" \
  > logs/bowtie_ncRNA_index.log 2>&1

STAR \
  --runThreadN 10 \
  --runMode genomeGenerate \
  --genomeDir "$STAR_INDEX" \
  --genomeFastaFiles "$GENOME" \
  --sjdbGTFfile "$GTF" \
  > logs/STAR_genomeGenerate.log 2>&1
