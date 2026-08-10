#!/usr/bin/env bash
set -euo pipefail
source ./00_config.sh
mkdir -p ref

[[ -f "$GENOME" ]] || { echo "[ERROR] Missing $GENOME" >&2; exit 1; }
[[ -f "$GFF" ]] || { echo "[ERROR] Missing $GFF" >&2; exit 1; }

seqkit faidx "$GENOME"

gffread -T -o ref/temp.gtf "$GFF"
gtftk convert_ensembl -i ref/temp.gtf > "$GTF"
rm -f ref/temp.gtf

gtfToGenePred "$GTF" ref/genes.genePred
genePredToBed ref/genes.genePred ref/genes.bed12
