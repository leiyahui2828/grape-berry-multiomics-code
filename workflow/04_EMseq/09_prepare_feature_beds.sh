#!/usr/bin/env bash
set -euo pipefail
source ./00_config.sh

awk '$3=="gene"' "$GTF" \
  | awk '{print $1"\t"$4"\t"$5"\t"$10"\t"$7}' \
  | sed 's/"//g' \
  | sed 's/;//' \
  > ref/genes.bed5

[[ -f "$TE_GFF" ]] || { echo "[ERROR] Missing $TE_GFF" >&2; exit 1; }

grep 'TE' "$TE_GFF" \
  | awk '{print $1"\t"$4"\t"$5"\t"$10"\t"$7}' \
  | sed 's/"//g' \
  | sed 's/Motif://' \
  > ref/TE.bed5
