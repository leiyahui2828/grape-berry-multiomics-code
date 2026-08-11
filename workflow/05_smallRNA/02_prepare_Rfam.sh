#!/usr/bin/env bash
set -euo pipefail
source ./00_config.sh

mkdir -p "$RFAM_DIR"

FAMILY="${FAMILY:-$RFAM_DIR/family.txt}"
FASTADIR="${FASTADIR:-$RFAM_DIR/fasta_files}"

[[ -f "$FAMILY" ]] || {
    echo "[ERROR] Missing Rfam family table: $FAMILY" >&2
    echo "[INFO] The historical workflow downloaded Rfam CURRENT." >&2
    exit 1
}

cut -f 1,2,19 "$FAMILY" > "$RFAM_DIR/family.txt1"

awk '$0 !~ /Gene; miRNA/' "$RFAM_DIR/family.txt1" \
  > "$RFAM_DIR/family.txt.rmMIR"

cut -f 1 "$RFAM_DIR/family.txt.rmMIR" | while read -r family; do
    zcat "$FASTADIR/${family}.fa.gz" | sed "s/^>/>${family}:/"
done > "$RFAM_DIR/Rfam_rmMIR.fasta"

bowtie-build \
  "$RFAM_DIR/Rfam_rmMIR.fasta" \
  "$RFAM_DIR/Rfam_rmMIR.fasta"

awk '$0 ~ /Gene; miRNA/' "$RFAM_DIR/family.txt1" \
  > "$RFAM_DIR/family.txt.MIR"

cut -f 1 "$RFAM_DIR/family.txt.MIR" | while read -r family; do
    zcat "$FASTADIR/${family}.fa.gz" | sed "s/^>/>${family}:/"
done > "$RFAM_DIR/Rfam_MIR.fasta"

makeblastdb \
  -in "$RFAM_DIR/Rfam_MIR.fasta" \
  -dbtype nucl
