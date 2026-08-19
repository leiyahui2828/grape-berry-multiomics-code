#!/usr/bin/env bash
set -euo pipefail
source ./00_config.sh

FILES=()
for f in ref/rRNA.fasta ref/tRNA.fasta ref/snRNA.fasta ref/snoRNA.fasta ref/miRNA.fasta; do
    [[ -f "$f" ]] && FILES+=("$f")
done

if [[ "${#FILES[@]}" -eq 0 ]]; then
    echo "[ERROR] No grape ncRNA FASTA files found under ref/." >&2
    exit 1
fi

cat "${FILES[@]}" | seqkit seq --rna2dna > "$NCRNA_FASTA"
