#!/usr/bin/env bash
set -euo pipefail
source ./00_config.sh

[[ -f scripts/STAR_log2stat.py ]] || {
    echo "[ERROR] Missing scripts/STAR_log2stat.py" >&2
    exit 1
}

mkdir -p "$OUTDIR/qc"

python scripts/STAR_log2stat.py \
  -o "$OUTDIR/qc/STAR_mapping.stat.tsv" \
  "$OUTDIR"/STAR_mapping/*/*Log.final.out
