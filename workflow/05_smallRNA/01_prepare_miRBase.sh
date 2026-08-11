#!/usr/bin/env bash
set -euo pipefail
source ./00_config.sh

mkdir -p "$MIRBASE_DIR"

MATURE="${MATURE:-$MIRBASE_DIR/mature.fa}"
HAIRPIN="${HAIRPIN:-$MIRBASE_DIR/hairpin.fa}"
SPC="$SPECIES_PREFIX"

for f in "$MATURE" "$HAIRPIN"; do
    [[ -f "$f" ]] || { echo "[ERROR] Missing miRBase file: $f" >&2; exit 1; }
done

grep "^>${SPC}" "$MATURE" | sed 's/^>//' | awk '{print $1}' \
  > "$MIRBASE_DIR/mature.${SPC}.id"

seqkit grep -f "$MIRBASE_DIR/mature.${SPC}.id" "$MATURE" \
  | awk '{print $1}' \
  > "$MIRBASE_DIR/mature.${SPC}.fa.tmp"

sed '/^>/!s/U/T/gi' "$MIRBASE_DIR/mature.${SPC}.fa.tmp" \
  > "$MIRBASE_DIR/mature.${SPC}.fa"

grep "^>${SPC}" "$HAIRPIN" | sed 's/^>//' | awk '{print $1}' \
  > "$MIRBASE_DIR/hairpin.${SPC}.id"

seqkit grep -f "$MIRBASE_DIR/hairpin.${SPC}.id" "$HAIRPIN" \
  | awk '{print $1}' \
  > "$MIRBASE_DIR/hairpin.${SPC}.fa.tmp"

sed '/^>/!s/U/T/gi' "$MIRBASE_DIR/hairpin.${SPC}.fa.tmp" \
  > "$MIRBASE_DIR/hairpin.${SPC}.fa"

bowtie-build \
  "$MIRBASE_DIR/hairpin.${SPC}.fa" \
  "$MIRBASE_DIR/hairpin.${SPC}"
