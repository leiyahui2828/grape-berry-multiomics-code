#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# RNA-seq differential expression analysis using Trinity run_DE_analysis.pl
#
# Input:
#   quanti/genes.counts.matrix
#   samples.txt
#   contrasts.txt
#
# Output:
#   DE_analysis/
#   DE_analysis/de_result.tsv
#
# Notes:
#   - This script preserves the analysis logic used in the original workflow:
#       Trinity v2.14.0 run_DE_analysis.pl + DESeq2
#   - Set TRINITY_HOME to your local Trinity installation before running.
###############################################################################

# -----------------------------
# 1. User-configurable paths
# -----------------------------
TRINITY_HOME="${TRINITY_HOME:-/path/to/trinityrnaseq-v2.14.0}"

COUNT_MATRIX="${COUNT_MATRIX:-quanti/genes.counts.matrix}"
SAMPLES_FILE="${SAMPLES_FILE:-samples.txt}"
CONTRASTS_FILE="${CONTRASTS_FILE:-contrasts.txt}"
OUTDIR="${OUTDIR:-DE_analysis}"

RUN_DE="${TRINITY_HOME}/Analysis/DifferentialExpression/run_DE_analysis.pl"

# -----------------------------
# 2. Check required inputs
# -----------------------------
for f in "$COUNT_MATRIX" "$SAMPLES_FILE" "$CONTRASTS_FILE"; do
    if [[ ! -f "$f" ]]; then
        echo "[ERROR] Required file not found: $f" >&2
        exit 1
    fi
done

if [[ ! -f "$RUN_DE" ]]; then
    echo "[ERROR] Trinity run_DE_analysis.pl not found:" >&2
    echo "        $RUN_DE" >&2
    echo "[INFO] Please set TRINITY_HOME, for example:" >&2
    echo "       export TRINITY_HOME=/home/user/trinityrnaseq-v2.14.0" >&2
    exit 1
fi

mkdir -p "$OUTDIR"

# -----------------------------
# 3. Differential expression
# -----------------------------
echo "[INFO] Running DESeq2 differential expression analysis..."

perl "$RUN_DE" \
    --matrix "$COUNT_MATRIX" \
    --method DESeq2 \
    --samples_file "$SAMPLES_FILE" \
    --contrasts "$CONTRASTS_FILE" \
    --output "$OUTDIR"

# -----------------------------
# 4. Merge all DESeq2 results
# -----------------------------
echo "[INFO] Merging DESeq2 comparison tables..."

mapfile -t DE_FILES < <(
    find "$OUTDIR" -maxdepth 1 -type f -name '*DE_results' | sort
)

if [[ "${#DE_FILES[@]}" -eq 0 ]]; then
    echo "[ERROR] No *DE_results files were generated in $OUTDIR" >&2
    exit 1
fi

awk 'FNR==1 && NR!=1 {next} {print}' "${DE_FILES[@]}" \
    > "${OUTDIR}/de_result.tsv"

# The original workflow added 'gene_id' to the first header column.
sed -i '1s/^/gene_id\t/' "${OUTDIR}/de_result.tsv"

echo "[DONE] Differential expression analysis completed."
echo "[OUTPUT] ${OUTDIR}/de_result.tsv"
