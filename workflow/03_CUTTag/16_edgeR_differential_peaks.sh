#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 ]]; then
    echo "Usage: bash $0 <H3K27ac|H3K4me3|H3K9ac|H3K27me3>" >&2
    exit 1
fi

MARK="$1"
source ./00_config.sh
PEAK_TYPE="$(get_peak_type "$MARK")"
OUT="$(result_dir "$MARK")"
mapfile -t SAMPLES < <(get_samples "$MARK")

if [[ "${#SAMPLES[@]}" -eq 0 ]]; then
    echo "[ERROR] No samples found for mark: $MARK" >&2
    exit 1
fi

mkdir -p "$OUT/diff_peaks" "$OUT/logs"

SAMPLES_TXT="$OUT/samples.txt"
CONTRASTS_TXT="$OUT/contrasts.txt"

: > "$SAMPLES_TXT"
for sample in "${SAMPLES[@]}"; do
    # e.g. H3K27ac_PN10_rep1 -> H3K27ac_PN10
    group="${sample%_rep*}"
    printf "%s\t%s\n" "$group" "$sample" >> "$SAMPLES_TXT"
done

groups=()
for stage in $(get_stages); do
    groups+=("${MARK}_${stage}")
done

: > "$CONTRASTS_TXT"
for ((i=0; i<${#groups[@]}; i++)); do
    for ((j=i+1; j<${#groups[@]}; j++)); do
        printf "%s %s\n" "${groups[$i]}" "${groups[$j]}" >> "$CONTRASTS_TXT"
    done
done

run_DE_analysis.pl \
    -m "$OUT/counts/merged_peaks.counts.matrix" \
    --method edgeR \
    -s "$SAMPLES_TXT" \
    --contrasts "$CONTRASTS_TXT" \
    -o "$OUT/diff_peaks" \
    --dispersion 0.1 \
    > "$OUT/logs/run_edgeR.log" 2>&1

pushd "$OUT/diff_peaks" >/dev/null

rm -f *count_matrix *.MA_n_Volcano.pdf *.Rscript 2>/dev/null || true

for file in merged_peaks.counts.matrix.*.DE_results; do
    [[ -e "$file" ]] || continue

    new="${file#merged_peaks.counts.matrix.}"

    printf "PeakID\tsampleA\tsampleB\tlog2FoldChange\tPValue\tPadj\n" > "$new"

    grep -v "^sampleA" "$file" \
        | awk 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$6,$7}' \
        >> "$new"

    rm "$file"
done

popd >/dev/null

echo "[DONE] edgeR differential peak analysis completed."
