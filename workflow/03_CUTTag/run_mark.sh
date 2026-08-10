#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 ]]; then
    echo "Usage: bash run_mark.sh <H3K27ac|H3K4me3|H3K9ac|H3K27me3>" >&2
    exit 1
fi

MARK="$1"

bash 01_fastp.sh "$MARK"

# Reference index is shared among all marks. Run only if absent.
if [[ ! -f ref/genome.1.bt2 ]]; then
    bash 02_bowtie2_build.sh
fi

bash 03_bowtie2_align.sh "$MARK"
bash 04_samtools_sort.sh "$MARK"
bash 05_filter_bam.sh "$MARK"
bash 06_qc.sh "$MARK"
bash 07_bigwig.sh "$MARK"
bash 08_deeptools_profiles.sh "$MARK"
bash 09_MACS3_peak_calling.sh "$MARK"
bash 10_filter_peaks_qscore.sh "$MARK"
bash 11_select_reproducible_peaks.sh "$MARK"
bash 12_peak_annotation_uropa.sh "$MARK"
bash 13_merge_consensus_peaks.sh "$MARK"
bash 14_featureCounts.sh "$MARK"
bash 15_build_count_matrix.sh "$MARK"
bash 16_edgeR_differential_peaks.sh "$MARK"
bash 17_extract_peak_sequence.sh "$MARK"

echo "[DONE] CUT&Tag workflow completed for $MARK"
