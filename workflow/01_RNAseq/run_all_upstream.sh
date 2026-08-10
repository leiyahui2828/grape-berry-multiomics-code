#!/usr/bin/env bash
set -euo pipefail

bash 01_fastp.sh
bash 02_build_hisat2_index.sh
bash 03_hisat2_align.sh
bash 04_samtools_sort_index.sh
bash 05_combine_qc_reports.sh
bash 06_featureCounts.sh
bash 07_build_expression_matrix.sh
bash 08_DESeq2.sh

echo "[DONE] RNA-seq upstream workflow completed."
