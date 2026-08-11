#!/usr/bin/env bash
set -euo pipefail

bash 01_prepare_miRBase.sh
bash 02_prepare_Rfam.sh
bash 03_prepare_reference.sh
bash 04_fastqc_trim.sh
bash 05_collapse_reads.sh
bash 06_length_statistics.sh
bash 07_map_genome.sh
bash 08_map_annotation_databases.sh
bash 09_merge_sRNA_annotation.sh
bash 10_prepare_novel_miRNA_reads.sh
bash 11_miRDP2_predict.sh
bash 12_rename_novel_miRNA.sh

echo "[DONE] Core small-RNA/miRNA processing completed."
echo "[INFO] Scripts 13–20 cover sequence analysis, family analysis,"
echo "       expression, differential analysis and target prediction."
