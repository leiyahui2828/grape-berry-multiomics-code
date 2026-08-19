#!/usr/bin/env bash
set -euo pipefail
bash 01_prepare_reference.sh
bash 02_prepare_ncRNA.sh
bash 03_build_indexes.sh
bash 04_fastp.sh
bash 05_remove_ncRNA.sh
bash 06_clean_reads_QC.sh
bash 07_STAR_mapping.sh
bash 08_STAR_mapping_QC.sh
bash 09_riboWaltz_RPF_QC.sh
bash 10_featureCounts_CDS.sh
bash 11_build_RPF_expression_matrix.sh
bash 12_DESeq2_RPF.sh
bash 13_merge_RPF_DE.sh
echo "[DONE] Core Ribo-seq workflow completed."
