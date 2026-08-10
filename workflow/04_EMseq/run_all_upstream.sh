#!/usr/bin/env bash
set -euo pipefail

bash 01_fastp.sh
bash 02_prepare_reference.sh
bash 03_bismark_genome_preparation.sh
bash 04_bismark_align.sh
bash 05_samtools_sort_index_coverage.sh
bash 06_batmeth2_index.sh
bash 07_batmeth2_calmeth.sh
bash 08_merge_calmeth_outputs.sh
bash 09_prepare_feature_beds.sh
bash 10_batmeth2_methyGff_genes.sh
bash 11_batmeth2_methyGff_TEs.sh
bash 12_merge_feature_methylation.sh

echo "[DONE] Core EM-seq methylation processing completed."
echo "[INFO] Differential methylation scripts 13–16 require:"
echo "       software/run_DM_analysis.R"
echo "       samples_info.txt"
echo "       contrasts.txt"
