#!/usr/bin/env bash
set -euo pipefail

bash 01_fastp.sh
bash 02_bwa_index.sh
bash 03_hicFindRestSite.sh
bash 04_bwa_align.sh
bash 05_sam2bam.sh
bash 06_hicBuildMatrix.sh
bash 07_prepare_stage_matrices.sh
bash 08_hicNormalize.sh
bash 09_cooler_zoomify.sh
bash 10_estimate_resolution.sh
bash 11_CALDER2.sh
bash 12_hicFindTADs.sh
bash 13_hicDifferentialTAD.sh

echo "[DONE] Hi-C upstream processing completed."
echo "[NOTE] Run 10_estimate_resolution.sh separately after restoring scripts/hicres.py."
