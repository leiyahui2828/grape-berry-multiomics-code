#!/usr/bin/env bash
set -euo pipefail
bash 14_prepare_RNA_Ribo_integration.sh
bash 15_Ribo_RNA_correlation.sh
bash 16_Ribo_RNA_DE_compare.sh
bash 17_translation_efficiency.sh
bash 18_differential_translation.sh
echo "[DONE] RNA/Ribo integration completed."
