#!/usr/bin/env bash
set -euo pipefail

[[ -f scripts/expression_draw1.R ]] || {
    echo "[ERROR] Missing scripts/expression_draw1.R" >&2
    exit 1
}

mkdir -p miRNA_exp/sr_exp_plot

ln -sf ../miRNAs_expressed.TPM.txt \
  miRNA_exp/sr_exp_plot/mirna_exp.TPM.txt

ln -sf ../sample.txt \
  miRNA_exp/sr_exp_plot/sample.txt

Rscript scripts/expression_draw1.R \
  miRNA_exp/sr_exp_plot/mirna_exp.TPM.txt \
  miRNA_exp/sr_exp_plot/sample.txt \
  miRNA_exp/sr_exp_plot/mirna_exp.TPM
