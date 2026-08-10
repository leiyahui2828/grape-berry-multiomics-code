#!/usr/bin/env bash
set -euo pipefail

mkdir -p quanti logs

HELPER="software/RunFeatureCounts/abundance_estimates_to_matrix.pl"

if [[ ! -f "$HELPER" ]]; then
    echo "[ERROR] Missing original helper script: $HELPER" >&2
    exit 1
fi

cat > quant_files.txt <<'EOF'
quanti/PN4_rep1.count
quanti/PN4_rep2.count
quanti/PN4_rep3.count
quanti/PN7_rep1.count
quanti/PN7_rep2.count
quanti/PN7_rep3.count
quanti/PN10_rep1.count
quanti/PN10_rep2.count
quanti/PN10_rep3.count
quanti/PN15_rep1.count
quanti/PN15_rep2.count
quanti/PN15_rep3.count
quanti/PN18_rep1.count
quanti/PN18_rep2.count
quanti/PN18_rep3.count
EOF

perl "$HELPER" \
    --est_method featureCounts \
    --out_prefix quanti/genes \
    --quant_files quant_files.txt \
    > logs/abundance_estimates_to_matrix.log 2>&1

# Add feature ID column label, exactly as in the supplied analysis.
sed -i '1s/^/gene_id/' quanti/genes.counts.matrix
sed -i '1s/^/gene_id/' quanti/genes.TMM.TPM.matrix
