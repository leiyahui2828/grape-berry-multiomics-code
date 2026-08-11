#!/usr/bin/env bash
set -euo pipefail

mkdir -p miRNA_family/phylogeny

cp miRNA_family/all.hairpin.fa miRNA_family/phylogeny/
cd miRNA_family/phylogeny

muscle -align all.hairpin.fa -output all.hairpin.mfa

trimal -in all.hairpin.mfa -out all.hairpin.phylip -phylip_paml
trimal -in all.hairpin.mfa -out all.hairpin.meg -mega

treebest nj all.hairpin.mfa -W -t jtt -b 100 \
  > all.hairpin.treebest.out

sed 's/:\([0-9.]\+\)\[&&NHX:B=\([0-9]\+\)\]/\2:\1/' \
  all.hairpin.treebest.out \
  | awk '{printf $0}' \
  > all.hairpin.treebest.nwk

fasttree \
  -out all.hairpin.fasttree.nwk \
  all.hairpin.mfa

iqtree \
  -s all.hairpin.mfa \
  -m MFP \
  -bb 1000 \
  -nt 10 \
  -pre all.hairpin.iqtree
