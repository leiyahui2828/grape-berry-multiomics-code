# EM-seq / DNA methylation analysis

This directory organizes the commands in the supplied
`甲基化组学数据上游分析(2).txt` as ordinary Bash scripts.

Design:

```text
PN4, PN7, PN10, PN15, PN18
3 biological replicates per stage
15 paired-end libraries in total
```

## Workflow

```text
FASTQ
  ↓
fastp
  ↓
Bismark alignment
  ↓
SAMtools sorting / indexing / coverage
  ↓
BatMeth2 index + calmeth
  ↓
site/window methylation tables
  ↓
BatMeth2 methyGff
  ├── genes/promoters
  └── transposable elements
  ↓
DMC / DMR / promoter-DMP analysis
```

## Files

```text
04_EMseq/
├── README.md
├── METHOD_AUDIT.md
├── 00_samples.tsv
├── 00_config.sh
├── 01_fastp.sh
├── 02_prepare_reference.sh
├── 03_bismark_genome_preparation.sh
├── 04_bismark_align.sh
├── 05_samtools_sort_index_coverage.sh
├── 06_batmeth2_index.sh
├── 07_batmeth2_calmeth.sh
├── 08_merge_calmeth_outputs.sh
├── 09_prepare_feature_beds.sh
├── 10_batmeth2_methyGff_genes.sh
├── 11_batmeth2_methyGff_TEs.sh
├── 12_merge_feature_methylation.sh
├── 13_DMC_analysis.sh
├── 14_DMR_analysis.sh
├── 15_DMP_promoter_analysis.sh
├── 16_merge_DMRs.sh
├── run_all_upstream.sh
├── software/
│   └── README.md
└── ref/
    └── README.md
```

## Historical key parameters

```text
fastp:
  --thread 5

Bismark:
  -p 4
  paired-end

BatMeth2 calmeth:
  -Q 20
  --remove_dup
  --coverage 4
  -nC 1
  --Regions 600
  --step 50000

BatMeth2 methyGff:
  --coverage 4
  -nC 1
  --body --promoter --GENE
  --distance 2000
  --step 0.01

DMR:
  window = 600 bp
  step = 200 bp
  contexts = CpG / CHG / CHH
```

## Required local files not present in the uploaded code text

```text
software/batmeth2merge.R
software/run_DM_analysis.R
samples_info.txt
contrasts.txt
```

Restore the exact historical versions before the GitHub/Zenodo release.
