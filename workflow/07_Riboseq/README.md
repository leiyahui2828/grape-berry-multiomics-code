# Ribo-seq analysis workflow

This directory is a **grape-adapted Ribo-seq workflow** organized from the
uploaded Arabidopsis Ribo-seq reference code.

Design:

```text
PN4, PN7, PN10, PN15, PN18
3 biological replicates per stage
15 single-end Ribo-seq libraries
```

Ribo-seq sample IDs are matched one-to-one with RNA-seq:

```text
PN4_rep1 ... PN18_rep3
```

## Important

Unlike the already organized RNA-seq, Hi-C, CUT&Tag and EM-seq modules, this is
not yet a reconstruction of an executed grape Ribo-seq log. It is an adaptation
of the supplied Arabidopsis reference workflow. Therefore adapter, strandness,
RPF-length filters, ncRNA database and other library-specific parameters must be
verified against the actual grape analysis before the final GitHub/Zenodo
release.

## Directory

```text
07_Riboseq/
├── README.md
├── METHOD_AUDIT.md
├── 00_config.sh
├── 00_samples.tsv
├── sample_info.txt
├── contrasts.txt
├── 01_prepare_reference.sh
├── 02_prepare_ncRNA.sh
├── 03_build_indexes.sh
├── 04_fastp.sh
├── 05_remove_ncRNA.sh
├── 06_clean_reads_QC.sh
├── 07_STAR_mapping.sh
├── 08_STAR_mapping_QC.sh
├── 09_riboWaltz_RPF_QC.sh
├── 10_featureCounts_CDS.sh
├── 11_build_RPF_expression_matrix.sh
├── 12_DESeq2_RPF.sh
├── 13_merge_RPF_DE.sh
├── 14_prepare_RNA_Ribo_integration.sh
├── 15_Ribo_RNA_correlation.sh
├── 16_Ribo_RNA_DE_compare.sh
├── 17_translation_efficiency.sh
├── 18_differential_translation.sh
├── run_core_upstream.sh
├── run_integration.sh
├── scripts/
├── ref/
└── software/
```

## Main workflow

```text
FASTQ
 ↓
fastp
 ↓
Bowtie ncRNA depletion
 ↓
STAR genome alignment
 ↓
riboWaltz
RPF length / P-site / periodicity QC
 ↓
featureCounts on CDS
 ↓
RPF count + TPM matrices
 ↓
DESeq2 RPF differential abundance
 ↓
RNA-seq + Ribo-seq integration
 ├── correlation
 ├── RNA/RPF DE comparison
 ├── TE ratio
 └── differential translation (DTEG)
```

## Reference-workflow defaults

```text
fastp:
  adapter = TGGAATTCTCGGGTGCCAA
  retained length = 20–40 nt

ncRNA depletion:
  Bowtie -v 2

STAR:
  --outFilterMultimapNmax 1
  --outFilterMismatchNmax 2
  --quantMode TranscriptomeSAM GeneCounts

riboWaltz:
  custom RPF length = 26–38 nt

featureCounts:
  feature = CDS
  strand = 1
  paired-end = FALSE
```

These settings are inherited from the Arabidopsis reference and are deliberately
kept configurable in `00_config.sh`.

## Run

Core Ribo-seq processing:

```bash
bash run_core_upstream.sh
```

RNA/Ribo integration:

```bash
export RNA_TPM=../01_RNAseq/quanti/genes.TMM.TPM.matrix
export RNA_COUNTS=../01_RNAseq/quanti/genes.counts.matrix
export RNA_DE=../01_RNAseq/DE_analysis/de_result.tsv

bash run_integration.sh
```
