# CUT&Tag method–code audit

This audit separates what is directly supported by the uploaded files from
what still requires confirmation.

## 1. Source coverage

The uploaded `cut_tag_code` dry-run contains H3K27ac only:

- PN4 rep1/rep2
- PN7 rep1/rep2
- PN10 rep1/rep2
- PN15 rep1/rep2
- PN18 rep1/rep2

The uploaded Snakefile is generic and supports both narrow and broad peaks,
but the exact config files for H3K4me3, H3K9ac and H3K27me3 were not included.

Therefore:

- H3K27ac command parameters in this directory are reconstructed directly from
  the supplied execution log.
- H3K4me3/H3K9ac/H3K27me3 are structurally prepared, including the correct
  narrowPeak/broadPeak filename logic, but their exact historical MACS3-extra
  parameters and replicate-selection config must still be checked.

## 2. Peak classes

The repository uses:

```text
H3K27ac  narrowPeak
H3K4me3  narrowPeak
H3K9ac   narrowPeak
H3K27me3 broadPeak
```

For H3K27me3, `09_MACS3_peak_calling.sh` adds `--broad` and expects
`*_peaks.broadPeak`.

For narrow marks, it expects `*_peaks.narrowPeak`.

`17_extract_peak_sequence.sh` only performs summit-centered extraction for
narrow peaks because the narrowPeak format contains the summit-offset column.

## 3. Blacklist issue that must be resolved

The H3K27ac dry-run explicitly used:

```text
--blackListFileName ref/hg38.blacklist.bed
```

The study organism is Vitis vinifera, so this file name is suspicious.

Before the final release:

```bash
head ref/hg38.blacklist.bed
cut -f1 ref/hg38.blacklist.bed | sort -u | head
```

If it is truly an hg38 human blacklist, it should not be used for grape.
If it is a grape blacklist that was simply given a misleading filename, rename
it and document its construction.

For safety, the publication-oriented scripts leave `BLACKLIST` empty by
default instead of silently applying the questionable file.

## 4. Duplicate handling differs from the current manuscript draft

The supplied H3K27ac code uses:

```text
alignmentSieve --ignoreDuplicates
```

The current manuscript draft has previously described duplicate removal with
Picard MarkDuplicates.

These are not the same workflow description. The final Methods must follow the
actual final analysis that generated the manuscript figures.

## 5. Peak caller format

The supplied H3K27ac run used:

```text
--format BAM
--nomodel --shift -100 --extsize 200
```

even though the libraries are paired-end.

This directory preserves that H3K27ac command because the goal here is to
organize the actual historical analysis, not silently replace it with a new
pipeline. If the final manuscript analysis is re-run using BAMPE, it should be
released as a new verified analysis version.

## 6. Two-stage peak thresholding

The H3K27ac source workflow first called MACS3 with:

```text
--qvalue 0.075
```

and then filtered the peak-file q-score column using:

```text
awk '$9>3'
```

Thus the final replicate-level cutoff files are more stringent than the initial
MACS3 call.

The H3K27ac IDR rule in the supplied Snakefile uses the raw MACS3 replicate peak
files as IDR input rather than the q-score-filtered replicate peak files.
This behavior is preserved in `11_select_reproducible_peaks.sh`.

## 7. Differential peak method

The supplied H3K27ac dry-run uses:

```text
run_DE_analysis.pl --method edgeR --dispersion 0.1
```

for all ten pairwise developmental-stage comparisons.

If the manuscript currently states that differential histone-mark peaks were
identified with DESeq2, that statement does not match this uploaded H3K27ac
workflow and should be reconciled with the final files actually used in the
figures.

## 8. Peak annotation method

The supplied workflow performs pipeline peak annotation with UROPA.

If ChIPseeker was later used for manuscript-level annotation/plots, keep that
as a separate downstream R script and distinguish it from the upstream UROPA
annotation.

## 9. No matched input/control in the H3K27ac dry-run

The MACS3 H3K27ac commands shown in the uploaded execution log do not contain a
`-c` matched control BAM.

This is consistent with an enrichment-based peak-calling workflow without
matched input/IgG controls.

## 10. Files still needed before release

Custom/local helpers referenced by the supplied workflow:

```text
scripts/combine_fastp_report.R
scripts/plot_macs_cutoff_analysis.py
software/RunFeatureCounts/run-featurecounts.R
software/RunFeatureCounts/abundance_estimates_to_matrix.pl
```

External executables/packages should be installed, not copied into this repo:

```text
fastp
Bowtie2
SAMtools
deepTools
MACS3
IDR
bedtools
UROPA
Trinity run_DE_analysis.pl
edgeR
```
