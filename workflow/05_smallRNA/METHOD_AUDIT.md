# small RNA / miRNA method–code audit

## 1. Study design

The source code contains ten libraries:

```text
PN4  rep1/rep2
PN7  rep1/rep2
PN10 rep1/rep2
PN15 rep1/rep2
PN18 rep1/rep2
```

All are single-end small-RNA libraries.

## 2. FastQC copy/paste typos in the historical text

The source FastQC block uses the rep1 FASTQ again for:

```text
PN15_rep2
PN18_rep2
```

while writing rep2 log names.

The organized `04_fastqc_trim.sh` loops over the sample table and therefore
uses each sample's own FASTQ. This is an intentional correction of an obvious
copy/paste typo and should be checked against the actual historical FastQC
outputs if those reports are cited.

## 3. Final length-filtering command

The source first uses:

```text
--length 18
```

without an upper-length threshold, then explicitly notes that reads >30 bp were
not removed and re-runs Trim Galore with:

```text
--length 18
--max_length 30
--quality 30
```

The organized workflow uses this corrected final block.

## 4. miRBase and Rfam versions

The source downloads miRBase `hairpin.fa` and `mature.fa`, but does not record
an exact miRBase release number.

The Rfam section states that the CURRENT database was re-downloaded around
September 2024, but `CURRENT` is not a reproducible release identifier.

Before publication release, record exact miRBase and Rfam versions/releases.

## 5. Bowtie versus Bowtie2 index inconsistency

The historical reference-preparation code uses:

```text
bowtie2-build repeat.fa repeat
bowtie2-build intron.fa intron
bowtie2-build exon.fa exon
```

but the downstream classification commands use the Bowtie 1 executable:

```text
bowtie ... ref/repeat
bowtie ... ref/intron
bowtie ... ref/exon
```

Bowtie1 and Bowtie2 use different index formats. This inconsistency must be
resolved from the actual files used on the server before `v1.0.0`.

For that reason, `03_prepare_reference.sh` preserves the historical build
commands, while `08_map_annotation_databases.sh` preserves the historical
mapping commands. Do not assume the pair is directly runnable until verified.

## 6. Exon/CDS reference is uncertain in the source

The source itself comments:

```text
"这里的外显子序列是cds不知道是否可行？"
```

Therefore the exact `ref/exon.fa` construction was not finalized in the supplied
text. The organized workflow requires an existing `ref/exon.fa` rather than
inventing how it was generated.

## 7. Novel-miRNA candidate definition

Novel-miRNA input reads are selected from annotations classified as:

```text
unknown
intron_sense
intron_antisense
```

This rule is preserved.

## 8. miRDP2 workflow

The historical path contains both `miRDP2-v1.1.5` and `1.1.4`, but the executed
file is explicitly named:

```text
miRDP2-v1.1.4_pipeline.bash
```

The repository therefore labels the executed pipeline as v1.1.4 unless later
server metadata shows otherwise.

## 9. Annotation pie-chart code was marked problematic

The source labels the pie-chart block:

```text
存在问题：需要重新进行修改
```

The organized core workflow does not automatically execute that plotting step.
The helper is still listed in `scripts/README.md` so the historical code can be
preserved separately.

## 10. Differential-expression block was also marked for later modification

The source uses DESeq2 through `script/run_DE_analysis.pl`, followed by
`process_files.sh`, and then comments that this section still required later
modification.

`18_DESeq2_miRNA.sh` preserves the historical logic but should be compared with
the final differential-miRNA result files used in the manuscript.

## 11. Target-prediction intersection is gene-level in the supplied code

The source parses psRobot, TargetFinder and psRNATarget predictions, maps
transcript IDs to gene IDs, then keeps only the first column before intersecting
the three methods.

Therefore the final `target_mirna.out` in this historical block represents genes
supported by all three prediction methods; it does not by itself guarantee that
the same miRNA ID was supported by all three tools.

Any later same-miRNA concordance filtering should be stored as a separate
downstream script rather than silently changing this historical step.

## 12. External online analysis

psRNATarget is an online platform. The final repository cannot reproduce the
web submission from command line unless the exact input/settings/output file
are archived.

At minimum preserve:

```text
input mature-miRNA FASTA
target transcript FASTA
psRNATarget parameter settings
downloaded psRNATarget output
```

or document them in the manuscript/repository.

## 13. Local helper scripts still required

See `scripts/README.md`. They include the custom annotation, statistics,
plotting, novel-miRNA renaming, DE, ID-mapping and target-prediction wrappers.

Do not reconstruct these helpers from memory for the final publication release;
copy the exact versions that generated the manuscript results.
