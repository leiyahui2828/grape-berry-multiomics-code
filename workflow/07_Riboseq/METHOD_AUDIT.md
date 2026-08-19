# Ribo-seq method–code audit

## Source status

The uploaded archive is an Arabidopsis workflow. The files here adapt its
structure to the grape five-stage / three-replicate design; they are not yet
claimed to be the exact executed grape workflow.

## Items requiring grape-specific verification

1. **Single-end libraries**  
   The reference pipeline is single-end and uses `RPF_IS_PAIRED=FALSE`.
   Confirm against the actual grape FASTQ files.

2. **Adapter**  
   The reference uses `TGGAATTCTCGGGTGCCAA`. Confirm against the grape
   Ribo-seq library-preparation kit.

3. **Length filters**  
   The reference keeps 20–40 nt after fastp but uses 26–38 nt in riboWaltz.
   These represent different analysis stages. Final manuscript wording should
   follow the actual grape data and final analysis.

4. **ncRNA depletion**  
   Do not reuse Arabidopsis ncRNA sequences. Supply a grape-specific rRNA/tRNA
   (and, if used, snRNA/snoRNA/miRNA) FASTA and document its source/version.

5. **STAR statistics filename**  
   The reference script searches for `*.final.out`; STAR produces
   `Log.final.out`. The organized code uses `Log.final.out`.

6. **riboWaltz versus RiboCode**  
   The uploaded reference explicitly uses `run_riboWaltz.R`. If an older
   vendor/method draft says RiboCode, the final repository should match the
   workflow that actually generated the grape manuscript results.

7. **featureCounts parameters**  
   Reference values are `CDS`, strand=1, single-end. Confirm grape library
   strandness before release.

8. **RPF DE vs differential translation**  
   DESeq2 on Ribo-seq counts alone tests differential RPF abundance. It is not
   differential translation efficiency. The latter is handled separately by
   RNA+Ribo integration with `DTEG_modified.R`.

9. **TE ratio vs DTEG**  
   `ribo_rna_TE.R` computes a descriptive normalized Ribo/RNA ratio.
   `DTEG_modified.R` uses RNA and Ribo count matrices for differential
   translational regulation. These should not be described as the same test.

10. **RNA/Ribo pairing**  
    Confirm identical biological sample IDs, gene IDs and annotation release
    for all 15 RNA-seq and Ribo-seq samples.

11. **Contrasts**  
    `contrasts.txt` currently includes all ten pairwise stage comparisons with
    later stage vs earlier stage. Remove unused comparisons if the final paper
    used a smaller set.

12. **Helper scripts**  
    The uploaded archive references many helper scripts outside its own
    directory. Their names are listed under `scripts/README.md`. Exact grape
    versions must be restored or validated before `v1.0.0`.
