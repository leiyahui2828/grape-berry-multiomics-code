# EM-seq method–code audit

1. The uploaded code contains all 15 libraries: five stages × three replicates.

2. The historical workflow uses **Bismark for alignment** and
   **BatMeth2 calmeth for methylation calling**. Methods should describe both
   steps if these files generated the manuscript results.

3. Duplicate handling is performed in the BatMeth2 calling step via
   `--remove_dup`; the supplied BAM-preparation section does not show Picard
   duplicate removal.

4. A minimum methylation coverage threshold of 4 is directly supported by both
   `calmeth` and `methyGff`.

5. Two separate window/step settings occur:
   - BatMeth2 calmeth: `--Regions 600 --step 50000`
   - DMR analysis: `--window 600 --step 200`
   They are different analysis stages and should not be conflated.

6. One isolated exploratory command involving `PN4_rep1.PN19.bam` was excluded
   because it is not part of the five-stage experiment.

7. The historical output directory is spelled `site_samll_window_ML`.
   The organized code retains this spelling for compatibility.

8. The code includes gene/promoter and TE methylation profiling. The final
   repository should document the exact PN40024 gene annotation and grape TE
   annotation versions.

9. Differential methylation requires the exact historical:
   `software/run_DM_analysis.R`, `samples_info.txt`, and `contrasts.txt`.

10. The large Singularity image `sifs/methtools20220515.sif` should not be
    committed to GitHub. Record the contained Bismark/BatMeth2 versions instead.
