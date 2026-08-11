# Historical helper scripts to restore

The uploaded analysis text calls the following local/custom helper scripts.
Copy the exact historical versions into this directory before publication
release.

```text
get_intron.pl
stat_srna.pl
draw_length_srna.R
draw_length_srna1.R
common_specific_reads.pl
stat_bwt.pl
srna_anno.pl
draw_srna_annPie1.R
novel_name_mirdp2.pl
mirna_base_content.pl
draw_base_content1.R
blastRfam2family.pl
gene_family/get_meme_pos.pl
draw_motif_seqlogo.R
expression_draw1.R
run_DE_analysis.pl
process_files.sh
map_data_ids
sum_seq.sh
```

The source also calls:

```text
targetfinder_threads.pl
combine_fastp_report.R
```

If these were local wrappers used for the manuscript analysis, archive the exact
versions as well.

Do not recreate these scripts from memory if their outputs were used in the
paper. Copy the versions that actually generated the manuscript results.
