#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(ggplot2)
  library(ggrepel)
  library(ggvenn)
})

# ============================================================
# 1. Parameters
# ============================================================

option_list <- list(
  make_option(
    c("--ribo_de"),
    type = "character",
    help = "Ribo-seq differential expression result"
  ),
  make_option(
    c("--rna_de"),
    type = "character",
    help = "RNA-seq differential expression result"
  ),
  make_option(
    c("-o", "--outdir"),
    type = "character",
    default = "DE_compare_out",
    help = "Output directory [default: %default]"
  ),
  make_option(
    c("--padj_cutoff"),
    type = "double",
    default = 0.05,
    help = "Adjusted P value cutoff [default: %default]"
  ),
  make_option(
    c("--log2fc_cutoff"),
    type = "double",
    default = 1,
    help = "Absolute log2FC cutoff [default: %default]"
  ),
  make_option(
    c("--label_top"),
    type = "integer",
    default = 20,
    help = "Number of top genes to label in scatter plot [default: %default]"
  ),
  make_option(
    c("--scatter_width"),
    type = "double",
    default = 7,
    help = "Scatter plot PDF width [default: %default]"
  ),
  make_option(
    c("--scatter_height"),
    type = "double",
    default = 6,
    help = "Scatter plot PDF height [default: %default]"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$ribo_de) || is.null(opt$rna_de)) {
  stop("Please provide --ribo_de and --rna_de")
}

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)


# ============================================================
# 2. Read DE result
#    Input file must contain GeneID column
# ============================================================

read_de_result <- function(file) {
  
  required_cols <- c(
    "GeneID",
    "sampleA", "sampleB",
    "baseMeanA", "baseMeanB", "baseMean",
    "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"
  )
  
  dat <- read.delim(
    file,
    header = TRUE,
    sep = "\t",
    check.names = FALSE,
    stringsAsFactors = FALSE,
    comment.char = "",
    quote = ""
  )
  
  missing_cols <- setdiff(required_cols, colnames(dat))
  
  if (length(missing_cols) > 0) {
    stop(
      "Missing columns in ",
      file,
      ": ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  dat <- dat %>%
    mutate(
      GeneID = as.character(GeneID),
      sampleA = as.character(sampleA),
      sampleB = as.character(sampleB),
      baseMeanA = as.numeric(baseMeanA),
      baseMeanB = as.numeric(baseMeanB),
      baseMean = as.numeric(baseMean),
      log2FoldChange = as.numeric(log2FoldChange),
      lfcSE = as.numeric(lfcSE),
      stat = as.numeric(stat),
      pvalue = as.numeric(pvalue),
      padj = as.numeric(padj)
    )
  
  return(dat)
}


# ============================================================
# 3. Add significance and direction
# ============================================================

add_sig_info <- function(dat, prefix, padj_cutoff, log2fc_cutoff) {
  
  dat <- dat %>%
    mutate(
      sig = !is.na(padj) &
        !is.na(log2FoldChange) &
        padj < padj_cutoff &
        abs(log2FoldChange) >= log2fc_cutoff,
      
      direction = case_when(
        sig & log2FoldChange >= log2fc_cutoff ~ "Up",
        sig & log2FoldChange <= -log2fc_cutoff ~ "Down",
        TRUE ~ "NS"
      )
    )
  
  dat <- dat %>%
    rename_with(
      .fn = ~ paste0(prefix, "_", .x),
      .cols = -c(GeneID, sampleA, sampleB)
    )
  
  return(dat)
}


# ============================================================
# 4. Safe file name
# ============================================================

safe_name <- function(x) {
  gsub("[^A-Za-z0-9_.-]+", "_", x)
}


# ============================================================
# 5. Venn plot
#    Use ggvenn and fixed colors
# ============================================================

plot_venn <- function(set1, set2, name1, name2, title, outfile) {
  
  set1 <- unique(na.omit(set1))
  set2 <- unique(na.omit(set2))
  
  color1 <- "#4DBBD5"
  color2 <- "#E64B35"
  
  venn_list <- list()
  venn_list[[name1]] <- set1
  venn_list[[name2]] <- set2
  
  if (length(unique(c(set1, set2))) == 0) {
    
    p <- ggplot() +
      annotate(
        "text",
        x = 0,
        y = 0,
        label = paste0(title, "\nNo genes"),
        size = 6
      ) +
      theme_void() +
      labs(title = title)
    
  } else {
    
    p <- ggvenn(
      venn_list,
      fill_color = c(color1, color2),
      fill_alpha = 0.45,
      stroke_color = "black",
      stroke_size = 0.7,
      set_name_size = 5,
      text_size = 5,
      show_percentage = FALSE
    ) +
      theme_void(base_size = 14) +
      theme(
        plot.title = element_text(
          hjust = 0.5,
          face = "bold",
          size = 16
        )
      ) +
      labs(title = title)
  }
  
  ggsave(
    filename = outfile,
    plot = p,
    width = 6.5,
    height = 5.5
  )
}


# ============================================================
# 6. Read and preprocess data
# ============================================================

ribo_de <- read_de_result(opt$ribo_de) %>%
  add_sig_info(
    prefix = "ribo",
    padj_cutoff = opt$padj_cutoff,
    log2fc_cutoff = opt$log2fc_cutoff
  )

rna_de <- read_de_result(opt$rna_de) %>%
  add_sig_info(
    prefix = "rna",
    padj_cutoff = opt$padj_cutoff,
    log2fc_cutoff = opt$log2fc_cutoff
  )

rna_comparisons <- rna_de %>%
  distinct(sampleA, sampleB) %>%
  mutate(in_rna = TRUE)

ribo_comparisons <- ribo_de %>%
  distinct(sampleA, sampleB) %>%
  mutate(in_ribo = TRUE)

comparison_status <- full_join(
  rna_comparisons,
  ribo_comparisons,
  by = c("sampleA", "sampleB")
) %>%
  mutate(
    in_rna = replace_na(in_rna, FALSE),
    in_ribo = replace_na(in_ribo, FALSE),
    comparison_status = case_when(
      in_rna & in_ribo ~ "Both",
      in_rna & !in_ribo ~ "RNA_only",
      !in_rna & in_ribo ~ "Ribo_only",
      TRUE ~ "Unknown"
    )
  )

write.table(
  comparison_status,
  file = file.path(opt$outdir, "comparison_status.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

valid_comparisons <- comparison_status %>%
  filter(comparison_status == "Both") %>%
  select(sampleA, sampleB)

if (nrow(valid_comparisons) == 0) {
  stop("No common sampleA/sampleB comparison found between RNA-seq and Ribo-seq results.")
}

rna_de <- rna_de %>%
  semi_join(valid_comparisons, by = c("sampleA", "sampleB"))

ribo_de <- ribo_de %>%
  semi_join(valid_comparisons, by = c("sampleA", "sampleB"))


# ============================================================
# 7. Merge RNA-seq and Ribo-seq DE results
#    Use full_join to keep RNA-only and Ribo-only genes
# ============================================================

merged_de <- full_join(
  rna_de,
  ribo_de,
  by = c("GeneID", "sampleA", "sampleB")
)

merged_de <- merged_de %>%
  mutate(
    rna_present = !is.na(rna_log2FoldChange) & !is.na(rna_padj),
    ribo_present = !is.na(ribo_log2FoldChange) & !is.na(ribo_padj),
    
    data_source = case_when(
      rna_present & ribo_present ~ "Both",
      rna_present & !ribo_present ~ "RNA_only",
      !rna_present & ribo_present ~ "Ribo_only",
      TRUE ~ "Unknown"
    ),
    
    rna_sig = replace_na(rna_sig, FALSE),
    ribo_sig = replace_na(ribo_sig, FALSE),
    rna_direction = replace_na(rna_direction, "NS"),
    ribo_direction = replace_na(ribo_direction, "NS"),
    
    scatter_class = case_when(
      rna_sig & !ribo_sig ~ "RNA_sig_only",
      !rna_sig & ribo_sig ~ "Ribo_sig_only",
      
      rna_sig & ribo_sig &
        sign(rna_log2FoldChange) == sign(ribo_log2FoldChange) ~ "Homodirection",
      
      rna_sig & ribo_sig &
        sign(rna_log2FoldChange) != sign(ribo_log2FoldChange) ~ "Opposite",
      
      TRUE ~ "Not_significant"
    )
  )

write.table(
  merged_de,
  file = file.path(opt$outdir, "merged_DE_result.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ============================================================
# 8. Plot and summarize by comparison
# ============================================================

scatter_class_levels <- c(
  "RNA_sig_only",
  "Ribo_sig_only",
  "Homodirection",
  "Opposite",
  "Not_significant"
)

scatter_colors <- c(
  "RNA_sig_only" = "#377EB8",
  "Ribo_sig_only" = "#E69F00",
  "Homodirection" = "#E41A1C",
  "Opposite" = "#984EA3",
  "Not_significant" = "grey70"
)

comparison_list <- merged_de %>%
  distinct(sampleA, sampleB) %>%
  drop_na()

all_class_count <- list()
all_de_count <- list()
all_source_count <- list()

for (i in seq_len(nrow(comparison_list))) {
  
  sample_a <- comparison_list$sampleA[i]
  sample_b <- comparison_list$sampleB[i]
  
  comparison_name <- paste0(sample_a, "_vs_", sample_b)
  comparison_safe <- safe_name(comparison_name)
  
  message("Processing: ", comparison_name)
  
  dat <- merged_de %>%
    filter(sampleA == sample_a, sampleB == sample_b) %>%
    mutate(
      scatter_class = factor(
        scatter_class,
        levels = scatter_class_levels
      )
    )
  
  # ------------------------------------------------------------
  # 8.1 Scatter classification count
  # ------------------------------------------------------------
  
  class_count <- dat %>%
    count(
      sampleA,
      sampleB,
      scatter_class,
      name = "Gene_number"
    ) %>%
    complete(
      sampleA = sample_a,
      sampleB = sample_b,
      scatter_class = factor(scatter_class_levels, levels = scatter_class_levels),
      fill = list(Gene_number = 0)
    )
  
  all_class_count[[comparison_name]] <- class_count
  
  
  # ------------------------------------------------------------
  # 8.2 Data source count
  # ------------------------------------------------------------
  
  source_count <- dat %>%
    count(
      sampleA,
      sampleB,
      data_source,
      name = "Gene_number"
    )
  
  all_source_count[[comparison_name]] <- source_count
  
  
  # ------------------------------------------------------------
  # 8.3 Scatter plot
  #     Only genes with both RNA and Ribo log2FC are plotted
  # ------------------------------------------------------------
  
  dat_plot <- dat %>%
    filter(
      !is.na(rna_log2FoldChange),
      !is.na(ribo_log2FoldChange)
    )
  
  dat_ns <- dat_plot %>%
    filter(scatter_class == "Not_significant")
  
  dat_sig <- dat_plot %>%
    filter(scatter_class != "Not_significant")
  
  legend_dummy <- tibble(
    rna_log2FoldChange = 0,
    ribo_log2FoldChange = 0,
    scatter_class = factor(
      scatter_class_levels,
      levels = scatter_class_levels
    )
  )
  
  label_dat <- dat_sig %>%
    filter(
      !is.na(rna_log2FoldChange),
      !is.na(ribo_log2FoldChange)
    ) %>%
    mutate(
      label_score = abs(rna_log2FoldChange) + abs(ribo_log2FoldChange)
    ) %>%
    arrange(desc(label_score)) %>%
    slice_head(n = opt$label_top)
  
  p_scatter <- ggplot() +
    geom_hline(
      yintercept = c(-opt$log2fc_cutoff, opt$log2fc_cutoff),
      linetype = "dashed",
      color = "grey70"
    ) +
    geom_vline(
      xintercept = c(-opt$log2fc_cutoff, opt$log2fc_cutoff),
      linetype = "dashed",
      color = "grey70"
    ) +
    geom_point(
      data = legend_dummy,
      aes(
        x = rna_log2FoldChange,
        y = ribo_log2FoldChange,
        color = scatter_class
      ),
      size = 1.6,
      alpha = 0,
      show.legend = TRUE
    ) +
    geom_point(
      data = dat_ns,
      aes(
        x = rna_log2FoldChange,
        y = ribo_log2FoldChange,
        color = scatter_class
      ),
      size = 1.2,
      alpha = 0.35,
      na.rm = TRUE
    ) +
    geom_point(
      data = dat_sig,
      aes(
        x = rna_log2FoldChange,
        y = ribo_log2FoldChange,
        color = scatter_class
      ),
      size = 1.6,
      alpha = 0.85,
      na.rm = TRUE
    ) +
    geom_text_repel(
      data = label_dat,
      aes(
        x = rna_log2FoldChange,
        y = ribo_log2FoldChange,
        label = GeneID
      ),
      size = 3,
      max.overlaps = Inf,
      show.legend = FALSE
    ) +
    scale_color_manual(
      values = scatter_colors,
      breaks = scatter_class_levels,
      limits = scatter_class_levels,
      drop = FALSE
    ) +
    guides(
      color = guide_legend(
        override.aes = list(
          alpha = 1,
          size = 3
        )
      )
    ) +
    theme_bw(base_size = 13) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold"
      )
    ) +
    labs(
      title =  comparison_name,
      x = "RNA-seq log2FoldChange",
      y = "Ribo-seq log2FoldChange",
      color = "Category"
    )
  
  ggsave(
    filename = file.path(
      opt$outdir,
      paste0(comparison_safe, "_log2FC_scatter.pdf")
    ),
    plot = p_scatter,
    width = opt$scatter_width,
    height = opt$scatter_height
  )
  
  
  # ------------------------------------------------------------
  # 8.4 Extract upregulated and downregulated genes
  #     Use dat, not dat_plot, to keep RNA-only and Ribo-only genes
  # ------------------------------------------------------------
  
  rna_up <- dat %>%
    filter(rna_direction == "Up") %>%
    pull(GeneID) %>%
    unique()
  
  ribo_up <- dat %>%
    filter(ribo_direction == "Up") %>%
    pull(GeneID) %>%
    unique()
  
  rna_down <- dat %>%
    filter(rna_direction == "Down") %>%
    pull(GeneID) %>%
    unique()
  
  ribo_down <- dat %>%
    filter(ribo_direction == "Down") %>%
    pull(GeneID) %>%
    unique()
  
  
  # ------------------------------------------------------------
  # 8.5 Venn plots
  # ------------------------------------------------------------
  
  plot_venn(
    set1 = rna_up,
    set2 = ribo_up,
    name1 = "RNA Up",
    name2 = "Ribo Up",
    title = paste0(comparison_name, " upregulated genes"),
    outfile = file.path(
      opt$outdir,
      paste0(comparison_safe, "_up_venn.pdf")
    )
  )
  
  plot_venn(
    set1 = rna_down,
    set2 = ribo_down,
    name1 = "RNA Down",
    name2 = "Ribo Down",
    title = paste0(comparison_name, " downregulated genes"),
    outfile = file.path(
      opt$outdir,
      paste0(comparison_safe, "_down_venn.pdf")
    )
  )
  
  
  # ------------------------------------------------------------
  # 8.6 Up/down gene count
  # ------------------------------------------------------------
  
  de_count <- tibble(
    sampleA = sample_a,
    sampleB = sample_b,
    Omics = c("RNA-seq", "RNA-seq", "Ribo-seq", "Ribo-seq"),
    Direction = c("Up", "Down", "Up", "Down"),
    Gene_number = c(
      length(rna_up),
      length(rna_down),
      length(ribo_up),
      length(ribo_down)
    )
  )
  
  all_de_count[[comparison_name]] <- de_count
  
  
  # ------------------------------------------------------------
  # 8.7 Up/down gene count barplot
  # ------------------------------------------------------------
  
  p_bar <- ggplot(
    de_count,
    aes(
      x = Omics,
      y = Gene_number,
      fill = Direction
    )
  ) +
    geom_col(
      position = position_dodge(width = 0.75),
      width = 0.65
    ) +
    geom_text(
      aes(label = Gene_number),
      position = position_dodge(width = 0.75),
      vjust = -0.3,
      size = 4
    ) +
    scale_fill_manual(
      values = c(
        "Up" = "#E41A1C",
        "Down" = "#377EB8"
      ),
      breaks = c("Up", "Down"),
      limits = c("Up", "Down"),
      drop = FALSE
    ) +
    theme_bw(base_size = 13) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold"
      )
    ) +
    labs(
      title = comparison_name,
      x = NULL,
      y = "Gene number",
      fill = "Direction"
    ) +
    expand_limits(
      y = max(de_count$Gene_number) * 1.15 + 1
    )
  
  ggsave(
    filename = file.path(
      opt$outdir,
      paste0(comparison_safe, "_DE_count_barplot.pdf")
    ),
    plot = p_bar,
    width = 5.5,
    height = 5
  )
}


# ============================================================
# 9. Write summary tables
# ============================================================

bind_rows(all_class_count) %>%
  write.table(
    file = file.path(opt$outdir, "scatter_classification_count.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

bind_rows(all_de_count) %>%
  write.table(
    file = file.path(opt$outdir, "DE_gene_count.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

bind_rows(all_source_count) %>%
  write.table(
    file = file.path(opt$outdir, "data_source_count.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

message("Done. Results were saved to: ", opt$outdir)