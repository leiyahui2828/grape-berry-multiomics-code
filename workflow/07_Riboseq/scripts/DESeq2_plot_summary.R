#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(ggplot2)
  library(EnhancedVolcano)
})

# ============================================================
# 1. Parameters
# ============================================================

option_list <- list(
  make_option(
    c("-i", "--input"),
    type = "character",
    help = "DESeq2 result table"
  ),
  make_option(
    c("-o", "--outdir"),
    type = "character",
    default = "DESeq2_plot_summary",
    help = "Output directory [default: %default]"
  ),
  make_option(
    c("--prefix"),
    type = "character",
    default = "DESeq2",
    help = "Output file prefix [default: %default]"
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
    help = "Absolute log2 fold change cutoff [default: %default]"
  ),
  make_option(
    c("--label_top"),
    type = "integer",
    default = 20,
    help = "Number of top significant genes to label in volcano plot [default: %default]"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$input)) {
  stop("Please provide --input")
}

if (opt$padj_cutoff <= 0 || opt$padj_cutoff >= 1) {
  stop("--padj_cutoff should be between 0 and 1")
}

if (opt$log2fc_cutoff < 0) {
  stop("--log2fc_cutoff should be non-negative")
}

if (opt$label_top < 0) {
  stop("--label_top should be non-negative")
}

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)


# ============================================================
# 2. Basic function
# ============================================================

safe_name <- function(x) {
  gsub("[^A-Za-z0-9_.-]+", "_", x)
}

prefix_safe <- safe_name(opt$prefix)


# ============================================================
# 3. Read DESeq2 result
# ============================================================

de_result <- read.delim(
  opt$input,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  comment.char = "",
  quote = ""
)

if (ncol(de_result) < 5) {
  stop("Input table should contain at least GeneID, baseMean, log2FoldChange, pvalue and padj columns.")
}

if (!"GeneID" %in% colnames(de_result)) {
  colnames(de_result)[1] <- "GeneID"
}

required_cols <- c(
  "GeneID",
  "baseMean",
  "log2FoldChange",
  "pvalue",
  "padj"
)

missing_cols <- setdiff(required_cols, colnames(de_result))

if (length(missing_cols) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

de_result <- de_result %>%
  mutate(
    GeneID = as.character(GeneID),
    baseMean = as.numeric(baseMean),
    log2FoldChange = as.numeric(log2FoldChange),
    pvalue = as.numeric(pvalue),
    padj = as.numeric(padj)
  )


# ============================================================
# 4. Annotate significant genes
# ============================================================

de_annotated <- de_result %>%
  mutate(
    minus_log10_padj = -log10(padj),
    minus_log10_padj = ifelse(
      is.finite(minus_log10_padj),
      minus_log10_padj,
      NA
    ),
    
    DEG_class = case_when(
      !is.na(padj) &
        padj < opt$padj_cutoff &
        !is.na(log2FoldChange) &
        log2FoldChange >= opt$log2fc_cutoff ~ "Up",
      
      !is.na(padj) &
        padj < opt$padj_cutoff &
        !is.na(log2FoldChange) &
        log2FoldChange <= -opt$log2fc_cutoff ~ "Down",
      
      TRUE ~ "Not_sig"
    ),
    
    DEG_class = factor(
      DEG_class,
      levels = c("Down", "Not_sig", "Up")
    )
  )

write.table(
  de_annotated,
  file = file.path(opt$outdir, paste0(prefix_safe, "_annotated.tsv")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ============================================================
# 5. Summary table
# ============================================================

summary_table <- de_annotated %>%
  summarise(
    total_genes = n(),
    genes_with_nonzero_baseMean = sum(baseMean > 0, na.rm = TRUE),
    genes_with_pvalue = sum(!is.na(pvalue)),
    genes_with_padj = sum(!is.na(padj)),
    up_genes = sum(DEG_class == "Up", na.rm = TRUE),
    down_genes = sum(DEG_class == "Down", na.rm = TRUE),
    not_sig_genes = sum(DEG_class == "Not_sig", na.rm = TRUE),
    padj_cutoff = opt$padj_cutoff,
    log2fc_cutoff = opt$log2fc_cutoff
  )

write.table(
  summary_table,
  file = file.path(opt$outdir, paste0(prefix_safe, "_summary.tsv")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ============================================================
# 6. Volcano plot
# ============================================================

volcano_dat <- de_annotated %>%
  filter(
    !is.na(log2FoldChange),
    !is.na(padj)
  )

volcano_labels <- rep("", nrow(volcano_dat))

if (opt$label_top > 0) {
  label_genes <- volcano_dat %>%
    filter(DEG_class != "Not_sig") %>%
    arrange(padj, desc(abs(log2FoldChange))) %>%
    slice_head(n = opt$label_top) %>%
    pull(GeneID)
  
  volcano_labels[volcano_dat$GeneID %in% label_genes] <- volcano_dat$GeneID[
    volcano_dat$GeneID %in% label_genes
  ]
}

volcano_colors <- ifelse(
  volcano_dat$DEG_class == "Up",
  "#E41A1C",
  ifelse(
    volcano_dat$DEG_class == "Down",
    "#377EB8",
    "grey70"
  )
)

names(volcano_colors) <- as.character(volcano_dat$DEG_class)

p_volcano <- EnhancedVolcano(
  volcano_dat,
  lab = volcano_labels,
  x = "log2FoldChange",
  y = "padj",
  title = paste0(prefix_safe, " volcano plot"),
  subtitle = NULL,
  caption = NULL,
  xlab = bquote(~Log[2]~ "fold change"),
  ylab = bquote(~-Log[10]~ "adjusted P value"),
  pCutoff = opt$padj_cutoff,
  FCcutoff = opt$log2fc_cutoff,
  pointSize = 1.5,
  labSize = 3,
  colCustom = volcano_colors,
  colAlpha = 0.8,
  legendPosition = "none",
  drawConnectors = TRUE,
  widthConnectors = 0.3,
  max.overlaps = Inf
)

ggsave(
  filename = file.path(opt$outdir, paste0(prefix_safe, "_volcano.pdf")),
  plot = p_volcano,
  width = 7,
  height = 6
)


# ============================================================
# 7. MA plot
# ============================================================

p_ma <- ggplot(
  de_annotated,
  aes(
    x = log10(baseMean + 1),
    y = log2FoldChange,
    color = DEG_class
  )
) +
  geom_hline(
    yintercept = 0,
    color = "grey50"
  ) +
  geom_hline(
    yintercept = c(-opt$log2fc_cutoff, opt$log2fc_cutoff),
    linetype = "dashed",
    color = "grey70"
  ) +
  geom_point(
    size = 1.2,
    alpha = 0.75,
    na.rm = TRUE
  ) +
  scale_color_manual(
    values = c(
      "Down" = "#377EB8",
      "Not_sig" = "grey70",
      "Up" = "#E41A1C"
    ),
    limits = c("Down", "Not_sig", "Up"),
    drop = FALSE
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    legend.position = "none"
  ) +
  labs(
    title = paste0(prefix_safe, " MA plot"),
    x = "log10(baseMean + 1)",
    y = "log2FoldChange"
  )

ggsave(
  filename = file.path(opt$outdir, paste0(prefix_safe, "_MAplot.pdf")),
  plot = p_ma,
  width = 7,
  height = 5.5
)


# ============================================================
# 8. Count barplot
# ============================================================

count_dat <- de_annotated %>%
  count(DEG_class, name = "Gene_number") %>%
  complete(
    DEG_class = factor(
      c("Down", "Not_sig", "Up"),
      levels = c("Down", "Not_sig", "Up")
    ),
    fill = list(Gene_number = 0)
  )

p_count <- ggplot(
  count_dat,
  aes(
    x = DEG_class,
    y = Gene_number,
    fill = DEG_class
  )
) +
  geom_col(width = 0.65) +
  geom_text(
    aes(label = Gene_number),
    vjust = -0.3,
    size = 4
  ) +
  scale_fill_manual(
    values = c(
      "Down" = "#377EB8",
      "Not_sig" = "grey70",
      "Up" = "#E41A1C"
    ),
    limits = c("Down", "Not_sig", "Up"),
    drop = FALSE
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    legend.position = "none"
  ) +
  labs(
    title = paste0(prefix_safe, " gene count"),
    x = NULL,
    y = "Gene number"
  ) +
  expand_limits(
    y = max(count_dat$Gene_number) * 1.1 + 1
  )

ggsave(
  filename = file.path(opt$outdir, paste0(prefix_safe, "_count_barplot.pdf")),
  plot = p_count,
  width = 6,
  height = 5
)


message("Done. Results were saved to: ", opt$outdir)