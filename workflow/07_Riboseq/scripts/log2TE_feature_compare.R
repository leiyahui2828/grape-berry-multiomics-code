#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(ggplot2)
  library(ggpubr)
})

# ============================================================
# 1. Parameters
# ============================================================

option_list <- list(
  make_option(
    c("--logte"),
    type = "character",
    help = "log2TE matrix/table"
  ),
  make_option(
    c("--feature"),
    type = "character",
    help = "CDS/UTR feature table with length/GC/NMFE columns"
  ),
  make_option(
    c("--te_col"),
    type = "character",
    default = "all",
    help = "Column name in logTE table. Use 'all' to calculate row mean [default: %default]"
  ),
  make_option(
    c("-o", "--outdir"),
    type = "character",
    default = "log2TE_feature_compare",
    help = "Output directory [default: %default]"
  ),
  make_option(
    c("--min_group_n"),
    type = "integer",
    default = 3,
    help = "Minimum number of genes per TE class for statistical test [default: %default]"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$logte) || is.null(opt$feature)) {
  stop("Please provide --logte and --feature")
}

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(opt$outdir, "boxplot"), recursive = TRUE, showWarnings = FALSE)


# ============================================================
# 2. Basic functions
# ============================================================

safe_name <- function(x) {
  gsub("[^A-Za-z0-9_.-]+", "_", x)
}

read_table_first_col_gene <- function(file) {
  
  dat <- read.delim(
    file,
    header = TRUE,
    sep = "\t",
    check.names = FALSE,
    stringsAsFactors = FALSE,
    comment.char = "",
    quote = ""
  )
  
  if (ncol(dat) < 2) {
    stop("Input table should contain at least two columns: ", file)
  }
  
  colnames(dat)[1] <- "GeneID"
  dat$GeneID <- as.character(dat$GeneID)
  
  if (anyDuplicated(dat$GeneID) > 0) {
    dup_gene <- dat$GeneID[duplicated(dat$GeneID)][1]
    stop("Duplicated GeneID found in ", file, ": ", dup_gene)
  }
  
  return(dat)
}


# ============================================================
# 3. Read log2TE table and generate logTE_used
# ============================================================

logte_dat <- read_table_first_col_gene(opt$logte)

if (opt$te_col != "all") {
  
  if (!opt$te_col %in% colnames(logte_dat)) {
    stop(
      "Column specified by --te_col not found in logTE table: ",
      opt$te_col
    )
  }
  
  logte_use <- logte_dat %>%
    transmute(
      GeneID,
      logTE_used = as.numeric(.data[[opt$te_col]])
    )
  
  te_source <- opt$te_col
  
} else {
  
  numeric_cols <- setdiff(colnames(logte_dat), "GeneID")
  
  logte_numeric <- logte_dat %>%
    select(all_of(numeric_cols)) %>%
    mutate(across(everything(), as.numeric))
  
  logte_use <- tibble(
    GeneID = logte_dat$GeneID,
    logTE_used = rowMeans(logte_numeric, na.rm = TRUE)
  )
  
  logte_use$logTE_used[is.nan(logte_use$logTE_used)] <- NA
  te_source <- "mean_logTE"
}


# ============================================================
# 4. Classify genes by log2TE
# ============================================================

te_class_levels <- c(
  "log2(TE)<=-1",
  "-1<log2(TE)<=0",
  "0<log2(TE)<1",
  "log2(TE)>=1"
)

logte_class <- logte_use %>%
  mutate(
    TE_class = case_when(
      logTE_used <= -1 ~ "log2(TE)<=-1",
      logTE_used > -1 & logTE_used <= 0 ~ "-1<log2(TE)<=0",
      logTE_used > 0 & logTE_used < 1 ~ "0<log2(TE)<1",
      logTE_used >= 1 ~ "log2(TE)>=1",
      TRUE ~ NA_character_
    ),
    TE_class = factor(
      TE_class,
      levels = te_class_levels
    )
  ) %>%
  filter(!is.na(logTE_used), !is.na(TE_class))

write.table(
  logte_class,
  file = file.path(opt$outdir, "logTE_classification.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ============================================================
# 5. Read feature table and merge
# ============================================================

feature_dat <- read_table_first_col_gene(opt$feature)

feature_dat <- feature_dat %>%
  mutate(
    across(
      .cols = -GeneID,
      .fns = as.numeric
    )
  )

feature_cols <- setdiff(colnames(feature_dat), "GeneID")

if (length(feature_cols) == 0) {
  stop("No feature columns found in feature table.")
}

merged_dat <- inner_join(
  logte_class,
  feature_dat,
  by = "GeneID"
)

if (nrow(merged_dat) == 0) {
  stop("No common GeneID found between logTE table and feature table.")
}

write.table(
  merged_dat,
  file = file.path(opt$outdir, "merged_logTE_feature.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ============================================================
# 6. Summary statistics by TE class
# ============================================================

feature_summary <- merged_dat %>%
  pivot_longer(
    cols = all_of(feature_cols),
    names_to = "feature",
    values_to = "value"
  ) %>%
  filter(!is.na(value)) %>%
  group_by(feature, TE_class) %>%
  summarise(
    n = n(),
    mean = mean(value, na.rm = TRUE),
    median = median(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    min = min(value, na.rm = TRUE),
    max = max(value, na.rm = TRUE),
    .groups = "drop"
  )

write.table(
  feature_summary,
  file = file.path(opt$outdir, "feature_summary_by_class.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ============================================================
# 7. Kruskal-Wallis test
# ============================================================

kruskal_result <- map_dfr(
  feature_cols,
  function(feat) {
    
    test_dat <- merged_dat %>%
      select(GeneID, TE_class, value = all_of(feat)) %>%
      filter(!is.na(value), !is.na(TE_class))
    
    class_n <- test_dat %>%
      count(TE_class)
    
    valid_class_number <- sum(class_n$n >= opt$min_group_n)
    
    if (valid_class_number < 2) {
      return(
        tibble(
          feature = feat,
          test = "Kruskal-Wallis",
          pvalue = NA_real_,
          valid_class_number = valid_class_number,
          note = "Less than two valid classes"
        )
      )
    }
    
    test_dat <- test_dat %>%
      group_by(TE_class) %>%
      filter(n() >= opt$min_group_n) %>%
      ungroup()
    
    kt <- kruskal.test(value ~ TE_class, data = test_dat)
    
    tibble(
      feature = feat,
      test = "Kruskal-Wallis",
      pvalue = kt$p.value,
      valid_class_number = valid_class_number,
      note = "OK"
    )
  }
)

write.table(
  kruskal_result,
  file = file.path(opt$outdir, "feature_kruskal_test.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ============================================================
# 8. Pairwise Wilcoxon rank-sum test without p-value correction
# ============================================================

pairwise_comparisons <- combn(
  te_class_levels,
  2,
  simplify = FALSE
)

pairwise_result <- map_dfr(
  feature_cols,
  function(feat) {
    
    test_dat <- merged_dat %>%
      select(GeneID, TE_class, value = all_of(feat)) %>%
      filter(!is.na(value), !is.na(TE_class))
    
    map_dfr(
      pairwise_comparisons,
      function(comp) {
        
        group1 <- comp[1]
        group2 <- comp[2]
        
        x <- test_dat %>%
          filter(TE_class == group1) %>%
          pull(value)
        
        y <- test_dat %>%
          filter(TE_class == group2) %>%
          pull(value)
        
        if (length(x) < opt$min_group_n || length(y) < opt$min_group_n) {
          return(
            tibble(
              feature = feat,
              group1 = group1,
              group2 = group2,
              test = "Wilcoxon rank-sum",
              n1 = length(x),
              n2 = length(y),
              pvalue = NA_real_,
              note = "Sample size lower than min_group_n"
            )
          )
        }
        
        wt <- wilcox.test(
          x,
          y,
          exact = FALSE
        )
        
        tibble(
          feature = feat,
          group1 = group1,
          group2 = group2,
          test = "Wilcoxon rank-sum",
          n1 = length(x),
          n2 = length(y),
          pvalue = wt$p.value,
          note = "OK"
        )
      }
    )
  }
)

write.table(
  pairwise_result,
  file = file.path(opt$outdir, "feature_pairwise_wilcox_test.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ============================================================
# 9. Boxplot with pairwise Wilcoxon p-values
# ============================================================

for (feat in feature_cols) {
  
  plot_dat <- merged_dat %>%
    select(GeneID, TE_class, value = all_of(feat)) %>%
    filter(!is.na(value), !is.na(TE_class)) %>%
    mutate(
      TE_class = factor(
        TE_class,
        levels = te_class_levels
      )
    )
  
  if (nrow(plot_dat) == 0) {
    next
  }
  
  class_n <- plot_dat %>%
    count(TE_class)
  
  valid_classes <- class_n %>%
    filter(n >= opt$min_group_n) %>%
    pull(TE_class) %>%
    as.character()
  
  valid_comparisons <- pairwise_comparisons[
    sapply(
      pairwise_comparisons,
      function(x) all(x %in% valid_classes)
    )
  ]
  
  p <- ggplot(
    plot_dat,
    aes(
      x = TE_class,
      y = value,
      fill = TE_class
    )
  ) +
    geom_boxplot(
      width = 0.5,
      outlier.size = 0.5,
      outlier.alpha = 0.5
    ) +
    #geom_jitter(
    #  width = 0.15,
    #  size = 0.4,
    #  alpha = 0.25
    #) +
    theme_bw(base_size = 13) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold"
      ),
      axis.text.x = element_text(
        angle = 30,
        hjust = 1
      ),
      legend.position = "none"
    ) +
    labs(
      #title = paste0(feat, " by log2TE class"),
      #subtitle = paste0("Condition: ", te_source, "; pairwise test: Wilcoxon rank-sum; "),
      x = "log2(TE) class",
      y = feat
    )
  
  if (length(valid_comparisons) > 0) {
    p <- p +
      stat_compare_means(
        comparisons = valid_comparisons,
        method = "wilcox.test",
        method.args = list(exact = FALSE),
        label = "p.format",
        hide.ns = FALSE,
        size = 4
      )
  }
  
  ggsave(
    filename = file.path(
      opt$outdir,
      "boxplot",
      paste0("boxplot_", safe_name(feat), ".pdf")
    ),
    plot = p,
    width = 7,
    height = 7
  )
}

message("Done. Results were saved to: ", opt$outdir)