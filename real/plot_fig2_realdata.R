# Load required packages
library(ggplot2)
library(scales)
library(tidyr)
library(dplyr)
library(tibble)
library(tidyverse)
library(data.table)
library(cowplot)
library(extrafont)
library(here)
library(ggrepel)
library(grid)
library(ggrepel)
library(RColorBrewer)
source(here('real', 'utils', "real-bench.R"))
source('F:/ctSVGbench/my_theme.R')

# List of datasets used in the analysis
datasets <- c(
  "MERFISH_hypothalamus",
  "SeqFish+_mouse_ob",
  "Slide-seq_tumor",
  "Slide-seqV2_hippocampus",
  "Slide-seqV2_melanoma",
  "Slide-seqV2_mouseOB",
  "ST_developmental_heart",
  "ST_PDAC",
  "StereoSeq_CBMSTA_Macaque",
  "StereoSeq_CBMSTA_Marmoset",
  "StereoSeq_MDESTA",
  "StereoSeq_mouseOB",
  "Visium_bladder",
  "Visium_intestine",
  "Visium_liver",
  "Visium_lymph_node",
  "Visium_melanoma",
  "Visium_mousebrain",
  "Visium_pancreas",
  "Visium_skin",
  "Visium_spleen",
  "Visium_tail",
  "VisiumHD_LUSC_AXB-5488-D1"
  )

##------------------
## Number of detected genes per method
##------------------
svnum.df <- do.call(rbind, lapply(datasets, function(dataset) {
  dat.pval.wide <- get_wide_pval(dataset)
  dat.pval.binary <- as.data.frame((dat.pval.wide < 0.05) * 1)
  method_counts <- colSums(dat.pval.binary)
  df <- data.frame(
    method = names(method_counts),
    count = as.numeric(method_counts),
    dataset = dataset
  )
  return(df)
}))
svnum.df$dataset <- gsub("VisiumHD_LUSC_AXB-5488-D1","VisiumHD_lung_cancer",svnum.df$dataset)

# Sort datasets based on STANCE counts
dataset_order <- svnum.df %>%
  filter(method == "STANCE") %>%
  arrange(count) %>%
  pull(dataset)

svnum.df$dataset <- factor(svnum.df$dataset, levels = dataset_order)

n <- length(unique(svnum.df$dataset))  # 23 datasets
my_colors <- colorRampPalette(brewer.pal(12, "Set3"))(n)
ds_levels <- levels(svnum.df$dataset)

# Background rectangles for dataset grouping
rect_df <- data.frame(
  x = seq_along(ds_levels),
  dataset = ds_levels,
  xmin = seq_along(ds_levels) - 0.4,
  xmax = seq_along(ds_levels) + 0.4,
  ymin = 0,
  ymax = 5,
  fill = my_colors[seq_along(ds_levels)]
)

merged_df <- merge(svnum.df, rect_df, by = "dataset")
merged_df$dataset <- factor(merged_df$dataset, levels = dataset_order)
fill_colors <- setNames(merged_df$fill[match(dataset_order, merged_df$dataset)], dataset_order)
# Plot: detected gene count by dataset and method
p1 <- ggplot() +
  geom_rect(data = merged_df,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = dataset),
            color = NA, inherit.aes = FALSE) +
  scale_fill_manual(name = "Dataset", values = fill_colors) +
  geom_segment(data = merged_df, aes(x = x, y = 0, yend = count, color = method), size = 0.3) +
  geom_point(data = merged_df, aes(x = x, y = count, color = method), size = 1) +
  facet_wrap(~method) +
  scale_color_manual(values = method_colors, guide = "none") +
  labs(x = "Datasets", y = "Detected gene count") +
  theme_minimal() +
  my_theme +
  scale_y_sqrt() +
  theme(
    legend.position = 'right',
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(2, 0, 0, 0),
    legend.key.size = unit(0.05, "in"),
    legend.key.width = unit(0.05, "in")
  ) +
  guides(fill = guide_legend(ncol = 1))

p1

##------------------
## Overlap proportion between methods
##------------------

alpha <- 0.05
calc_overlap_by_method_partition <- function(pmat, method, alpha = 0.05) {
  sig <- pmat < alpha
  sig[is.na(sig)] <- FALSE
  
  # Count how many methods identify each gene as significant
  k <- rowSums(sig, na.rm = TRUE)
  
  # Get genes detected by the given method (for denominator)
  mask <- sig[, method]
  denom <- sum(mask, na.rm = TRUE)
  
  # If no significant genes, return zeros
  if (denom == 0) {
    return(tibble(
      overlap = factor(c("1 method", "2 methods", "3 methods", "4 methods"),
                       levels = c("1 method", "2 methods", "3 methods", "4 methods")),
      prop = 0,
      method = method
    ))
  }
  
  # Count overlap levels (only within selected method)
  tab <- table(factor(k[mask], levels = 1:4))
  prop <- as.numeric(tab) / denom
  
  tibble(
    overlap = factor(c("1 method", "2 methods", "3 methods", "4 methods"),
                     levels = c("1 method", "2 methods", "3 methods", "4 methods")),
    prop = prop,
    method = method
  )
}

# Combine overlap data for all datasets
df_plot <- bind_rows(lapply(datasets, function(ds) {
  pmat <- get_wide_pval(ds)
  bind_rows(lapply(colnames(pmat), function(m) {
    calc_overlap_by_method_partition(pmat, m, alpha) %>% mutate(dataset = ds)
  }))
}))
df_plot$dataset <- gsub("VisiumHD_LUSC_AXB-5488-D1","VisiumHD_lung_cancer",df_plot$dataset)

df_plot$dataset <- factor(df_plot$dataset, levels = dataset_order)
df_plot$method <- factor(df_plot$method, levels = c("C-SIDE","spVC", "Celina", "STANCE"))
# Plot: overlap proportions
p2 <- ggplot(df_plot, aes(x = dataset, y = prop, fill = overlap)) +
  geom_col(position = "stack") +
  facet_wrap(~ method, nrow = 1) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.02))) +
  scale_fill_manual(values = c("1 method" = "#c6dbef",
                               "2 methods" = "#9ecae1",
                               "3 methods" = "#6baed6",
                               "4 methods" = "#2171b5")) +
  labs(x = "Datasets", y = "Proportion (within method)", fill = "Overlap size") +
  theme_minimal() +
  my_theme +
  theme(
    legend.position = 'bottom',
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(2, 0, 0, 0),
    legend.key.size = unit(0.05, "in"),
    legend.key.width = unit(0.05, "in")
  ) +
  guides(fill = guide_legend(nrow = 1))
p2

##------------------
## Concordance between methods
##------------------

dts <- c(
  'MERFISH_hypothalamus',
  "ST_PDAC",
  "Visium_spinal",
  "Slide-seqV2_melanoma",
  "Visium_skin"
)
conc.res <- do.call(rbind, lapply(datasets, function(dataset) {
  dat.pval.wide <- get_wide_pval(dataset)
  data <- get_conc(dat.pval.wide, dataset = dataset)
  return(data)
}))
conc.res$dataset <- gsub("VisiumHD_LUSC_AXB-5488-D1","VisiumHD_lung_cancer",conc.res$dataset)

conc.res$method1 <- factor(conc.res$method1,levels = c("C-SIDE","spVC", "Celina", "STANCE"))
conc.res1 <- subset(conc.res, dataset %in% dts)
conc.res2 <- subset(conc.res, !(dataset %in% dts))

# Plot: top concordance curves
p3 <- ggplot(conc.res1, aes(x = rank, y = conc, color = method2)) +
  geom_line(linewidth = 0.3) +
  facet_grid(method1 ~ dataset) +
  theme_minimal() +
  my_theme +
  scale_color_manual(values = method_colors) +
  labs(title = "", y = "Concordance at the top", x = 'Rank', color = "") +
  theme(
    strip.text = element_text(color = "black", size = 7),
    axis.title.y = element_text(size = 7),
    axis.title.x = element_text(size = 7),
    legend.position = "none",
    plot.margin = margin(2, 0, 0, 0)
  )
p3


##------------------
## Rank correlation between methods (Top 100 genes)
##------------------

library(dplyr)

get_top100_rank_corr <- function(dataset) {
  # Load the p-value matrix (genes × methods)
  dat.pval.wide <- get_wide_pval(dataset)
  
  # Compute rank for each method
  rank_mat <- apply(dat.pval.wide, 2, rank, ties.method = "average")
  
  # Extract the top 100 genes for each method
  top100_genes <- apply(rank_mat, 2, function(r) {
    names(sort(r))[1:100]
  })
  
  # Prepare a correlation matrix (method × method)
  methods <- colnames(rank_mat)
  cor_mat <- matrix(NA, nrow = length(methods), ncol = length(methods),
                    dimnames = list(methods, methods))
  
  # Calculate pairwise Spearman correlation between methods
  for (i in seq_along(methods)) {
    for (j in seq_along(methods)) {
      g1 <- top100_genes[, i]
      g2 <- top100_genes[, j]
      common <- intersect(g1, g2)
      cor_mat[i, j] <- cor(rank_mat[common, i], rank_mat[common, j], method = "spearman")
    }
  }
  
  # Return long-format data frame
  return(as.data.frame(as.table(cor_mat)) %>%
           rename(method1 = Var1, method2 = Var2, correlation = Freq) %>%
           mutate(dataset = dataset))
}

# Run correlation calculation for all datasets
res.list <- lapply(datasets, get_top100_rank_corr)
long_cor <- do.call(rbind, res.list)
long_cor$dataset <- gsub("VisiumHD_LUSC_AXB-5488-D1","VisiumHD_lung_cancer",long_cor$dataset)
long_cor$method1 <- factor(long_cor$method1,levels = c("C-SIDE","spVC", "Celina", "STANCE"))
long_cor$method2 <- factor(long_cor$method2,levels = c("C-SIDE","spVC", "Celina", "STANCE"))

# Boxplot: pairwise correlation between methods
p4 <- as.data.frame(long_cor) %>% 
  filter(method1 != method2) %>%
  ggplot(aes(x = method1, y = correlation, fill = method2)) +
  geom_boxplot(alpha = 1, width = 0.6, color = "grey30") +
  scale_fill_manual(values = method_colors) +
  theme_test() +
  my_theme +
  theme(
    plot.margin = margin(2, 2, 0, 0),
    axis.title.y = element_text(size = 7),
    legend.position = 'bottom'
  ) +
  labs(x = "", y = "Rank correlation", fill = "")
p4
ggsave(sprintf("./Fig/s/S_real_cor1.pdf"), width = 6.69, height = 3)

# Heatmap: correlation matrix for each dataset
p <- ggplot(long_cor, aes(x = method1, y = method2, fill = correlation)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                       name = "Spearman\nCorrelation") +
  theme_minimal() +
  my_theme +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    panel.grid = element_blank()
  ) +
  coord_fixed() +
  facet_wrap(~ dataset, ncol = 4)
p
ggsave(sprintf("./Fig/s/S_real_cor2.pdf"), width = 6.69, height = 9.7)


##------------------
## Jaccard similarity between methods
##------------------

compute_jaccard <- function(dat.pval, mode = c("pval", "top"), alpha = 0.05, topn = 100) {
  mode <- match.arg(mode)
  methods <- colnames(dat.pval)
  
  # Identify significant gene sets for each method
  sig.list <- lapply(methods, function(m) {
    pvals <- dat.pval[, m]
    names(pvals) <- rownames(dat.pval)
    
    if (mode == "pval") {
      # Genes with p-value < alpha
      genes <- names(pvals)[which(pvals < alpha)]
    } else if (mode == "top") {
      # Top-N ranked genes
      genes <- names(sort(pvals))[1:min(topn, length(pvals))]
    }
    return(genes)
  })
  names(sig.list) <- methods
  
  # Compute Jaccard index between all method pairs
  n <- length(methods)
  jaccard.mat <- matrix(NA, n, n, dimnames = list(methods, methods))
  
  for (i in 1:n) {
    for (j in 1:n) {
      A <- sig.list[[i]]
      B <- sig.list[[j]]
      jaccard.mat[i, j] <- length(intersect(A, B)) / length(union(A, B))
    }
  }
  return(jaccard.mat)
}

# Run Jaccard analysis for all datasets and significance criteria
jaccard.df <- do.call(rbind, lapply(datasets, function(dataset) {
  dat.pval.wide <- get_wide_pval(dataset)
  
  # Define parameter sets for computation
  param.list <- list(
    list(mode = "pval", alpha = 0.05, topn = NA, label = "adjusted p < 0.05"),
    list(mode = "pval", alpha = 0.01, topn = NA, label = "adjusted p < 0.01"),
    list(mode = "top", alpha = NA, topn = 20, label = "top 20"),
    list(mode = "top", alpha = NA, topn = 50, label = "top 50")
  )
  
  do.call(rbind, lapply(param.list, function(par) {
    jmat <- compute_jaccard(
      dat.pval.wide,
      mode = par$mode,
      alpha = ifelse(is.na(par$alpha), 0.05, par$alpha),
      topn = ifelse(is.na(par$topn), 100, par$topn)
    )
    jdf <- reshape2::melt(jmat, varnames = c("method1", "method2"), value.name = "Jaccard")
    jdf$dataset <- dataset
    jdf$mode <- par$label
    return(jdf)
  }))
}))

# Remove self-comparisons
jaccard.df <- subset(jaccard.df, method1 != method2)
jaccard.df <- rename(jaccard.df, method = method1)
jaccard.df$dataset <- gsub("VisiumHD_LUSC_AXB-5488-D1","VisiumHD_lung_cancer",jaccard.df$dataset)

# Boxplot: Jaccard index distribution across datasets
p5 <- ggplot(jaccard.df, aes(x = method2, y = Jaccard, fill = method)) +
  geom_boxplot(alpha = 1, width = 0.6, linewidth = 0.1, outlier.shape = NA) +
  scale_fill_manual(values = method_colors) +
  facet_wrap(~ mode, ncol = 2) +
  theme_minimal() +
  my_theme +
  ylim(0, 0.4) +
  labs(x = "", y = "Jaccard index", fill = "") +
  theme(
    axis.title.y = element_text(size = 7),
    axis.text.x = element_text(size = 7),
    legend.position = "bottom",
    plot.margin = margin(2, 1.5, 0, 2.5),
    legend.key.size = unit(0.05, "in"),
    legend.key.width = unit(0.05, "in")
  ) +
  guides(guide_legend(nrow = 1))
p5

# Combine multiple panels into a composite figure
p45 <- plot_grid(
  plotlist = list(p2, p5),
  rel_widths = c(1, 1),
  ncol = 2,
  labels = c("B", "C"),
  label_x = -0.01,
  label_y = 1.01
)

plot_grid(
  plotlist = list(p1, p45, p3),
  rel_heights = c(1.1, 1, 1.4),
  nrow = 3,
  labels = c("A", "", "D"),
  label_x = -0.005,
  label_y = 1.00
)
ggsave("./Fig/Fig2.pdf", width = 6.69, height = 9)

# Save workspace
save.image("./plot/real_plot.rda")

##------------------
## Supplemental concordance plots (paginated)
##------------------

i <- 1
start <- 1

while (i <= 4) {
  pages <- ifelse(i == 4, 4, 5)
  
  subset_data <- conc.res2 %>% 
    subset(dataset %in% unique(conc.res2$dataset)[start:(start + pages - 1)])
  
    pS3 <- ggplot(subset_data, aes(x = rank, y = conc, color = method2)) +
    geom_line() +
    facet_wrap(dataset ~ method1, nrow = pages, ncol = 4) +
    theme_minimal() +
    my_theme +
    scale_color_manual(values = method_colors) +
    labs(title = "", y = "Concordance at the top", color = "") +
    guides(color = guide_legend(override.aes = list(size = 4))) +
    theme(
      strip.text = element_text(color = "black", size = 7),
      legend.position = "none"
    )
  
  ggsave(sprintf("./Fig/s/conc_methods%s.pdf", i), width = 6.9, height = 8)
  
    start <- start + pages
  i <- i + 1
}


