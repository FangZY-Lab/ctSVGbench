library(ggplot2)
library(dplyr)
library(cowplot)
library(tidyverse)
library(readxl)
library(hrbrthemes)
library(ggsci)
library(ggpubr)
library(ggtext)
library(ggprism)
library(ggpmisc)
library(pals)
library(viridis)
library(viridisLite)
library(scico)
library(reshape2)
library(ggrepel)
library(here)
library(tibble)
library(pROC)
library(aplot)
library(RColorBrewer)

# Load custom theme and simulation functions
source('F:/ctSVGbench/my_theme.R')
source('./sim/utils/sim-bench.R')

# Define datasets and patterns
dts <- c(
  "ST_PDAC",
  "Visium_liver",
  "Visium_mousebrain",
  "StereoSeq_MDESTA",
  "StereoSeq_CBMSTA_Macaque",
  "Visium_spleen",
  "Slide-seqV2_melanoma",
  "SeqFish+_mouse_ob",
  "StereoSeq_CBMSTA_Marmoset",
  "Slide-seq_tumor"
)
patterns <- c("pathology","hotspot", "stripe", "gradient",  "periodic", "neighbor")

# ---- Prepare data frame for metric calculation ----
prepare_df <- function(df, score_name, extract_info = TRUE) {
  if (extract_info) {
    df <- df %>%
      mutate(
        pattern = str_extract(dataset, "(?<=-)[^-]+(?=-)")
      )
  }
  
  df <- df %>%
    mutate(!!score_name := score) %>%
    group_by(dataset) %>%
    mutate(
      rank = rank(!!sym(score_name), ties.method = "average"),
      is_top = (rank == max(rank))
    ) %>%
    ungroup()
  
  return(df)
}

# ---- Calculate metrics (AUC or Jaccard) for all datasets ----
get_metric_df <- function(datasets, svg_id, metric = c("auc", "jaccard")) {
  metric <- match.arg(metric)
  
  calc_fn <- if (metric == "auc") {
    function(pred, label) roc(label, pred, levels = c(0, 1), direction = ">")$auc
  } else {
    function(pred, label) {
      pred <- as.logical(pred < 0.05)
      intersection <- sum(pred & label)
      union <- sum(pred | label)
      if (union == 0) return(NA)
      intersection / union
    }
  }
  
  do.call(rbind, lapply(datasets, function(dataset) {
    res <- get_pvalue_wide(dataset, svg_id)
    dat.pval.wide <- res$dat.pval.wide
    label <- res$label[rownames(dat.pval.wide)]
    
    scores <- sapply(dat.pval.wide, function(pred) calc_fn(pred, label))
    
    data.frame(
      score = scores,
      methods = names(scores),
      dataset = dataset
    )
  }))
}

# ---- Calculate AUC for all datasets and patterns ----
paramset = 'P1'
auc <- do.call(rbind, lapply(dts, function(dt) {
  do.call(rbind, lapply(patterns, function(pt) {
    dataset <- sprintf("sim_%s-%s-%s-rep1", dt, pt , paramset)
    auc_raw <- get_metric_df(dataset, svg_id, metric = "auc")
    
    auc_raw <- auc_raw %>%
      mutate(
        dt = dt,
        pattern = pt
      )
    
    prepare_df(auc_raw, score_name = "auc", extract_info = FALSE)
  }))
}))

auc <- as.data.frame(auc)

# Order patterns by mean Celina AUC
celina_auc <- auc %>%
  filter(methods == "Celina") %>%
  group_by(pattern) %>%
  summarise(mean_auc = mean(auc, na.rm = TRUE))

auc <- auc %>%
  mutate(pattern = factor(pattern, 
                          levels = celina_auc %>% arrange(desc(mean_auc)) %>% pull(pattern)))

auc$dt <- recode(auc$dt, "SeqFish+_mouse_ob" = "SeqFish+_mouse_OB")
auc$methods <- factor(auc$methods,levels = c("C-SIDE","spVC", "Celina", "STANCE"))
auc$pattern <- stringr::str_to_sentence(auc$pattern)
# ---- Plot AUC heatmap ----
p1 <- ggplot(auc, aes(x = methods, y = dt, fill = auc)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(
    colours = colorRampPalette(c("#008bd0", "#eeeeee", "#ffa61d"))(100),
    guide = guide_colorbar(
      title.position = "top",
      barwidth = unit(0.1, "in"),
      barheight = unit(1, "in"),
      label.hjust = 0.5
    )
  ) +
  facet_wrap(~pattern, ncol = 6) +
  theme_minimal() +
  my_theme +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.title.x = element_blank(),
    plot.margin = margin(0, 2, 0, 1),
    legend.key.size = unit(0.05, "in"), 
    legend.key.width = unit(0.05, "in"),
    legend.position = "right"
  ) +
  labs(x = "Method", y = "Simulated datasets", fill = "AUC")
p1

# ---- Load required packages for density plots ----
library(reshape2)
library(ggplot2)
library(dplyr)
library(patchwork)
library(tidyr)

# Define dataset, pattern, parameter set
dt <- "Slide-seqV2_melanoma"
pt <- "gradient"
paramset <- 'P1'

# Source simulation functions
source('./sim/utils/sim-bench.R')

dataset <- sprintf("sim_%s-%s-%s-rep1", dt, pt, paramset)
res <- get_pvalue_wide(dataset, svg_id)
dat.pval.wide <- res$dat.pval.wide

# ---- Prepare rank and p-value long tables ----
prepare_plot_density <- function(dat.pval.wide,
                                 celltype = 4,
                                 sel1.range = 1:200,
                                 sel2.range = 201:1200,
                                 rank_scope = c("global", "within")) {
  rank_scope <- match.arg(rank_scope)
  gene.prefix <- paste0("celltype", celltype, "gene")
  rows <- grep(paste0("^", gene.prefix), rownames(dat.pval.wide))
  if (length(rows) == 0) stop("No rows found for the requested celltype prefix.")
  dat_sub <- dat.pval.wide[rows, , drop = FALSE]
  
  # Define gene groups
  grp1_genes <- paste0(gene.prefix, sel1.range)
  grp2_genes <- paste0(gene.prefix, sel2.range)
  
  pval_long <- reshape2::melt(as.matrix(dat_sub))
  colnames(pval_long) <- c("gene", "method", "pval")
  pval_long$pval <- as.numeric(pval_long$pval)
  
  pval_long$group <- ifelse(pval_long$gene %in% grp1_genes,
                            paste0(min(sel1.range), "–", max(sel1.range)),
                            ifelse(pval_long$gene %in% grp2_genes,
                                   paste0(min(sel2.range), "–", max(sel2.range)),
                                   "other"))
  pval_long <- pval_long[pval_long$group != "other", , drop = FALSE]
  
  # Prepare rank long table
  if (rank_scope == "global") {
    all.rank <- apply(dat.pval.wide, 2, rank, ties.method = "average")
    rank_sub <- all.rank[rownames(dat_sub), , drop = FALSE]
  } else {
    rank_sub <- apply(dat_sub, 2, rank, ties.method = "average")
  }
  rank_long <- reshape2::melt(as.matrix(rank_sub))
  colnames(rank_long) <- c("gene", "method", "rank")
  
  rank_long$group <- ifelse(rank_long$gene %in% grp1_genes,
                            paste0(min(sel1.range), "–", max(sel1.range)),
                            ifelse(rank_long$gene %in% grp2_genes,
                                   paste0(min(sel2.range), "–", max(sel2.range)),
                                   "other"))
  rank_long <- rank_long[rank_long$group != "other", , drop = FALSE]
  
  list(rank_df = as.data.frame(rank_long),
       pval_df = as.data.frame(pval_long))
}

# ---- Set plotting parameters ----
my_colors <- c("#008bd0", "#ffa61d")
base_font_size <- 7
title_font_size <- 7

# ---- Prepare data for celltype 4 ----
prep <- prepare_plot_density(dat.pval.wide,
                             celltype = 4,
                             sel1.range = 1:200,
                             sel2.range = 201:1200,
                             rank_scope = "global")
pval_df <- prep$pval_df

# Rank density plot for celltype 4
p_rank <- ggplot(pval_df, aes(x = pval, fill = group, color = group)) +
  geom_density(alpha = 0.35, size = 0.6) +
  facet_wrap(~method, scales = "free", ncol = 4) +
  theme_minimal(base_size = base_font_size) +
  scale_fill_manual(values = my_colors) +
  scale_color_manual(values = my_colors) +
  labs(x = "P value", y = "Density", fill = "gene id", color = "gene id") +
  my_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right",
        legend.key.size = unit(0.05, "in"),
        legend.key.width = unit(0.05, "in"),
        plot.margin = margin(4, 2, 4, 0))
p_rank

# ---- Prepare data for celltype 6 ----
prep_ct6 <- prepare_plot_density(dat.pval.wide,
                                 celltype = 6,
                                 sel1.range = 1:75,
                                 sel2.range = 201:1200,
                                 rank_scope = "global")
pval_df_ct6 <- prep_ct6$pval_df

# Rank density plot for celltype 6
p_rank_ct6 <- ggplot(pval_df_ct6, aes(x = pval, fill = group, color = group)) +
  geom_density(alpha = 0.35, size = 0.6) +
  facet_wrap(~method, scales = "free", ncol = 4) +
  theme_minimal(base_size = base_font_size) +
  scale_fill_manual(values = my_colors) +
  scale_color_manual(values = my_colors) +
  labs(x = "P value", y = "Density", fill = "gene id", color = "gene id") +
  my_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right",
        legend.key.size = unit(0.05, "in"),
        legend.key.width = unit(0.05, "in"),
        plot.margin = margin(4, 2, 4, 0))

# ---- Calculate TPR and FDR metrics for all datasets ----
alpha <- 0.05
datasets.all <- expand.grid(dt = dts, pattern = patterns, paramset = c("P1"), stringsAsFactors = FALSE) |>
  transform(dataset = sprintf("sim_%s-%s-%s-rep1", dt, pattern, paramset))

metric.df <- do.call(rbind, lapply(seq_len(nrow(datasets.all)), function(i) {
  dataset <- datasets.all$dataset[i]
  pattern <- datasets.all$pattern[i]
  dt <- datasets.all$dt[i]
  
  res <- get_pvalue_wide(dataset, svg_id)
  dat.pval.wide <- res$dat.pval.wide
  label <- res$label[rownames(dat.pval.wide)]
  
  do.call(rbind, lapply(names(dat.pval.wide), function(name) {
    pval <- dat.pval.wide[[name]]
    sig <- (pval < alpha)
    TP <- sum(sig & label == 1)
    FP <- sum(sig & label == 0)
    FN <- sum(!sig & label == 1)
    TPR <- ifelse((TP + FN) > 0, TP / (TP + FN), NA)
    FDR <- ifelse((TP + FP) > 0, FP / (TP + FP), NA)
    NUM <- sum(sig)
    data.frame(TPR = TPR, FDR = FDR, NUM = NUM, methods = name, dataset = dataset, pattern = pattern, dt = dt, stringsAsFactors = FALSE)
  }))
}))

# ---- Plot sensitivity boxplots ----
metric.df$methods <- factor(metric.df$methods,levels = c("C-SIDE","spVC", "Celina", "STANCE"))
metric.df$pattern <- stringr::str_to_sentence(metric.df$pattern)

metric.df$pattern <- factor(metric.df$pattern,levels = c("Pathology","Stripe","Hotspot","Gradient","Periodic","Neighbor" ))

p4 <- metric.df %>%
  ggplot(aes(x = pattern, y = TPR, fill = methods)) +
  geom_boxplot(alpha = 1, width = 0.6, linewidth = 0.1, outlier.shape = NA, position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = method_colors) +
  my_theme +
  labs(fill = "", x = 'Pattern', y = "Sensitivity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top",
        plot.margin = margin(0, 1.5, 0, 1),
        legend.key.size = unit(0.05, "in"),
        legend.key.width = unit(0.05, "in"))
p4
# ---- Save workspace and combine plots ----
save.image('./plot/sim.rda')
p5 <- readRDS('./plot/sim_null.rds')
p5 <- p5 + theme(plot.margin = margin(16, 1.5, 0, 5))

p45 <- plot_grid(p4, p5, nrow = 1, rel_widths = c(1.75, 1), labels = c("C", "D"))

# Add left margin to heatmap plot
p1_with_margin <- ggdraw() +
  draw_plot(p1, x = 1/6, y = 0, width = 5/6, height = 1)

# Combine plots into final figure
final_plot <- plot_grid(p1_with_margin, p45, p_rank, p_rank_ct6,
                        nrow = 4,
                        rel_heights = c(0.9, 1,0.7, 0.7),
                        labels = c("B", "","E", "F"))
final_plot
ggsave('./Fig/Fig3-1.pdf', final_plot, width = 6.69, height = 7, units = "in")

# ---- Summarize metrics for exporting ----
summary_df_acc1 <- auc %>%
  group_by(methods, pattern) %>%
  summarise(auc = mean(score, na.rm = TRUE)) %>%
  pivot_wider(names_from = pattern, values_from = auc, names_glue = "{pattern}_auc") %>%
  column_to_rownames('methods')

summary_df_acc2 <- metric.df %>%
  group_by(methods,pattern) %>%
  summarise(sensitivity = mean(TPR, na.rm = TRUE)) %>%
  pivot_wider(names_from = pattern, values_from = sensitivity, names_glue = "{pattern}_sensitivity") %>%
  column_to_rownames('methods')

summary_df <- read.csv("metrics_summary.csv", row.names = 1)
summary_df$gradient_sensitivity <- summary_df_acc2[rownames(summary_df),"gradient_sensitivity"]
summary_df$hotspot_sensitivity <- summary_df_acc2[rownames(summary_df),"hotspot_sensitivity"]
summary_df$neighbor_sensitivity <- summary_df_acc2[rownames(summary_df),"neighbor_sensitivity"]
summary_df$pathology_sensitivity <- summary_df_acc2[rownames(summary_df),"pathology_sensitivity"]
summary_df$periodic_sensitivity <- summary_df_acc2[rownames(summary_df),"periodic_sensitivity"]
summary_df$stripe_sensitivity <- summary_df_acc2[rownames(summary_df),"stripe_sensitivity"]
  
summary_df$pathology_auc <- summary_df_acc1[rownames(summary_df),"pathology_auc"]
summary_df$stripe_auc <- summary_df_acc1[rownames(summary_df),"stripe_auc"]
summary_df$hotspot_auc <- summary_df_acc1[rownames(summary_df),"hotspot_auc"]
summary_df$gradient_auc <- summary_df_acc1[rownames(summary_df),"gradient_auc"]
summary_df$periodic_auc <- summary_df_acc1[rownames(summary_df),"periodic_auc"]
summary_df$neighbor_auc <- summary_df_acc1[rownames(summary_df),"neighbor_auc"]

write.csv(summary_df, "metrics_summary.csv", row.names = TRUE)
