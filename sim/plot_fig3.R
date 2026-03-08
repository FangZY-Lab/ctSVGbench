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
method_levels <- c("spVC_1","spVC_2","C-SIDE", "Celina", "STANCE","CTSV","ctSVG")

# Load custom theme and simulation functions
# source('F:/ctSVGbench/my_theme.R')
source('./my_theme.R')
source('./sim/utils/sim-bench-sc.R')

# Define datasets and patterns
dts_sc <- c(
  "StereoSeq_CBMSTA_Macaque1_T110",
  "StereoSeq_CBMSTA_Macaque1_T42",
  "StereoSeq_CBMSTA_Marmoset1_T478",
  "StereoSeq_CBMSTA_Marmoset1_T514",
  "StereoSeq_CBMSTA_Mouse1_T167",
  "StereoSeq_CBMSTA_Mouse1_T169",
  "StereoSeq_CBMSTA_Mouse1_T171",
  "StereoSeq_CBMSTA_Mouse1_T176",
  "StereoSeq_CBMSTA_Mouse1_T185",
  "StereoSeq_CBMSTA_Mouse1_T189",
  "StereoSeq_CBMSTA_Mouse2_T349",
  "VisiumHD_LUAD_2431", 
  "VisiumHD_LUAD_6123", 
  "VisiumHD_LUAD_6976", 
  "VisiumHD_LUSC_5488", 
  "VisiumHD_LUSC_7437", 
  "VisiumHD_LUSC_7941",
  "SeqFish+_cortex",
  "Stereoseq_mosta_E16.5_E1S3_whole_brain",
  "Stereoseq_mosta_Dorsal_midbrain"    
)#20

dts_sp <- c(
  "ST_PDAC",
  "Visium_liver",
  "Visium_mousebrain",
  "StereoSeq_MDESTA",
  "Visium_spleen",
  "SeqFish+_mouse_ob",
  "Slide-seq_tumor",
  "Slide-seqV2_hippocampus",
  "Slide-seqV2_mouseOB", 
  "Slide-seqV2_melanoma_GSM6025935_MBM05_rep1",
  "Slide-seqV2_melanoma_GSM6025936_MBM05_rep2",
  "Slide-seqV2_melanoma_GSM6025937_MBM05_rep3",
  "Slide-seqV2_melanoma_GSM6025938_MBM06",
  "Slide-seqV2_melanoma_GSM6025939_MBM07",
  "Slide-seqV2_melanoma_GSM6025940_MBM08",
  "Slide-seqV2_melanoma_GSM6025949_ECM08",
  "Slide-seqV2_melanoma_GSM6025950_ECM10"  
)#17

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

# ---- Calculate AUC for all datasets and patterns ----
paramset = 'P1'
auc_sc <- do.call(rbind, lapply(dts_sc, function(dt) {
  do.call(rbind, lapply(patterns, function(pt) {
    dataset <- sprintf("sim_%s-%s-%s-rep1", dt, pt, paramset)
    print(dataset)
    res <- get_pvalue_wide(dataset, svg_id,cell.level=T)
    dat.pval.wide <- res$dat.pval.wide
    label <- res$label[rownames(dat.pval.wide)]
    scores <- sapply(dat.pval.wide, function(pred) {
      roc(label, pred, levels = c(0, 1), direction = ">")$auc
    })
    
    auc_raw <- data.frame(
      score = scores,
      methods = names(scores),
      dataset = dataset
    )
    
    
    auc_raw <- auc_raw %>%
      mutate(
        dt = dt,
        pattern = pt
      )
    prepare_df(auc_raw, score_name = "auc", extract_info = FALSE)
  }))
}))

auc_sp <- do.call(rbind, lapply(dts_sp, function(dt) {
  do.call(rbind, lapply(patterns, function(pt) {
    dataset <- sprintf("sim_%s-%s-%s-rep1", dt, pt , paramset)
    print(dataset)
    res <- get_pvalue_wide(dataset, svg_id)
    dat.pval.wide <- res$dat.pval.wide
    label <- res$label[rownames(dat.pval.wide)]
    scores <- sapply(dat.pval.wide, function(pred) {
      roc(label, pred, levels = c(0, 1), direction = ">")$auc
    })
    
    auc_raw <- data.frame(
      score = scores,
      methods = names(scores),
      dataset = dataset
    )
    
    auc_raw <- auc_raw %>%
      mutate(
        dt = dt,
        pattern = pt
      )
    prepare_df(auc_raw, score_name = "auc", extract_info = FALSE)
  }))
}))
auc_sc <- as.data.frame(auc_sc)
auc_sp <- as.data.frame(auc_sp)
auc_sc$resolution  <- 'cell'
auc_sp$resolution  <- 'spot'
auc <- rbind(auc_sc,auc_sp)

# Order patterns by mean Celina AUC
celina_auc <- auc %>%
  filter(methods == "Celina") %>%
  group_by(pattern) %>%
  summarise(mean_auc = mean(auc, na.rm = TRUE))

auc <- auc %>%
  mutate(pattern = factor(pattern,levels=patterns))

auc$dt <- recode(auc$dt, "SeqFish+_mouse_ob" = "SeqFish+_mouse_OB")
table(auc$methods)

auc$methods <- factor(auc$methods,levels =method_levels)
auc$pattern <- stringr::str_to_sentence(auc$pattern)

# dataset_unique <- auc %>% 
#   distinct(dt, .keep_all = TRUE)

# dataset_unique <- dataset_unique %>%
#   mutate(
#     tech_platform = str_extract(dt, "^[^_]+"),
#      ) %>%
#   select(
#     dt,
#     tech_platform,  
#     resolution     
#   )

# write.csv(
#   dataset_unique,
#   file = "sim_data_info.csv",  
#   row.names = FALSE,            
#   na = ""                       
# )

source('./sim/pre_fig3B.R')

p1
# ---- Load required packages for density plots ----
library(reshape2)
library(ggplot2)
library(dplyr)
library(patchwork)
library(tidyr)


# ---- Load required packages for density plots ----
paramset='P1'
datasets <- expand.grid(dt = dts_sp, pt = patterns, stringsAsFactors = FALSE) %>%
  mutate(dataset = sprintf("sim_%s-%s-%s-rep1", dt, pt, paramset))

## 1) 
pval_long <- do.call(rbind, lapply(seq_len(nrow(datasets)), function(i) {
  ds <- datasets$dataset[i]
  
  res <- get_pvalue_wide(ds, svg_id,cell.level=F)   # wide: gene x method
  
  # wide -> long
  tmp <- reshape2::melt(as.matrix(res$dat.pval.wide))
  colnames(tmp) <- c("gene", "method", "pval")
  tmp$pval <- as.numeric(tmp$pval)
  
  tmp$dataset <- ds
  tmp
}))

## 2) 
pval_long <- pval_long %>%
  mutate(
    celltype = as.integer(str_extract(gene, "(?<=celltype)\\d+")),
    gene_id  = as.integer(str_extract(gene, "(?<=gene)\\d+"))
  ) %>% filter(celltype==4)

## pval_long 
pval_long <- pval_long %>%
  mutate(
    neg_class = case_when(
      celltype == 4 & gene_id >= 1   & gene_id <= 200 ~ "non-target ctSVG",
      # celltype == 5 & gene_id >= 76  & gene_id <= 150 ~ "affected",
      # celltype == 6 & gene_id >= 1   & gene_id <= 75  ~ "affected",
      celltype %in% c(4, 5, 6) & gene_id >= 200 & gene_id <= 1200 ~ "non-ctSVG",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(neg_class))

head(pval_long)
pval_long$method <- factor(pval_long$method,levels =method_levels)


plot_dendity <- pval_long %>%
  ggplot(aes(x = pval, color = neg_class, fill = neg_class)) +
  geom_density(alpha = 0.3) +
  labs(
    x = "p-value",
    y = "Density",
    color = "",
    fill = ""
  ) +
  theme_minimal() +
  facet_grid(~ method) +
  coord_cartesian(ylim = c(0,20)) +
  scale_x_continuous(
    labels = scales::label_number(accuracy = 0.1)
  )+
  scale_y_continuous(
    breaks = c(0,5,10,20),
    trans = pseudo_log_trans(base = 10)
  ) +
  my_theme +
  theme(legend.key.size = unit(0.05, "in"), 
        legend.key.width = unit(0.05, "in"),
        legend.key.height = unit(0.05, "in") )+
  scale_fill_manual(values=c("#008bd0", "#ffa61d")) +
  scale_color_manual(values=c("#008bd0", "#ffa61d"))

plot_dendity

# ---- Calculate TPR and FDR metrics for all datasets ----
alpha <- 0.05
datasets.all.sp <- expand.grid(dt = dts_sp, pattern = patterns, paramset = c("P1"), stringsAsFactors = FALSE) |>
  transform(dataset = sprintf("sim_%s-%s-%s-rep1", dt, pattern, paramset))

metric.df.sp <- do.call(rbind, lapply(seq_len(nrow(datasets.all.sp)), function(i) {
  dataset <- datasets.all.sp$dataset[i]
  pattern <- datasets.all.sp$pattern[i]
  dt <- datasets.all.sp$dt[i]
  
  res <- get_pvalue_wide(dataset, svg_id,cell.level=F)
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

datasets.all.sc <- expand.grid(dt = dts_sc, pattern = patterns, paramset = c("P1"), stringsAsFactors = FALSE) |>
  transform(dataset = sprintf("sim_%s-%s-%s-rep1", dt, pattern, paramset))

metric.df.sc <- do.call(rbind, lapply(seq_len(nrow(datasets.all.sc)), function(i) {
  dataset <- datasets.all.sc$dataset[i]
  pattern <- datasets.all.sc$pattern[i]
  dt <- datasets.all.sc$dt[i]
  
  res <- get_pvalue_wide(dataset, svg_id,cell.level=T)
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

metric.df <- rbind(metric.df.sc,metric.df.sp)
# ---- Plot sensitivity boxplots ----
metric.df$pattern <- stringr::str_to_sentence(metric.df$pattern)

metric.df$pattern <- factor(metric.df$pattern,levels = c("Pathology","Stripe","Hotspot","Gradient","Periodic","Neighbor" ))
metric.df$methods <- factor(metric.df$methods,levels =method_levels)

p4 <- metric.df %>%
  ggplot(aes(x = pattern, y = TPR, fill = methods)) +
  geom_boxplot(alpha = 1, width = 0.6, linewidth = 0.1, outlier.shape = NA, position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = method_colors) +
  my_theme +
  labs(fill = "", x = 'Pattern', y = "Sensitivity") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = "right",
        plot.margin = margin(0, 0.5, 0, 5),
        legend.key.size = unit(0.05, "in"),
        legend.key.width = unit(0.05, "in"),
        # legend.direction = "horizontal",
        # legend.box = "horizontal",     
        legend.margin = margin(1, 0, 0, 0))+
  guides(
    color = guide_legend(
      ncol = 1,                
      byrow = TRUE,            
      title.hjust = 0.5    ))
p4
# ---- Save workspace and combine plots ----

p5 <- p5 + theme(plot.margin = margin(0, 2, 1, 3),
                 legend.position = "none")

p45 <- plot_grid(p4, p5, nrow = 1, rel_widths = c(2, 1), labels = c("C", "D"),
                 label_x = -0.02,label_y = 1.06)

# Add left margin to heatmap plot
p1_with_margin <- ggdraw() +
  draw_plot(p1, x = 0.95/6, y = 0, width = 5.05/6, height = 1)

# Combine plots into final figure
final_plot <- plot_grid(p1_with_margin,  p45,plot_dendity,
                        nrow = 3,
                        rel_heights = c(1.7, 0.7,0.6),
                        labels = c("B", "","E"))
final_plot
ggsave('./Fig/Fig3-bcde.pdf', final_plot, width = 6.69, height = 6.95, units = "in")

# ---- Summarize metrics for exporting ----
summary_df_acc1 <- auc %>%
  # filter(resolution=="spot") %>% 
  group_by(methods, pattern) %>%
  summarise(auc = mean(score, na.rm = TRUE)) %>%
  pivot_wider(names_from = pattern, values_from = auc, names_glue = "{pattern}_auc") %>%
  column_to_rownames('methods')

summary_df_acc2 <- metric.df %>%
  group_by(methods,pattern) %>%
  summarise(sensitivity = median(TPR, na.rm = TRUE)) %>%
  pivot_wider(names_from = pattern, values_from = sensitivity, names_glue = "{pattern}_sensitivity") %>%
  column_to_rownames('methods')

summary_df <- read.csv("metrics_summary.csv", row.names = 1)

summary_df$gradient_sensitivity <- summary_df_acc2[rownames(summary_df),"Gradient_sensitivity"]
summary_df$hotspot_sensitivity <- summary_df_acc2[rownames(summary_df),"Hotspot_sensitivity"]
summary_df$neighbor_sensitivity <- summary_df_acc2[rownames(summary_df),"Neighbor_sensitivity"]
summary_df$pathology_sensitivity <- summary_df_acc2[rownames(summary_df),"Pathology_sensitivity"]
summary_df$periodic_sensitivity <- summary_df_acc2[rownames(summary_df),"Periodic_sensitivity"]
summary_df$stripe_sensitivity <- summary_df_acc2[rownames(summary_df),"Stripe_sensitivity"]

summary_df$pathology_auc <- summary_df_acc1[rownames(summary_df),"Pathology_auc"]
summary_df$stripe_auc <- summary_df_acc1[rownames(summary_df),"Stripe_auc"]
summary_df$hotspot_auc <- summary_df_acc1[rownames(summary_df),"Hotspot_auc"]
summary_df$gradient_auc <- summary_df_acc1[rownames(summary_df),"Gradient_auc"]
summary_df$periodic_auc <- summary_df_acc1[rownames(summary_df),"Periodic_auc"]
summary_df$neighbor_auc <- summary_df_acc1[rownames(summary_df),"Neighbor_auc"]

write.csv(summary_df, "metrics_summary.csv", row.names = TRUE)
