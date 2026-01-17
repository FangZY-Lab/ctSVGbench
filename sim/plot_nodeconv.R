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
source('./my_theme.R')
source('./sim/utils/sim-bench.R')
source('./sim/utils/sim-bench-nodeconv.R')


# Define datasets and patterns
dts <- c(
  "ST_PDAC",
  "Visium_mousebrain",
  "StereoSeq_MDESTA",
  "Slide-seq_tumor",
  "Slide-seqV2_hippocampus",
  "StereoSeq_CBMSTA_Macaque"
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

# ---- Calculate AUC for all datasets and patterns ----
paramset = 'P1'
auc_nodcv <- do.call(rbind, lapply(dts, function(dt) {
  do.call(rbind, lapply(patterns, function(pt) {
    dataset <- sprintf("sim_%s-%s-%s-rep1-noRCTD", dt, pt, paramset)
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

auc_dcv <- do.call(rbind, lapply(dts, function(dt) {
  do.call(rbind, lapply(patterns, function(pt) {
    dataset <- sprintf("sim_%s-%s-%s-rep1", dt, pt, paramset)
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

auc_nodcv <- as.data.frame(auc_nodcv)
auc_dcv <- as.data.frame(auc_dcv)
auc_nodcv$dcv  <- 'no'
auc_dcv$dcv  <- 'yes'
auc <- rbind(auc_nodcv,auc_dcv)

# Order patterns by mean Celina AUC
celina_auc <- auc %>%
  filter(methods == "Celina") %>%
  group_by(pattern) %>%
  summarise(mean_auc = mean(auc, na.rm = TRUE))

auc <- auc %>%
  mutate(pattern = factor(pattern, 
                          levels = celina_auc %>% arrange(desc(mean_auc)) %>% pull(pattern)))

auc$dt <- recode(auc$dt, "StereoSeq_CBMSTA_Macaque" = "StereoSeq_CBMSTA_Maca")
auc$methods <- factor(auc$methods,levels = c("C-SIDE","spVC", "Celina", "STANCE","CTSV"))
auc$pattern <- stringr::str_to_sentence(auc$pattern)
# ---- Plot AUC heatmap ----

ggplot(auc, aes(x = dcv, y = rank, color = methods, group = methods)) +
  geom_line(alpha = 0.5) +
  geom_point() +
  labs(x="Using RCTD deconvolution",y="Rank(higher = better)",fill="Methods")+
  scale_color_manual(values = method_colors)+
  facet_wrap(dt~pattern, ncol = 6) +
  theme_minimal() +
  my_theme +
  theme(strip.text = element_text(size = 5))

ggsave('Fig/s/noRCTD.pdf',width = 6.69,height = 9)


source('./sim/utils/calc_false_positive_rate.R')

dts_sim <- c(
  "ST_PDAC",
  "Visium_mousebrain",
  "StereoSeq_MDESTA",
  "Slide-seq_tumor",
  "Slide-seqV2_hippocampus",
  "StereoSeq_CBMSTA_Macaque"
)

fpr_df_norctd <- do.call(rbind, lapply(dts_sim, function(dt) {
  do.call(rbind, lapply(patterns, function(pt) {
    dataset <- sprintf("sim_%s-%s-%s-rep1-noRCTD", dt, pt, paramset)
    res <- get_pvalue_wide(dataset, svg_id)
    fpr <- calc_false_positive_rate(res$dat.pval.wide, alpha = 0.05)
    fpr$dataset <- dataset
    fpr$pattern <-pt 
    fpr
  }))
}))

fpr_df_rctd <- do.call(rbind, lapply(dts_sim, function(dt) {
  do.call(rbind, lapply(patterns, function(pt) {
    dataset <- sprintf("sim_%s-%s-%s-rep1", dt, pt, paramset)
    res <- get_pvalue_wide(dataset, svg_id)
    fpr <- calc_false_positive_rate(res$dat.pval.wide, alpha = 0.05)
    fpr$dataset <- dataset
    fpr$pattern <-pt 
    fpr
  }))
}))
fpr_df_rctd$dcv  <- 'yes'
fpr_df_norctd$dcv  <- 'no'
fpr_df_dcv <- rbind(fpr_df_rctd,fpr_df_norctd)
library(ggplot2)
library(dplyr)
library(stringr)

fpr_df_dcv <- fpr_df_dcv %>%
  mutate(pattern = factor(pattern,levels=patterns))

head(fpr_df_dcv)
library(dplyr)
library(ggplot2)
library(tidyr)
fpr_plot_dcv <- fpr_df_dcv

ggplot(fpr_plot_dcv, aes(x = method, y = fpr, fill = dcv, group = interaction(method, dcv))) +
  geom_boxplot() +
  facet_wrap(~pattern, ncol = 3) +  
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +  
  labs(y = "FPR", x = "Methods", fill = "Using RCTD deconvolution") +
  scale_fill_manual(values = c("#FDB462", "#B3DE69")) +
  theme_minimal() + 
  my_theme+
  theme(
  axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 5)
  ) +
  coord_cartesian(ylim = c(0, 0.2)) +  
  theme(legend.position = "right")

ggsave('Fig/s/affect_dcv_boxplot.pdf', width = 6.69, height = 3, dpi = 300)
