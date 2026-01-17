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
# source('F:/ctSVGbench/my_theme.R')
source('./my_theme.R')
source('./sim/utils/sim-bench-sc.R')


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
  "stereoseq_mosta_E16.5_E1S3_whole_brain",
  "stereoseq_mosta_Dorsal_midbrain"    
) # 20

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
) # 17

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

auc$methods <- factor(auc$methods,levels = c("C-SIDE", "spVC", "Celina", "STANCE", "CTSV", "ctSVG"))
auc$pattern <- stringr::str_to_sentence(auc$pattern)

data_info=read.csv('./data_info.csv')
data_info <- data_info %>% 
  select(
    !resolution     
  )
merged_data <- merge(auc, data_info, by = "dt", all.x = TRUE)

library(ggplot2)
library(ggnewscale)  
library(RColorBrewer) 
library(dplyr)        

# ========== 1. get NA and arrange  ==========
merged_data <- merged_data %>%
  mutate(method_pattern = paste(pattern, methods, sep = ": ")) %>% 
  complete(dt, method_pattern) %>%
  group_by(dt) %>%
  fill(resolution, tech_platform, species, tissue, .direction = "downup") %>%
  ungroup() %>%
  mutate(
    auc = as.numeric(auc),
    pattern = ifelse(
      is.na(pattern),
      sub(":.*", "", method_pattern),
      pattern
    ),
    methods = ifelse(
      is.na(methods),
      sub(".*:\\s*", "", method_pattern),
      methods
    )
  )
merged_data <- merged_data %>%
  arrange(pattern, methods) %>% 
  mutate(
    method_pattern = factor(method_pattern, levels = unique(method_pattern)),
    pattern_x = method_pattern
  ) %>% 
  arrange(resolution, tech_platform, species, tissue,dt) %>% 
  mutate(dt = factor(dt, levels = unique(dt)))

merged_data <- merged_data %>%
  mutate(pattern = factor(pattern,levels=str_to_sentence(patterns))) %>% 
  arrange(pattern,methods,method_pattern) %>% 
  mutate(method_pattern = factor(method_pattern, levels = unique(method_pattern)))


merged_data <- merged_data %>%
  arrange(resolution, tech_platform, species, tissue, dt) %>% 
  mutate(dt = factor(dt, levels = unique(dt)))


# 1.2 legend prepare
dt_annotations <- merged_data %>%
  distinct(dt, resolution, tech_platform, species, tissue)

pattern_annotations <- merged_data %>%
  distinct(method_pattern, pattern, pattern_x)

# ========== 2. paltette prepare ==========
palette_mapping <- list(
  resolution = "Set3",       
  tech_platform = "Set2",    
  species = "Paired",        
  tissue = "Dark2",          
  pattern = "Set1"           
)

get_dim_color_palette <- function(categories, palette) {
  n <- length(categories)
  max_col <- brewer.pal.info[palette, "maxcolors"]
  
  if (n > max_col) {
    colors <- colorRampPalette(brewer.pal(max_col, palette))(n)
  } else {
    colors <- brewer.pal(n, palette)
  }
  setNames(colors, categories)
}


res_colors <- get_dim_color_palette(unique(dt_annotations$resolution), palette_mapping$resolution)
tech_colors <- get_dim_color_palette(unique(dt_annotations$tech_platform), palette_mapping$tech_platform)
species_colors <- get_dim_color_palette(unique(dt_annotations$species), palette_mapping$species)
tissue_colors <- get_dim_color_palette(unique(dt_annotations$tissue), palette_mapping$tissue)
pattern_colors <- get_dim_color_palette(unique(pattern_annotations$pattern), palette_mapping$pattern)

# ========== 3. plot heatmap ==========
annotation_width <- 1
offset_base <- -annotation_width
pattern_anno_height <- 0.6
pattern_offset_base <- -pattern_anno_height
anno_width <- 1
anno_height <- 1
na_df <- merged_data %>% filter(is.na(auc)) #crosses and grey for NA
na_df <- as.data.frame(na_df)

p1 <- ggplot() +
  geom_tile(
    data = merged_data,
    aes(
      x = method_pattern,
      y = dt,
      fill = auc
    ),
    color = "white",
    width = 1,
    height = 1
  ) +
  scale_fill_gradientn(
    colours = colorRampPalette(c("#008bd0", "#eeeeee", "#ffa61d"))(100),
    na.value = "#e0e0e0",
    limits = c(0.3, max(merged_data$auc, na.rm = TRUE)),
    breaks = c(0.3, 1),
    guide = guide_colorbar(
      title.position = "top",
      barwidth = unit(0.05, "in"),
      barheight = unit(0.5, "in"),
      label.hjust = 0.5
    ),
    name = "AUC"
  ) +
  
  new_scale_fill() +
  geom_tile(
    data = dt_annotations,
    aes(x = offset_base, y = dt, fill = resolution),
    width = annotation_width,
    color = NA,
    linewidth = 0,
    height = 1,
    position = "identity"
  ) +
  scale_fill_manual(values = res_colors, name = "Resolution") +
  guides(
    fill = guide_legend(
      nrow = 3,
      byrow = TRUE,
      order = 3  #legend order
    )
  ) +
  
  new_scale_fill() +
  geom_tile(
    data = dt_annotations,
    aes(x = offset_base - annotation_width, y = dt, fill = tech_platform),
    width = annotation_width,
    color = NA,
    linewidth = 0,
    height = 1,
    position = "identity"
  ) +
  scale_fill_manual(values = tech_colors, name = "Tech Platform") +
  guides(
    fill = guide_legend(
      nrow = 4,
      byrow = TRUE,
      order = 4  
    )
  ) +
  
  new_scale_fill() +
  geom_tile(
    data = dt_annotations,
    aes(x = offset_base - 2 * annotation_width, y = dt, fill = species),
    width = annotation_width,
    color = NA,
    linewidth = 0,
    height = 1,
    position = "identity"
  ) +
  scale_fill_manual(values = species_colors, name = "Species") +
  guides(
    fill = guide_legend(
      nrow = 3,
      byrow = TRUE,
      order = 5  
    )
  ) +
  
  new_scale_fill() +
  geom_tile(
    data = dt_annotations,
    aes(x = offset_base - 3 * annotation_width, y = dt, fill = tissue),
    width = annotation_width,
    color = NA,
    linewidth = 0,
    height = 1,
    position = "identity"
  ) +
  scale_fill_manual(values = tissue_colors, name = "Tissue") +
  guides(
    fill = guide_legend(
      ncol = 2,
      byrow = TRUE,
      order = 6  
    )
  ) +
  
  new_scale_fill() +
  geom_tile(
    data = pattern_annotations,
    aes(x = pattern_x, y = pattern_offset_base, fill = pattern),
    height = pattern_anno_height, color = NA
  ) +
  scale_fill_manual(values = pattern_colors, name = "Pattern") +
  
  geom_segment(
    data = na_df,
    aes(
      x = as.integer(method_pattern) - 0.4,
      xend = as.integer(method_pattern) + 0.4,
      y = as.integer(dt) - 0.4,
      yend = as.integer(dt) + 0.4
    ),
    colour = "#808080",
    linewidth = 0.25,
    alpha = 0.6
  ) +
  
  geom_segment(
    data = na_df,
    aes(
      x = as.integer(method_pattern) - 0.4,
      xend = as.integer(method_pattern) + 0.4,
      y = as.integer(dt) + 0.4,
      yend = as.integer(dt) - 0.4
    ),
    colour = "#808080",
    linewidth = 0.25,
    alpha = 0.6
  ) +
  
  scale_x_discrete(
    labels = function(x) {
      sub(".*:\\s*", "", x)
    },
    expand = c(0, 0)
  ) +
  scale_y_discrete(
    expand = c(0, 0)
  ) +
  theme_minimal() +
  my_theme +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    ),
    axis.title.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "right",
    legend.key.size = unit(0.06, "in"),
    legend.key.width = unit(0.06, "in"),
    legend.spacing.x = unit(0.01, "in"),
    legend.spacing.y = unit(0.01, "in"),
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 7),
    legend.box = "vertical",
    legend.margin = margin(1, 1, 0, 0),
    plot.margin = margin(3, 1, 0, 1),
    panel.grid = element_blank(),
    panel.border = element_blank()
  ) +
  labs(y = "Simulated datasets")

p1

p1_with_margin <- ggdraw() +
  draw_plot(p1, x = 1/6, y = 0, width = 5/6, height = 1)
ggsave('./Fig/Fig3-1.pdf',p1_with_margin , width = 6.69, height = 4, units = "in")

# ---- fig3C density ----
patterns <- c("pathology","hotspot", "stripe", "gradient",  "periodic", "neighbor")
paramset='P1'
datasets <- expand.grid(dt = dts_sp, pt = patterns, stringsAsFactors = FALSE) %>%
  mutate(dataset = sprintf("sim_%s-%s-%s-rep1", dt, pt, paramset))

## 1) 
pval_long <- do.call(rbind, lapply(seq_len(nrow(datasets)), function(i) {
  ds <- datasets$dataset[i]
  
  res <- get_pvalue_wide(ds, svg_id)   # wide: gene x method
  
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
plot_dendity <- pval_long %>%
  ggplot(aes(x = pval, color = neg_class, fill = neg_class)) +
  geom_density(alpha = 0.3) +
  labs(
    x = "p-value",
    y = "Density",
    color = "",
    fill = ""
  ) +
  theme_minimal()+
  facet_grid(~method)+
  ylim(0,10)+
  guides(
    color = guide_legend(
      label.position = "bottom", 
      label.hjust = 0.5,       
      title.position = "top",     
      title.hjust = 0.5           
    ),
    fill = guide_legend(
      label.position = "bottom", 
      label.hjust = 0.5,       
      title.position = "top",     
      title.hjust = 0.5           
    )    
  )+
  scale_x_continuous(
    labels = function(x) sprintf("%.1f", x)  #
  ) +  
  my_theme+
  theme(legend.key.size = unit(0.05, "in"), 
        legend.key.width = unit(0.05, "in"),
        legend.key.height = unit(0.05, "in") )+
  scale_fill_manual(values=c("#008bd0", "#ffa61d"))+
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
metric.df$methods <- factor(metric.df$methods,levels = c("C-SIDE","spVC", "Celina", "STANCE","CTSV","ctSVG"))
metric.df$pattern <- stringr::str_to_sentence(metric.df$pattern)

metric.df$pattern <- factor(metric.df$pattern,levels = c("Pathology","Stripe","Hotspot","Gradient","Periodic","Neighbor" ))

p4 <- metric.df %>%
  ggplot(aes(x = pattern, y = TPR, fill = methods)) +
  geom_boxplot(alpha = 1, width = 0.6, linewidth = 0.1, outlier.shape = NA, position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = method_colors) +
  my_theme +
  labs(fill = "", x = 'Pattern', y = "Sensitivity") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = "none",
        plot.margin = margin(0, 1.5, 0, 1),
        legend.key.size = unit(0.05, "in"),
        legend.key.width = unit(0.05, "in"),
          legend.direction = "horizontal",
    legend.box = "horizontal",     
    legend.margin = margin(1, 0, 0, 0))
p4
# ---- Save workspace and combine plots ----

p5 <- p5 + theme(plot.margin = margin(16, 1.5, 0, 5))

p45 <- plot_grid(p4, p5, nrow = 1, rel_widths = c(1.75, 1), labels = c("D", "E"))

# Add left margin to heatmap plot
p1_with_margin <- ggdraw() +
  draw_plot(p1, x = 1/6, y = 0, width = 5/6, height = 1)

# Combine plots into final figure
final_plot <- plot_grid(p1_with_margin, plot_dendity, p45,
                        nrow = 3,
                        rel_heights = c(1.7, 0.6,0.7),
                        labels = c("B", "C",""))
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

fpr_df <- rbind(fpr_df_sp,fpr_df_sc)
summary_df_spe_aff <- fpr_df %>%
  group_by(pattern, method) %>%
  summarise(specificity = (1-mean(fpr, na.rm = TRUE))) %>%
  ungroup()

spe_aff_wide <- summary_df_spe_aff %>%
  select(-starts_with("affected_specificity_")) %>%  
  pivot_wider(
    names_from   = pattern,
    values_from  = specificity,
    names_prefix = "affected_specificity_"
  )
spe_aff_wide <- as.data.frame(spe_aff_wide)
rownames(spe_aff_wide) <- spe_aff_wide$method
spe_aff_wide$method <- NULL
cols_to_add <- grep("^affected_specificity_", colnames(spe_aff_wide), value = TRUE)

summary_df[rownames(summary_df), cols_to_add] <- spe_aff_wide[rownames(summary_df), cols_to_add]

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

