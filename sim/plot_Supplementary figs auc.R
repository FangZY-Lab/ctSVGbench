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
paramset = 'P2'
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


ggsave(sprintf('./Fig/s/sim-auc-%s.pdf',paramset),p1, width = 6.69, height = 4, units = "in")
