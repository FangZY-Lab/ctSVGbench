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
source('./my_theme.R')

# List of datasets used in the analysis
datasets.sc <- c(
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
  "MERFISH_hypothalamus",
  "SeqFish+_cortex"  
)


##------------------
## Number of detected genes per method
##------------------
svnum.df.sc <- do.call(rbind, lapply(datasets.sc, function(dataset) {
  dat.pval.wide <- get_wide_pval(dataset,cell.level=T)
  dat.pval.binary <- as.data.frame((dat.pval.wide < 0.05) * 1)
  method_counts <- colSums(dat.pval.binary)
  df <- data.frame(
    method = names(method_counts),
    count = as.numeric(method_counts),
    dataset = dataset
  )
  return(df)
}))

datasets.sp <- c(
  "Slide-seq_tumor",
  "Slide-seqV2_hippocampus",
  "Slide-seqV2_mouseOB",
  "ST_developmental_heart",
  "ST_PDAC",
  "StereoSeq_MDESTA",
  "StereoSeq_mouseOB",
  "Visium_bladder",
  "Visium_intestine",
  "Visium_liver",
  "Visium_lymph_node",
  "Visium_mousebrain",
  "Visium_pancreas",
  "Visium_skin",
  "Visium_spleen",
  "Visium_tail",
  "Slide-seqV2_melanoma_GSM6025944_MBM13"
)


svnum.df.sp <- do.call(rbind, lapply(datasets.sp, function(dataset) {
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
svnum.df.sp$resolution <- "spot"
svnum.df.sc$resolution <- "cell"

svnum.df <- rbind(svnum.df.sc,svnum.df.sp)

dataset_order <- svnum.df %>%
  arrange(resolution,count) %>%
  pull(dataset) %>%
  unique()

# 3. 将dataset列转换为因子，使用新的排序
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
merged_df$method <- factor(merged_df$method, levels = 
                             c("spVC_1","spVC_2","Celina","STANCE","CTSV","C-SIDE","ctSVG"))

fill_colors <- setNames(merged_df$fill[match(dataset_order, merged_df$dataset)], dataset_order)
svnum.df <- svnum.df %>%
  mutate(dataset_order_number = as.integer(dataset))
write.csv(svnum.df,file = "Fig/fig2A.csv",row.names = F)
# Plot: detected gene count by dataset and method

p1 <- ggplot() +
  geom_rect(data = merged_df,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = dataset),
            color = NA, inherit.aes = FALSE) +
  scale_fill_manual(name = "Dataset", values = fill_colors,guide = "none") +
  geom_segment(data = merged_df, aes(x = x, y = 0, yend = count, linetype = resolution, color = method), size = 0.3) +
  geom_point(data = merged_df, aes(x = x, y = count, shape = resolution,color = method), size = 1) +
  facet_wrap(
    ~method,  
    ncol = 2, 
    labeller = label_value  
  ) +
  scale_color_manual(values = method_colors, guide = "none") +
  labs(x = "Datasets", y = "Detected gene count",shape="Dataset \nresolution",linetype="Dataset \nresolution") +
  theme_minimal() +
  my_theme +
  scale_y_sqrt() +
  theme(
    strip.background = element_blank(), 
    strip.placement = "outside", 
    legend.position = 'right',
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(2, 3, 0, 0),
    legend.key.size = unit(0.05, "in"),
    legend.key.width = unit(0.05, "in")
  )
p1
ggsave('./Fig/P1.pdf',width = 6.69,height = 3.5)

dataset_order.sc <- svnum.df.sc %>%
  arrange(resolution,count) %>%
  pull(dataset) %>%
  unique()

svnum.df.sp$dataset <- recode(svnum.df.sp$dataset,"Slide-seqV2_melanoma_GSM6025944_MBM13"="Slide-seqV2_melanoma_MBM13")
dataset_order.sp <- svnum.df.sp %>%
  arrange(resolution,count) %>%
  pull(dataset) %>%
  unique()

datasets <- c(datasets.sp,datasets.sc)

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
      overlap = factor(c("1 method", "2 methods", "3 methods", "4 methods","5 methods","6 methods"),
                       levels = c("1 method", "2 methods", "3 methods", "4 methods","5 methods","6 methods")),
      prop = 0,
      method = method
    ))
  }
  
  # Count overlap levels (only within selected method)
  tab <- table(factor(k[mask], levels = 1:6))
  prop <- as.numeric(tab) / denom
  
  tibble(
    overlap = factor(c("1 method", "2 methods", "3 methods", "4 methods","5 methods","6 methods"),
                     levels = c("1 method", "2 methods", "3 methods", "4 methods","5 methods","6 methods")),
    prop = prop,
    method = method
  )
}
df_plot_sp <- bind_rows(lapply(datasets.sp, function(ds) {
  pmat <- get_wide_pval(ds)
  bind_rows(lapply(colnames(pmat), function(m) {
    calc_overlap_by_method_partition(pmat, m, alpha) %>% mutate(dataset = ds)
  }))
}))
df_plot_sc <- bind_rows(lapply(datasets.sc, function(ds) {
  pmat <- get_wide_pval(ds,cell.level=T)
  bind_rows(lapply(colnames(pmat), function(m) {
    calc_overlap_by_method_partition(pmat, m, alpha) %>% mutate(dataset = ds)
  }))
}))
df_plot <- rbind(df_plot_sp,df_plot_sc)
library(ggplot2)
library(ggnewscale)

df_plot$dataset <- factor(df_plot$dataset, levels = c(dataset_order.sc, dataset_order.sp))
df_plot <- df_plot[df_plot$method != 'spVC_2',]
df_plot$method <- factor(df_plot$method, levels = c("spVC_1","spVC_2","Celina","STANCE","CTSV","C-SIDE","ctSVG"))

# ===== 底部两块色带的数据（连续两段）=====
n_sc <- length(dataset_order.sc)
n_sp <- length(dataset_order.sp)
x_min_sc <- 0.5
x_max_sc <- n_sc + 0.5
x_min_sp <- n_sc + 0.5
x_max_sp <- n_sc + n_sp + 0.5

band_df <- data.frame(
  xmin  = c(x_min_sc, x_min_sp),
  xmax  = c(x_max_sc, x_max_sp),
  ymin  = c(-0.015, -0.015),
  ymax  = c(-0.0, -0.0),
  level = c("Single-cell resolution", "Spot-level resolution")
)

# ===== 主图 =====
p2 <- ggplot(df_plot, aes(x = dataset, y = prop, fill = overlap)) +
  geom_col(position = "stack") +
  facet_wrap(~ method, nrow = 1) +
  scale_y_continuous(limits = c(-0.015, 1), expand = expansion(mult = c(0, 0.02))) +
  scale_fill_manual(
    values = c(
      "1 method"   = "#deebf7",
      "2 methods"  = "#c6dbef",
      "3 methods"  = "#9ecae1",
      "4 methods"  = "#6baed6",
      "5 methods"  = "#3182bd",
      "6 methods"  = "#08519c"
    ),
    name = "",
    guide = guide_legend(order = 1, nrow = 1)   # ✅ overlap 放左边
  ) +
  ggnewscale::new_scale_fill() +
  geom_rect(
    data = band_df,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = level),
    inherit.aes = FALSE,
    color = NA
  ) +
  scale_fill_manual(
    values = c(
      "Single-cell resolution" = "#FFD6A5",
      "Spot-level resolution"  = "#A8E6CF"
    ),
    name = "Data resolution",
    guide = "none" 
    # guide = guide_legend(order = 2, nrow = 1)   # ✅ resolution 放右边
  ) +
  labs(x = "Datasets", y = "Proportion (within method)", fill = "") +
  theme_minimal() +
  my_theme +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(2, 0, 0, 0),
    legend.key.size = unit(0.02, "in"),
    legend.key.width = unit(0.02, "in"),
    legend.spacing.x = unit(0.02, "in"),  
    legend.margin = margin(0, 0, 0, 0),   
    legend.position = "bottom",
    panel.grid        = element_blank(),
    panel.border      = element_rect(color = "black", linewidth = 0.2, fill = NA),
    
    # plot.margin = margin(t = 5.5, r = 5.5, b = 25, l = 5.5)
  ) +
  coord_cartesian(clip = "off")

p2

##------------------
## Rank correlation between methods (Top 100 genes)
##------------------

library(dplyr)

get_top100_rank_corr <- function(dataset,cell.level=F) {
  # Load the p-value matrix (genes × methods)
  dat.pval.wide <- get_wide_pval(dataset,cell.level=cell.level)
  
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
res.list.sc <- lapply(datasets.sc, get_top100_rank_corr, cell.level = TRUE)
long_cor.sc <- do.call(rbind, res.list.sc)
long_cor.sc$method1 <- factor(long_cor.sc$method1,levels = 
                                c("spVC_1","spVC_2","Celina","STANCE","CTSV","C-SIDE","ctSVG"))
long_cor.sc$method2 <- factor(long_cor.sc$method2,levels = 
                                c("spVC_1","spVC_2","Celina","STANCE","CTSV","C-SIDE","ctSVG"))

long_cor.sc$dataset_short <- case_when(
  long_cor.sc$dataset == "StereoSeq_CBMSTA_Macaque1_T110" ~ "cbmsta_Macaque1_T110",
  long_cor.sc$dataset == "StereoSeq_CBMSTA_Macaque1_T42" ~ "cbmsta_Macaque1_T42",
  long_cor.sc$dataset == "StereoSeq_CBMSTA_Marmoset1_T478" ~ "cbmsta_Marm1_T478",
  long_cor.sc$dataset == "StereoSeq_CBMSTA_Marmoset1_T514" ~ "cbmsta_Marm1_T514",
  long_cor.sc$dataset == "StereoSeq_CBMSTA_Mouse1_T167" ~ "cbmsta_Mouse1_T167",
  long_cor.sc$dataset == "StereoSeq_CBMSTA_Mouse1_T169" ~ "cbmsta_Mouse1_T169",
  long_cor.sc$dataset == "StereoSeq_CBMSTA_Mouse1_T171" ~ "cbmsta_Mouse1_T171",
  long_cor.sc$dataset == "StereoSeq_CBMSTA_Mouse1_T176" ~ "cbmsta_Mouse1_T176",
  long_cor.sc$dataset == "StereoSeq_CBMSTA_Mouse1_T185" ~ "cbmsta_Mouse1_T185",
  long_cor.sc$dataset == "StereoSeq_CBMSTA_Mouse1_T189" ~ "cbmsta_Mouse1_T189",
  long_cor.sc$dataset ==  "MERFISH_hypothalamus" ~ "MERFISH_Hypotha",   
  long_cor.sc$dataset == "VisiumHD_LUAD_2431" ~ "visiumHD_2431",
  long_cor.sc$dataset == "VisiumHD_LUAD_6123" ~ "visiumHD_6123",
  long_cor.sc$dataset == "VisiumHD_LUAD_6976" ~ "visiumHD_6976",
  long_cor.sc$dataset == "VisiumHD_LUSC_5488" ~ "visiumHD_5488",
  long_cor.sc$dataset == "VisiumHD_LUSC_7437" ~ "visiumHD_7437",
  long_cor.sc$dataset == "VisiumHD_LUSC_7941" ~ "visiumHD_7941",
  long_cor.sc$dataset == "SeqFish+_cortex" ~ "seqfish+_cortex",
  TRUE ~ str_sub(long_cor.sc$dataset, 1, 10) 
)
# Heatmap: correlation matrix for each dataset
p <- ggplot(long_cor.sc, aes(x = method1, y = method2, fill = correlation)) +
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
  facet_wrap(~ dataset_short, ncol = 4)
p
ggsave(sprintf("./Fig/s/S_real_cor2_scell.pdf"), width = 6.69, height = 9)

res.list.sp <- lapply(datasets.sp, get_top100_rank_corr)
long_cor.sp <- do.call(rbind, res.list.sp)
long_cor.sp$method1 <- factor(long_cor.sp$method1,levels = 
                                c("spVC_1","spVC_2","Celina","STANCE","CTSV","C-SIDE","ctSVG"))
long_cor.sp$method2 <- factor(long_cor.sp$method2,levels = 
                                c("spVC_1","spVC_2","Celina","STANCE","CTSV","C-SIDE","ctSVG"))

long_cor.sp$dataset_short <- case_when(
  long_cor.sp$dataset == "Slide-seq_tumor" ~ "slideseq_tumor",
  long_cor.sp$dataset == "Slide-seqV2_hippocampus" ~ "slideseqV2_hip",
  long_cor.sp$dataset == "Slide-seqV2_mouseOB" ~ "slideseqV2_OB",
  long_cor.sp$dataset == "ST_developmental_heart" ~ "ST_dev_heart",
  long_cor.sp$dataset == "ST_PDAC" ~ "ST_PDAC",
  long_cor.sp$dataset == "StereoSeq_MDESTA" ~ "cbmsta_MDESTA",
  long_cor.sp$dataset == "StereoSeq_mouseOB" ~ "cbmsta_mouseOB",
  long_cor.sp$dataset == "Visium_bladder" ~ "visium_bladder",
  long_cor.sp$dataset == "Visium_intestine" ~ "visium_intes",
  long_cor.sp$dataset == "Visium_liver" ~ "visium_liver",
  long_cor.sp$dataset == "Visium_lymph_node" ~ "visium_lymph",
  long_cor.sp$dataset == "Visium_mousebrain" ~ "visium_brain",
  long_cor.sp$dataset == "Visium_pancreas" ~ "visium_panc",
  long_cor.sp$dataset == "Visium_skin" ~ "visium_skin",
  long_cor.sp$dataset == "Visium_spleen" ~ "visium_spleen",
  long_cor.sp$dataset == "Visium_tail" ~ "visium_tail",
  long_cor.sp$dataset == "Slide-seqV2_melanoma_GSM6025944_MBM13" ~ "slideseqV2_MBM13",
  TRUE ~ str_sub(long_cor.sp$dataset, 1, 10)  
)
# Heatmap: correlation matrix for each dataset
p <- ggplot(long_cor.sp, aes(x = method1, y = method2, fill = correlation)) +
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
  facet_wrap(~ dataset_short, ncol = 4)
p
ggsave(sprintf("./Fig/s/S_real_cor2_spot.pdf"), width = 6.69, height = 9)


long_cor <- rbind(long_cor.sc,long_cor.sp)

# Boxplot: pairwise correlation between methods
as.data.frame(long_cor) %>% 
  filter(method1 != method2) %>%
  ggplot(aes(x = method1, y = correlation, fill = method2)) +
  geom_boxplot(alpha = 1, width = 0.6, color = "grey30") +
  scale_fill_manual(values = method_colors) +
  theme_test() +
  my_theme +
  theme(
    # plot.margin = margin(2, 2, 0, 0),
    axis.title.y = element_text(size = 7),
    legend.position = 'bottom'
  ) +
  guides(fill = guide_legend(nrow = 1))+
  labs(x = "", y = "Rank correlation", fill = "")
ggsave(sprintf("./Fig/s/S_real_cor1.pdf"), width = 6.69, height = 3)



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
jaccard.df.sc <- do.call(rbind, lapply(datasets.sc, function(dataset) {
  dat.pval.wide <- get_wide_pval(dataset,cell.level=T)
  
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

jaccard.df.sp <- do.call(rbind, lapply(datasets.sp, function(dataset) {
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
jaccard.df <- rbind(jaccard.df.sc,jaccard.df.sp)
# Remove self-comparisons
jaccard.df <- subset(jaccard.df, method1 != method2)
jaccard.df <- rename(jaccard.df, method = method1)
jaccard.df$method <- factor(jaccard.df$method,levels = 
                              c("spVC_1","spVC_2","Celina","STANCE","CTSV","C-SIDE","ctSVG"))
jaccard.df$method2 <- factor(jaccard.df$method2,levels = 
                               c("spVC_1","spVC_2","Celina","STANCE","CTSV","C-SIDE","ctSVG"))



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
    axis.text.x = element_text(size = 7,angle = 45,hjust = 1),
    legend.position = "bottom",
    plot.margin = margin(2, 3, 0, 2.5),
    legend.key.size = unit(0.01, "in"),
    legend.key.width = unit(0.03, "in"),
    legend.spacing.x = unit(0.02, "in"),    
  ) +
  guides(fill = guide_legend(nrow = 1))

p5

dts <- c( 
  "StereoSeq_CBMSTA_Mouse1_T189",
  "StereoSeq_CBMSTA_Mouse2_T349",
  "MERFISH_hypothalamus")

conc.res1 <- do.call(rbind, lapply(dts, function(dataset) {
  dat.pval.wide <- get_wide_pval(dataset,cell.level=T)
  
  data <- get_conc(dat.pval.wide, dataset = dataset)
  return(data)
}))
conc.res1=conc.res1[conc.res1$method1 != 'spVC_2',]

conc.res1$method1 <- factor(conc.res1$method1,levels = 
                              c("spVC_1","spVC_2","Celina","STANCE","CTSV","C-SIDE","ctSVG"))


conc.res1 <- conc.res1 %>%
  mutate(
    dataset_short = case_when(
      dataset == "MERFISH_hypothalamus" ~ "merfish-Hypo",  
      dataset == "StereoSeq_CBMSTA_Mouse1_T189" ~ "cbmsta_T189",
      dataset == "StereoSeq_CBMSTA_Mouse2_T349" ~ "cbmsta_T349",
      TRUE ~ str_sub(dataset, 1, 10)  
    )
  )

p3 <- ggplot(conc.res1, aes(x = rank, y = conc, color = method2)) +
  geom_line(linewidth = 0.3) +
  facet_grid(dataset_short ~ method1) + 
  theme_minimal() +
  my_theme +
  scale_color_manual(values = method_colors) +
  labs(title = "", y = "Concordance at the top", x = 'Rank', color = "") +
  theme(
    strip.text = element_text(color = "black", size = 7),
    axis.title.y = element_text(size = 7),
    axis.title.x = element_text(size = 7),
    legend.position = "none",
    plot.margin = margin(0, 0, 0, 0)
  )+
  guides(color = guide_legend(nrow = 1))
p3

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
  rel_heights = c(1.5, 1, 1.1),
  nrow = 3,
  labels = c("A", "", "D"),
  label_x = -0.005,
  label_y = 1.00
)
ggsave("./Fig/Fig2.pdf", width = 6.69, height = 9)


##------------------
## Supplemental concordance plots (paginated)
##------------------

conc.res.sc <- do.call(rbind, lapply(datasets.sc, function(dataset) {
  dat.pval.wide <- get_wide_pval(dataset,cell.level=T)
  data <- get_conc(dat.pval.wide, dataset = dataset)
  return(data)
}))

conc.res.sp <- do.call(rbind, lapply(datasets.sp, function(dataset) {
  dat.pval.wide <- get_wide_pval(dataset)
  data <- get_conc(dat.pval.wide, dataset = dataset)
  return(data)
}))

conc.res.sc$method1 <- factor(conc.res.sc$method1,levels = 
                                c("spVC_1","spVC_2","Celina","STANCE","CTSV","C-SIDE","ctSVG"))
conc.res.sc <- subset(conc.res.sc, !(dataset %in% dts))

conc.res.sp$method1 <- factor(conc.res.sp$method1,levels = 
                                c("spVC_1","spVC_2","Celina","STANCE","CTSV","C-SIDE","ctSVG"))

conc.res.sc$dataset_short <- case_when(
  conc.res.sc$dataset == "StereoSeq_CBMSTA_Macaque1_T110" ~ "cbmsta_Macaq1_T110",
  conc.res.sc$dataset == "StereoSeq_CBMSTA_Macaque1_T42" ~ "cbmsta_Macaq1_T42",
  conc.res.sc$dataset == "StereoSeq_CBMSTA_Marmoset1_T478" ~ "cbmsta_Marm1_T478",
  conc.res.sc$dataset == "StereoSeq_CBMSTA_Marmoset1_T514" ~ "cbmsta_Marm1_T514",
  conc.res.sc$dataset == "StereoSeq_CBMSTA_Mouse1_T167" ~ "cbmsta_Mouse1_T167",
  conc.res.sc$dataset == "StereoSeq_CBMSTA_Mouse1_T169" ~ "cbmsta_Mouse1_T169",
  conc.res.sc$dataset == "StereoSeq_CBMSTA_Mouse1_T171" ~ "cbmsta_Mouse1_T171",
  conc.res.sc$dataset == "StereoSeq_CBMSTA_Mouse1_T176" ~ "cbmsta_Mouse1_T176",
  conc.res.sc$dataset == "StereoSeq_CBMSTA_Mouse1_T185" ~ "cbmsta_Mouse1_T185",
  conc.res.sc$dataset == "StereoSeq_CBMSTA_Mouse1_T189" ~ "cbmsta_Mouse1_T189",
  conc.res.sc$dataset ==  "MERFISH_hypothalamus" ~ "MERFISH_Hypotha",
  conc.res.sc$dataset == "VisiumHD_LUAD_2431" ~ "visiumHD_2431",
  conc.res.sc$dataset == "VisiumHD_LUAD_6123" ~ "visiumHD_6123",
  conc.res.sc$dataset == "VisiumHD_LUAD_6976" ~ "visiumHD_6976",
  conc.res.sc$dataset == "VisiumHD_LUSC_5488" ~ "visiumHD_5488",
  conc.res.sc$dataset == "VisiumHD_LUSC_7437" ~ "visiumHD_7437",
  conc.res.sc$dataset == "VisiumHD_LUSC_7941" ~ "visiumHD_7941",
  conc.res.sc$dataset == "SeqFish+_cortex" ~ "seqfish+_cortex",
  TRUE ~ str_sub(conc.res.sc$dataset, 1, 10) 
)

datasets <- rev(c(
  "StereoSeq_CBMSTA_Macaque1_T110",
  "StereoSeq_CBMSTA_Mouse1_T176",
  "VisiumHD_LUSC_7437",
  "StereoSeq_CBMSTA_Marmoset1_T478",
  "VisiumHD_LUSC_5488",
  "VisiumHD_LUSC_7941",
  "SeqFish+_cortex",
  "StereoSeq_CBMSTA_Macaque1_T42",
  "StereoSeq_CBMSTA_Marmoset1_T514",
  "StereoSeq_CBMSTA_Mouse1_T167",
  "StereoSeq_CBMSTA_Mouse1_T169",
  "StereoSeq_CBMSTA_Mouse1_T171",
  "StereoSeq_CBMSTA_Mouse1_T185",
  "VisiumHD_LUAD_2431",
  "VisiumHD_LUAD_6123",
  "VisiumHD_LUAD_6976"
))
library(dplyr)

conc.res.sc <- conc.res.sc %>%
  mutate(dataset = factor(dataset, levels = datasets)) %>%
  arrange(dataset)


i <- 1
start <- 1
while (i <= 3) {
  pages <- ifelse(i == 3, 6, 5)
  subset_data <- conc.res.sc %>% 
    subset(dataset %in% unique(conc.res.sc$dataset)[start:(start + pages - 1)])%>% 
    subset(method1 !="spVC")
  
  pS3 <- ggplot(subset_data, aes(x = rank, y = conc, color = method2)) +
    geom_line() +
    facet_wrap(dataset_short~method1, ncol = 6, nrow = pages) +
    theme_minimal() +
    my_theme +
    scale_color_manual(values = method_colors) +
    labs(title = "", y = "Concordance at the top", color = "") +
    guides(color = guide_legend(override.aes = list(size = 4))) +
    theme(
      strip.text = element_text(color = "black", size = 7),
      legend.position = "bottom"
    )+
    guides(color = guide_legend(nrow = 1))
  
  
  ggsave(sprintf("./Fig/s/conc_methods_sc%s.pdf", i), width = 6.9, height = 8)
  
  start <- start + pages
  i <- i + 1
}


conc.res.sp$dataset_short <- case_when(
  conc.res.sp$dataset == "Slide-seq_tumor" ~ "slideseq_tumor",
  conc.res.sp$dataset == "Slide-seqV2_hippocampus" ~ "slideseqV2_hip",
  conc.res.sp$dataset == "Slide-seqV2_mouseOB" ~ "slideseqV2_OB",
  conc.res.sp$dataset == "ST_developmental_heart" ~ "ST_dev_heart",
  conc.res.sp$dataset == "ST_PDAC" ~ "ST_PDAC",
  conc.res.sp$dataset == "StereoSeq_MDESTA" ~ "cbmsta_MDESTA",
  conc.res.sp$dataset == "StereoSeq_mouseOB" ~ "cbmsta_mouseOB",
  conc.res.sp$dataset == "Visium_bladder" ~ "visium_bladder",
  conc.res.sp$dataset == "Visium_intestine" ~ "visium_intes",
  conc.res.sp$dataset == "Visium_liver" ~ "visium_liver",
  conc.res.sp$dataset == "Visium_lymph_node" ~ "visium_lymph",
  conc.res.sp$dataset == "Visium_mousebrain" ~ "visium_brain",
  conc.res.sp$dataset == "Visium_pancreas" ~ "visium_panc",
  conc.res.sp$dataset == "Visium_skin" ~ "visium_skin",
  conc.res.sp$dataset == "Visium_spleen" ~ "visium_spleen",
  conc.res.sp$dataset == "Visium_tail" ~ "visium_tail",
  conc.res.sp$dataset == "Slide-seqV2_melanoma_GSM6025944_MBM13" ~ "slideseqV2_MBM13",
  TRUE ~ str_sub(conc.res.sp$dataset, 1, 10)  
)

i <- 1
start <- 1
while (i <= 3) {
  pages <- ifelse(i == 1, 5, 6)
  
  subset_data <- conc.res.sp %>% 
    subset(dataset %in% unique(conc.res.sp$dataset)[start:(start + pages - 1)])%>% 
    subset(method1 !="spVC")
  
  pS3 <- ggplot(subset_data, aes(x = rank, y = conc, color = method2)) +
    geom_line() +
    facet_wrap(dataset_short~method1 , ncol = 5, nrow = pages) +
    theme_minimal() +
    my_theme +
    scale_color_manual(values = method_colors) +
    labs(title = "", y = "Concordance at the top", color = "") +
    guides(color = guide_legend(override.aes = list(size = 4))) +
    theme(
      strip.text = element_text(color = "black", size = 7),
      legend.position = "bottom"
    )+
    guides(color = guide_legend(nrow = 1))
  
  
  ggsave(sprintf("./Fig/s/conc_methods_sp%s.pdf", i), width = 6.9, height = 8)
  
  start <- start + pages
  i <- i + 1
}


