library(RhpcBLASctl)
blas_set_num_threads(1)
library(spacexr)
library(ggplot2)
library(scales)
library(tidyr)
library(dplyr)
library(tibble)
library(purrr)
library(RColorBrewer)
library(stringr)
library(tidyverse)
library(data.table)
library(cowplot)
library(ggsci)
library(ggpubr)
library(ggtext)
library(ggprism)
library(ggpmisc)

# Specialized and utility packages
library(SpatialExperiment)
library(Matrix)
library(devtools)
library(reshape2)
library(here)
library(Seurat)
library(patchwork)

datasets <- c(
  "VisiumHD_LUAD_2431", 
  "VisiumHD_LUAD_6123", 
  "VisiumHD_LUAD_6976", 
  "VisiumHD_LUSC_5488", 
  "VisiumHD_LUSC_7437", 
  "VisiumHD_LUSC_7941"
)
source('./my_theme.R')

# # Step 1: Spatial Cell Type Mapping Visualization
# plots <- list()
# for(dataset in datasets){
#   data_dir <- paste0('./st_lusc+luad/', substr(dataset,15,25), "/square_008um")
#   file <- sprintf('myRCTD_%s.rds', dataset)
#   prop <- readRDS(here('real','prop', file))
  
#   # Load raw Visium data and create Seurat object
#   counts <- Read10X(file.path(data_dir, "filtered_feature_bc_matrix"))
#   obj <- CreateSeuratObject(counts, assay = "Spatial")
#   image.data <- Seurat:::Read10X_Image(
#     image.dir = file.path(data_dir, "spatial"),
#     filter.matrix = TRUE
#   )
#   obj[["slice1"]] <- image.data  
#   obj <- NormalizeData(obj, assay = "Spatial")
  
#   # Assign cell types
#   celltype_assign <- colnames(prop)[max.col(prop, ties.method = "first")]
#   names(celltype_assign) <- rownames(prop)
#   spots_to_plot <- names(celltype_assign)
  
#   # Subset Seurat object
#   obj_subset <- obj[, spots_to_plot]
#   obj_subset <- AddMetaData(obj_subset, metadata = celltype_assign, col.name = "celltype")
  
#   # Define cell type colors
#   cols <- RColorBrewer::brewer.pal(5, "Set2")
#   names(cols) <- c("Cancer", "Endothelial", "Epithelial", "Immune", "Stroma")
#   cols <- c(cols, "Unassigned" = NA)
  
#   # Plot spatial map
#   p <- SpatialDimPlot(
#     obj_subset,
#     group.by = "celltype",
#     pt.size.factor = 15,
#     label = FALSE,
#     crop = TRUE 
#   ) +
#     scale_fill_manual(values = cols, na.value = "lightgrey") +
#     theme(plot.title = element_text(hjust = 0.5, size = 7),
#           legend.position = "bottom",
#           text = element_text(size = 7),
#           axis.title = element_text(size = 7),
#           axis.text = element_text(size = 7),
#           legend.key.size = unit(0.2, "in")) +
#     ggtitle(substr(dataset,10,40)) +
#     labs(fill="Cell type\n(NA = grey)")
  
#   plots[[dataset]] <- p  
# }

# # Combine all tissue plots
# combined_tissue <- wrap_plots(plots, ncol = 3) + plot_layout(guides = "collect") &
#   theme(legend.position = "bottom",
#         legend.key.width = unit(0.1, "in"),
#         legend.text = element_text(size = 7))


# Step 2: Identify Significant Genes per Method

sig <- list()
for(dataset in datasets){
  methods <- c('C-SIDE','spVC','CELINA','STANCE',"CTSV","ctsvg")
  file <- sprintf('myRCTD_%s.rds', dataset)
  puck <- readRDS(here('real','puck', sprintf('myRCTD_%s.rds', dataset)))
  ctsvg=readRDS(here('real','res',sprintf('%s-ctsvg.rds',dataset)))
  res.ctsvg <- split(ctsvg, ctsvg$cluster)
  res.ctsvg <- lapply(res.ctsvg, \(df)
                      data.frame(
                        pval = df$pval,
                        row.names = df$gene
                      ))  
  res.ctsv=readRDS(here('real','res',sprintf('%s-CTSV.rds',dataset)))                       
  res.cside <- readRDS(here('real','res', sprintf('%s-C-SIDE.rds', dataset)))
  res.celina <- readRDS(here('real','res', sprintf('%s-CELINA.rds', dataset)))
  res.stance <- readRDS(here('real','res', sprintf('%s-STANCE.rds', dataset)))
  
  prop <- readRDS(here('real','prop', sprintf('myRCTD_%s.rds', dataset))) 
  pos <- readRDS(here('real','pos', sprintf('myRCTD_%s.rds', dataset))) 
  spVC=readRDS(here('real','res',sprintf('%s-spVC.rds',dataset)))
  if(is.null(spVC)){
    res.spVC=NULL
  }else {
    idx=match(names(res.celina),colnames(prop))
    genes.v=names(spVC$results.varying)
    res.spVC <- lapply(idx,function(ct){
      pval=sapply(spVC$results.varying[genes.v],function(x){
        x$p.value[paste0("gamma_X", ct)]
      })
      names(pval)=sapply(strsplit(names(pval),"\\."),"[[",1)
      data.frame(pval = na.omit(pval))
    })
    
    names(res.spVC) <- names(res.celina)
  }
  
  # Combine all method results
  all_lists <- list(CSIDE = res.cside,
                    spVC = if (is.null(spVC)) NULL else res.spVC,
                    Celina = res.celina, 
                    STANCE = res.stance,
                    ctSVG = res.ctsvg,
                    CTSV = res.ctsv)
  all_genes <- unique(unlist(lapply(all_lists, function(lst) {
    unlist(lapply(lst, rownames))
  })))
  
  # Construct a unified p-value data frame
  dat.pval <- do.call(rbind, lapply(names(all_lists), function(method){
    lst <- all_lists[[method]]
    n <- min(3, length(names(lst)))
    list <- lapply(names(lst)[1:n], function(i){  
      df <- lst[[i]]
      if(nrow(df) > 0){ 
        pval_cols <- grep("p[_\\.]?val|pvalue|pval", colnames(df), ignore.case = TRUE) 
        if (length(pval_cols) != 1) stop(sprintf("Data frame %s[[%d]] missing pval column", method, i))
        pvals <- p.adjust(df[,pval_cols], method = "BH")
        names(pvals) <- rownames(df)
        pvals_full <- rep(1, length(all_genes))
        names(pvals_full) <- all_genes
        pvals_full[names(pvals)] <- pvals
        pvals_full <- data.frame(pval = pvals_full)
        pvals_full$gene <- paste0(make.names(i), '_', rownames(pvals_full))
        return(pvals_full)
      } else {
        pvals_full <- rep(1, length(all_genes))
        names(pvals_full) <- all_genes
        pvals_full <- data.frame(pval = pvals_full)
        pvals_full$gene <- paste0(make.names(i), '_', rownames(pvals_full))
        return(pvals_full)
      }
    })
    list <- list[!sapply(list, is.null)]
    if(length(list) > 0){
      df.pval <- do.call(rbind, list)
      df.pval$method <- method
      return(df.pval)
    } else return(NULL)
  }))
  
  dat.pval$method <- ifelse(dat.pval$method == "CSIDE", "C-SIDE", dat.pval$method)
  alpha <- 0.05
  
  dat.pval <- dat.pval %>%
    mutate(gene = sub("\\.", " ", gene)) %>%
    separate(gene, into = c("cell_type", "gene_name"), sep = "_")
  df.sig <- dat.pval %>% filter(pval < alpha)
  sig[[dataset]] <- df.sig
}
# write.csv(sig_df,'Fig/TableS3_lung_ctSVG_gene.csv',row.names = F)


# Step 3: Compute Jaccard Similarity Across Patients
sig_df <- bind_rows(sig, .id = "sample_id") %>%
  separate(sample_id, into = c("platform", "cancer_sample"), sep = "_", extra = "merge", remove = FALSE) %>%
  separate(cancer_sample, into = c("cancer", "sample"), sep = "_", extra = "merge", remove = FALSE)

# Define color palette
cols <- RColorBrewer::brewer.pal(5, "Set2")
names(cols) <- c("Cancer", "Endothelial", "Epithelial", "Immune", "Stroma")

# Compute Jaccard index between patients for each method and cancer type
gene_sets2 <- sig_df %>%
  group_by(cancer, cell_type, cancer_sample, method) %>%
  summarise(genes = list(unique(gene_name)), .groups = "drop") %>% 
  subset(cell_type == "Stroma")

jaccard_index <- function(a, b) length(intersect(a, b)) / length(union(a, b))

jar_results <- gene_sets2 %>%
  group_by(cancer, method, cell_type) %>%
  filter(n_distinct(cancer_sample) >= 2) %>%
  reframe({
    samples <- unique(cancer_sample)
    pairs <- combn(samples, 2, simplify = FALSE)
    map_dfr(pairs, function(p) {
      genes1 <- genes[cancer_sample == p[1]][[1]]
      genes2 <- genes[cancer_sample == p[2]][[1]]
      tibble(
        sample1 = p[1],
        sample2 = p[2],
        jaccard = jaccard_index(genes1, genes2)
      )
    })
  }) %>%
  ungroup()


summary_df <- jar_results %>%
  group_by(method, cancer, cell_type) %>%
  summarise(mean_jaccard = mean(jaccard), sd_jaccard = sd(jaccard), .groups = "drop")
method_levels <- c("C-SIDE", "spVC", "Celina", "STANCE","CTSV","ctSVG")
summary_df$method <- factor(summary_df$method, levels = method_levels)
jar_results$method  <- factor(jar_results$method,  levels = method_levels)

jac_plot <- ggplot(summary_df, aes(x = cancer, y = mean_jaccard, fill = cancer)) +
  geom_col(alpha = 1) +
  geom_errorbar(aes(ymin = mean_jaccard - sd_jaccard, ymax = mean_jaccard + sd_jaccard), width = 0.1) +
  geom_jitter(data = jar_results, aes(x = cancer, y = jaccard), width = 0.13, size = 0.1, color = "black") +
  facet_wrap(~method, ncol = 5) +
  theme_minimal(base_size = 7) +
  my_theme+
  scale_fill_manual(values = c("#33A02C", "#1F78B4")) +
  labs(y = "Mean Jaccard Index", x = "") +
  theme(legend.position = "none",
   axis.text.x = element_text(angle = 90, hjust = 0.5),
    plot.margin = margin(0, 2, 0, 1),
    panel.spacing.x = unit(0.01, "in"))
jac_plot
ggsave("Fig/jac.pdf",width=1)

# Step 4: Expression vs. Significance Correlation

source('./real/utils/real_expr_bench.R')

df <- do.call(rbind, lapply(datasets, function(dataset){
  dat.pval.wide <- get_pval_Stroma(dataset)
  file <- sprintf('myRCTD_%s.rds', dataset)
  puck <- readRDS(here('real','puck', file))
  
  counts <- puck@counts
  counts <- counts[rowSums(counts) > 0, colSums(counts) > 0]
  expr_normalized <- scater::normalizeCounts(counts)
  gene_expr <- Matrix::rowMeans(expr_normalized)
  
  expr_df <- data.frame(gene = names(gene_expr), expr = gene_expr)
  dat_long <- dat.pval.wide %>%
    rownames_to_column("gene") %>%
    pivot_longer(-gene, names_to = "method", values_to = "pval") %>%
    mutate(logp = -log10(pval + 1e-300))
  
  plot_df <- left_join(dat_long, expr_df, by = "gene")
  plot_df$dataset <- dataset
  return(plot_df)
}))
df$dataset <- substr(df$dataset, 10, 30)
df <- df %>%
  mutate(cancer = substr(dataset, 1, 4))

df$method <- factor(df$method,levels = c("C-SIDE","spVC", "Celina", "STANCE","CTSV","ctSVG")) 

df_stat_LUSC <- df %>%
  filter(cancer == "LUSC") %>%
  group_by(dataset, method) %>%
  filter(
    sd(expr, na.rm = TRUE) > 0,
    sd(logp, na.rm = TRUE) > 0
  ) %>%
  ungroup()


lusc_cor <- ggplot(
  data = subset(df, cancer == "LUSC"),
  aes(x = log10(expr + 0.5), y = logp)
) +
  geom_point(alpha = 0.4) +

  facet_grid(dataset ~ method) +

  stat_cor(
    data = df_stat_LUSC,         
    method = "spearman",
    label.x.npc = "left",
    label.y.npc = "top",
    size = 1.5
  ) +

  labs(
    x = expression(log[10]("normalized mean expression")),
    y = expression(-log[10](pval))
  ) +
  theme_minimal() +
  my_theme +
  theme(plot.margin = margin(0, 2, 0, 2),
    axis.text.x = element_text(size=4))




df_stat_LUAD <- df %>%
  filter(cancer == "LUAD") %>%
  group_by(dataset, method) %>%
  filter(
    sd(expr, na.rm = TRUE) > 0,
    sd(logp, na.rm = TRUE) > 0
  ) %>%
  ungroup()


luad_cor <- ggplot(
  data = subset(df, cancer == "LUAD"),
  aes(x = log10(expr + 0.5), y = logp)
) +
  geom_point(alpha = 0.4) +
  
  facet_grid(dataset ~ method) +
  
  stat_cor(
    data = df_stat_LUAD,         
    method = "spearman",
    label.x.npc = "left",
    label.y.npc = "top",
    size = 1.5
  ) +
  
  labs(
    x = expression(log[10]("normalized mean expression")),
    y = expression(-log[10](pval))
  ) +
  theme_minimal() +
  my_theme +
  theme(plot.margin = margin(0, 2, 0, 2),
  axis.text.x = element_text(size=5))

ggsave('./Fig/s/luad_cor.pdf', width = 6.69, height = 3.69)

# Compute Spearman correlations
cor_res <- df %>%
  group_by(method, dataset,cancer) %>%
  summarise(cor = cor(expr, logp, method = "spearman", use = "complete.obs"))
cor_res <- cor_res %>% filter(!is.na(cor))


# Violin plot summarizing correlations
plot_vl <- ggplot(cor_res, aes(x = method, y = cor)) +
  geom_violin(aes(fill = method),alpha = 0.6,show.legend = FALSE) +
  geom_jitter(
    aes(color = cancer),   
    width = 0.2,
    size = 2,
    alpha = 0.8
  )+
  theme_bw() +
  scale_color_manual(values = c("#33A02C", "#1F78B4")) +
  scale_fill_manual(values = method_colors)+
  labs(x = "Methods", y = "Spearman correlation", color = "Cancer") +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.margin = margin(0, 3, 0, 0)

  )+ 
  my_theme

plot_vl

# Step 5: Combine and Save Final Figures
p1 <- ggdraw() +
  draw_plot(jac_plot, x = 4/6.29, y = 0, width = 2.29/6.29, height = 1)

p2 <- plot_grid(lusc_cor,plot_vl, ncol = 2,rel_widths = c(1,0.4),
                labels = c('C',"D"),label_x = -0.01,label_y = 1.01)

final_plot <- plot_grid(p1, p2, nrow = 2, rel_heights = c(3, 3.7))
ggsave('./Fig/Fig6.1.pdf', width = 6.69, height = 6.7)


blank_plot <- ggplot() + 
  theme_void() 

p1 <- ggdraw() +
  draw_plot(jac_plot, x = 4/6.29, y = 0, width = 2.29/6.29, height = 1)

p2 <- plot_grid(
  blank_plot,  
  plot_vl, 
  ncol = 2,
  rel_widths = c(1,0.4),
  labels = c('C',"D"),  
  label_x = -0.01,
  label_y = 1.01
)
