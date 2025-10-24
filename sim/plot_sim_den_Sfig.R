# Required packages
library(reshape2)
library(ggplot2)
library(dplyr)
library(patchwork)
library(tidyr)
source('F:/ctSVGbench/my_theme.R')
# List of datasets
dts <- c(
  "ST_PDAC",
  # "Visium_liver", # excluded due to lack of celltype4
  "Visium_mousebrain",
  "StereoSeq_MDESTA",
  "StereoSeq_CBMSTA_Macaque",
  "Visium_spleen",
  "Slide-seqV2_melanoma",
  "SeqFish+_mouse_ob",
  "StereoSeq_CBMSTA_Marmoset",
  "Slide-seq_tumor"
)

# Define spatial patterns
patterns <- c("pathology", "hotspot", "stripe", "gradient", "periodic", "neighbor")

# Source utility functions for simulation benchmarking
source('./sim/utils/sim-bench.R')

# ---- Data preparation function ---
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
# ---- Plot density for celltype 4 ----
for(dt in dts) {
  for (paramset in c('P1','P3','P7')) {
    pval.allpt <- do.call(rbind, lapply(patterns, function(pt) {
      dataset <- sprintf("sim_%s-%s-%s-rep1", dt, pt, paramset)
      res <- get_pvalue_wide(dataset, svg_id)
      dat.pval.wide <- res$dat.pval.wide
      prep <- prepare_plot_density(dat.pval.wide,
                                   celltype = 4,
                                   sel1.range = 1:200,
                                   sel2.range = 201:1200,
                                   rank_scope = "global")
      pval_df <- prep$pval_df
      pval_df$pattern <- pt
      return(pval_df)
    }))
    pval.allpt$pattern <- stringr::str_to_sentence(pval.allpt$pattern)
    # Generate density plot
    ggplot(pval.allpt, aes(x = pval, fill = group, color = group)) +
      geom_density(alpha = 0.35, size = 0.6) +
      facet_wrap(pattern ~ method, scales = "free", ncol = 4) +
      theme_minimal(base_size = base_font_size) +
      scale_fill_manual(values = my_colors) +
      scale_color_manual(values = my_colors) +
      labs(x = "P value",
           y = "Density",
           fill = "gene id",
           color = "gene id") +
      my_theme +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right",
        legend.key.size = unit(0.05, "in"), 
        legend.key.width = unit(0.05, "in"),
        plot.margin = margin(4, 1, 4, 0)
      )
    
    ggsave(sprintf('./Fig/s/sim-density-%s-%s.pdf', dt, paramset), width = 6.69, height = 8.2, units = "in")  
  }
}

# ---- Plot density for celltype 6 ----
for(dt in dts) {
  for (paramset in c('P1','P3','P7')) {
    pval.allpt <- do.call(rbind, lapply(patterns, function(pt) {
      dataset <- sprintf("sim_%s-%s-%s-rep1", dt, pt, paramset)
      res <- get_pvalue_wide(dataset, svg_id)
      dat.pval.wide <- res$dat.pval.wide
      prep <- prepare_plot_density(dat.pval.wide,
                                   celltype = 6,
                                   sel1.range = 1:75,
                                   sel2.range = 201:1200,
                                   rank_scope = "global")
      pval_df <- prep$pval_df
      pval_df$pattern <- pt
      return(pval_df)
    }))
    pval.allpt$pattern <- stringr::str_to_sentence(pval.allpt$pattern)
    # Generate density plot
    ggplot(pval.allpt, aes(x = pval, fill = group, color = group)) +
      geom_density(alpha = 0.35, size = 0.6) +
      facet_wrap(pattern ~ method, scales = "free", ncol = 4) +
      scale_y_continuous(labels = function(x) sprintf("%.4f", x)) +  
      theme_minimal(base_size = base_font_size) +
      scale_fill_manual(values = my_colors) +
      scale_color_manual(values = my_colors) +
      labs(x = "P value",
           y = "Density",
           fill = "gene id",
           color = "gene id") +
      my_theme +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right",
        legend.key.size = unit(0.05, "in"), 
        legend.key.width = unit(0.05, "in"),
        plot.margin = margin(4, 1, 4, 0)
      )
    
    ggsave(sprintf('./Fig/s/sim-density-ct6-%s-%s.pdf', dt, paramset), width = 6.69, height = 8.2, units = "in")  
  }
}
