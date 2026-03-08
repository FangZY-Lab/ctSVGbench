library(dplyr)
library(purrr)
library(ggplot2)
library(scales)
library(here)
library(tidyr)
library(dplyr)
library(tibble)
library(tidyverse)
library(data.table)
library(cowplot)
library(ggsci)
library(ggpubr)
library(ggtext)
library(ggprism)
library(ggpmisc)
library(clusterProfiler)
library(org.Hs.eg.db)
library(VennDiagram)
library(UpSetR)
library(ggvenn)
library(DOSE)
library(Matrix)
library(RColorBrewer)
library(patchwork)
library(dplyr)

library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)

genes_list_sym <- ctsvg %>%
  filter(Group %in% c("LUAD_STANCE", "LUSC_STANCE")) %>%
  dplyr::select(Group, SYMBOL) %>%
  distinct() %>%
  filter(!is.na(SYMBOL)) %>%
  group_by(Group) %>%
  summarise(genes = list(SYMBOL), .groups = "drop")

# SYMBOL -> ENTREZID
genes_list_entrez <- lapply(genes_list_sym$genes, function(x) {
  bitr(x,
       fromType = "SYMBOL",
       toType   = "ENTREZID",
       OrgDb    = org.Hs.eg.db)$ENTREZID |> unique()
})

names(genes_list_entrez) <- genes_list_sym$Group
ckegg <- compareCluster(
  geneCluster   = genes_list_entrez,
  fun           = "enrichKEGG",
  organism      = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2
)

ckegg <- setReadable(ckegg,
                     OrgDb   = org.Hs.eg.db,
                     keyType = "ENTREZID")


dotplot(ckegg, showCategory = 15) +
  scale_size(range = c(1, 5)) +  
  theme_minimal(base_size = 7) +
  my_theme +
  theme(
    axis.text.y = element_text(size = 7),
    axis.text.x = element_text(size = 7))

ggsave('Fig/s/lung_stance.pdf',width = 6.69, height = 4.5)



#consvered gene
gene_intersection <- gene_sets2 %>% 
  filter(cancer=='LUAD') %>%
  group_by(cancer, method) %>%
  summarise(
    n_samples = n_distinct(cancer_sample),
    intersect_genes = list(purrr::reduce(genes, intersect)),
    n_intersect = length(intersect_genes[[1]]),
    .groups = "drop"
  )
gene_sets_sig <- gene_intersection %>%
  dplyr::select(genes = intersect_genes, method)


ego_results <- list()

for (i in seq_len(nrow(gene_sets_sig))) {
  method_i <- gene_sets_sig$method[i]
  gene_symbols <- gene_sets_sig$genes[[i]]
  
  gene_ids <- bitr(gene_symbols,
                   fromType = "SYMBOL",
                   toType = "ENTREZID",
                   OrgDb = org.Hs.eg.db)
  
  if (nrow(gene_ids) == 0) next  
  
  ego <- enrichGO(gene = gene_ids$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  ont = "BP",
                  pAdjustMethod = "BH",
                  qvalueCutoff = 0.05,
                  readable = TRUE)
  
  ego_simple <- simplify(ego,
                         cutoff = 0.3,
                         by = "p.adjust",
                         select_fun = min,
                         measure = "Wang")
  
  ego_simple <- filter(ego_simple,p.adjust<0.05)
  ego_results[[paste(method_i, sep = " ")]] <- ego_simple
}

all_terms <- unique(unlist(lapply(ego_results, function(x) x$ID)))

library(dplyr)
results_unsim <- list()
library(clusterProfiler)
for (i in seq_len(nrow(gene_sets_sig))) {
  method_i <- gene_sets_sig$method[i]
  gene_symbols <- gene_sets_sig$genes[[i]]
  
  # SYMBOL → ENTREZID
  gene_ids <- bitr(gene_symbols,
                   fromType = "SYMBOL",
                   toType = "ENTREZID",
                   OrgDb = org.Hs.eg.db)
  
  if (nrow(gene_ids) == 0) next  
  
  # GO enrichment
  ego <- enrichGO(gene = gene_ids$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  ont = "BP",
                  pAdjustMethod = "BH",
                  qvalueCutoff = 0.05,
                  readable = TRUE)
  
  results_unsim[[paste(method_i, sep = " ")]] <- ego
}
res <- do.call(rbind, lapply(names(results_unsim), function(method_name) {
  result_df <- results_unsim[[method_name]][, c('ID','Description', 'Count', 'p.adjust')]
  result_df$Method <- strsplit(method_name," ")[[1]][1]
  return(result_df)
}))

res_filtered <- res %>% 
  filter(ID %in% all_terms) 

pgo_p <- ggplot(res_filtered, aes(x = Method, y = reorder(Description, Count))) +
  geom_point(aes(size = Count, color = -log10(p.adjust))) +
  labs(x = "Method", y = "", title = "") +
  theme_bw() +
  # my_theme+
  scale_color_gradientn(
    colors = rev(RColorBrewer::brewer.pal(11, "Spectral")),
    guide = guide_colorbar(
      title.position = "top",
      barwidth = unit(1.5, "in"),
      barheight = unit(0.05, "in"),
      label.hjust = 0.5
    )
  )+
  scale_size_continuous(
    name = "Gene Count", 
    range = c(0.8, 4),
    guide = guide_legend(
      title.position = "top"
    )
  ) +theme(legend.position = "bottom",
           legend.spacing.x = unit(0.1, "in"),
           plot.margin = margin(0, 15, 0, 0),
           legend.title = element_text(),
           legend.box.margin = margin(1, 1, 1, 1),
           legend.margin     = margin(1, 1, 1, 1),)+
  guides(fill= guide_legend(nrow = 1),
         size  = guide_legend(nrow = 1))
pgo_p

library(cowplot)
legend_size <- ggpubr::get_legend(
  pgo_p + guides(color = "none") +  
    theme(legend.position = "bottom",
          legend.spacing.x = unit(0.01, "in"),
    )+my_theme+
    scale_size_continuous(
      name = "Gene Count", 
      range = c(0.8, 4),
      guide = guide_legend(
        title.position = "top",
        label.position = "bottom",  
        direction = "horizontal" ,   
        label.hjust = 0.5
        
      )
    )
)

legend_color <- ggpubr::get_legend(
  pgo_p + guides(size = "none") +  
    my_theme+
    theme(legend.position = "bottom",
          plot.margin = margin(0, 0, 0, 0),
          legend.box.margin = margin(0, 0, 0, 0),
          legend.margin     = margin(0, 0, 0, 0))
)

pgo_l <- plot_grid(legend_size, legend_color, nrow = 1, rel_widths = c(1, 2))
pgo <- plot_grid(pgo_p+theme(
  legend.position = "none",
  axis.text.x = element_text(angle = 45, hjust = 1)
  ),
  pgo_l,nrow=2,rel_heights = c(1,0.15))

ggsave('Fig/s/lung_conserv.pdf',width = 6.69, height = 7)


#sig gene
gene_sets_sig <- gene_sets2 %>%
  filter(cancer == 'LUAD') %>%  
  mutate(method = paste(cancer, method, sep = "_")) %>%
  group_by(method) %>%  
  summarise(genes_u = list(Reduce(union, genes))) %>% 
  dplyr::select(genes = genes_u, method)

ego_results <- list()

for (i in seq_len(nrow(gene_sets_sig))) {
  method_i <- gene_sets_sig$method[i]
  gene_symbols <- gene_sets_sig$genes[[i]]
  
  gene_ids <- bitr(gene_symbols,
                   fromType = "SYMBOL",
                   toType = "ENTREZID",
                   OrgDb = org.Hs.eg.db)
  
  if (nrow(gene_ids) == 0) next  
  
  ego <- enrichGO(gene = gene_ids$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  ont = "BP",
                  pAdjustMethod = "BH",
                  qvalueCutoff = 0.05,
                  readable = TRUE)
  
  ego_simple <- simplify(ego,
                         cutoff = 0.3,
                         by = "p.adjust",
                         select_fun = min,
                         measure = "Wang")
  
  ego_simple <- filter(ego_simple,p.adjust<0.05)
  ego_results[[paste(method_i, sep = " ")]] <- ego_simple
}

all_terms <- unique(unlist(lapply(ego_results, function(x) x$ID)))

library(dplyr)
results_unsim <- list()
library(clusterProfiler)
for (i in seq_len(nrow(gene_sets_sig))) {
  method_i <- gene_sets_sig$method[i]
  gene_symbols <- gene_sets_sig$genes[[i]]
  
  # SYMBOL → ENTREZID
  gene_ids <- bitr(gene_symbols,
                   fromType = "SYMBOL",
                   toType = "ENTREZID",
                   OrgDb = org.Hs.eg.db)
  
  if (nrow(gene_ids) == 0) next  
  
  # GO enrichment
  ego <- enrichGO(gene = gene_ids$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  ont = "BP",
                  pAdjustMethod = "BH",
                  qvalueCutoff = 0.05,
                  readable = TRUE)
  
  results_unsim[[paste(method_i, sep = " ")]] <- ego
}
res <- do.call(rbind, lapply(names(results_unsim), function(method_name) {
  result_df <- results_unsim[[method_name]][, c('ID','Description', 'Count', 'p.adjust')]
  result_df$Method <- strsplit(method_name," ")[[1]][1]
  return(result_df)
}))

res_filtered <- res %>% 
  filter(ID %in% all_terms) 

pgo_p <- ggplot(res_filtered, aes(x = Method, y = reorder(Description, Count))) +
  geom_point(aes(size = Count, color = -log10(p.adjust))) +
  labs(x = "Method", y = "", title = "") +
  theme_bw() +
  # my_theme+
  scale_color_gradientn(
    colors = rev(RColorBrewer::brewer.pal(11, "Spectral")),
    guide = guide_colorbar(
      title.position = "top",
      barwidth = unit(1.5, "in"),
      barheight = unit(0.05, "in"),
      label.hjust = 0.5
    )
  )+
  scale_size_continuous(
    name = "Gene Count", 
    range = c(0.8, 4),
    guide = guide_legend(
      title.position = "top"
    )
  ) +theme(legend.position = "bottom",
           legend.spacing.x = unit(0.1, "in"),
           plot.margin = margin(0, 15, 0, 0),
           legend.title = element_text(),
           legend.box.margin = margin(1, 1, 1, 1),
           legend.margin     = margin(1, 1, 1, 1),)+
  guides(fill= guide_legend(nrow = 1),
         size  = guide_legend(nrow = 1))
pgo_p

library(cowplot)
legend_size <- ggpubr::get_legend(
  pgo_p + guides(color = "none") +  
    theme(legend.position = "bottom",
          legend.spacing.x = unit(0.01, "in"),
    )+my_theme+
    scale_size_continuous(
      name = "Gene Count", 
      range = c(0.8, 4),
      guide = guide_legend(
        title.position = "top",
        label.position = "bottom",  
        direction = "horizontal" ,   
        label.hjust = 0.5
        
      )
    )
)

legend_color <- ggpubr::get_legend(
  pgo_p + guides(size = "none") +  
    my_theme+
    theme(legend.position = "bottom",
          plot.margin = margin(0, 0, 0, 0),
          legend.box.margin = margin(0, 0, 0, 0),
          legend.margin     = margin(0, 0, 0, 0))
)

pgo_l <- plot_grid(legend_size, legend_color, nrow = 1, rel_widths = c(1, 2))
pgo <- plot_grid(pgo_p+theme(
  legend.position = "none",
  axis.text.x = element_text(angle = 45, hjust = 1)
  ),
  pgo_l,nrow=2,rel_heights = c(1,0.15))

ggsave('Fig/s/lung_sig.pdf',width = 6.69, height = 7)


