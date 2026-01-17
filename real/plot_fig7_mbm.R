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
dataset <- c("Slide-seqV2_melanoma_GSM6025944_MBM13") 

methods <- c('C-SIDE','spVC','CELINA','STANCE','CTSV')
file <- sprintf('myRCTD_%s.rds',dataset)
puck <- readRDS(here('real','puck',sprintf('myRCTD_%s.rds',dataset)))

res.cside=readRDS(here('real','res',sprintf('%s-C-SIDE.rds',dataset)))
res.celina=readRDS(here('real','res',sprintf('%s-CELINA.rds',dataset)))
res.stance=readRDS(here('real','res',sprintf('%s-STANCE.rds',dataset)))
res.ctsv=readRDS(here('real','res',sprintf('%s-CTSV.rds',dataset)))

spVC=readRDS(here('real','res',sprintf('%s-spVC.rds',dataset)))
prop <- readRDS(here('real','prop',sprintf('myRCTD_%s.rds',dataset))) 
pos <- readRDS(here('real','pos_subset',sprintf('myRCTD_%s.rds',dataset))) 

genes.v=names(spVC$results.varying)
idx=match(names(res.celina),colnames(prop))

res.spvc <- lapply(idx,function(ct){
  pval=sapply(spVC$results.varying[genes.v],function(x){
    x$p.value[paste0("gamma_X", ct)]
  })
  names(pval)=sapply(strsplit(names(pval),"\\."),"[[",1)
  data.frame(pval = na.omit(pval))
})

names(res.spvc) <- names(res.celina)

all_lists <- list(
  CSIDE = res.cside, 
  spVC = res.spvc, 
  Celina = res.celina, 
  STANCE = res.stance,
  CTSV = res.ctsv
)

all_genes <- unique(unlist(lapply(all_lists, function(lst) {
  unlist(lapply(lst, rownames))
})))

dat.pval <- do.call(rbind,lapply(names(all_lists),function(method){
  lst <- all_lists[[method]]
  n <- min(3, length(names(lst)))
  list=lapply(names(lst)[1:n],function(i){  
    df <- lst[[i]]
    if(nrow(df)>0){ 
      pval_cols <- grep("p[_\\.]?val|pvalue|pval", colnames(df), ignore.case = TRUE) 
      
      if (length(pval_cols)!=1) {
        stop(sprintf("Data frame %s[[%d]] does not contain column pval ", method, i))
      }
      pvals <- p.adjust(df[,pval_cols],method = "BH")
      names(pvals) <- rownames(df)
      
      pvals_full <- rep(1, length(all_genes))
      names(pvals_full) <- all_genes
      
      pvals_full[names(pvals)] <- pvals
      pvals_full <- data.frame(pval=pvals_full)
      pvals_full$gene=paste0(make.names(i),'_',rownames(pvals_full))
      return(pvals_full)
      
    }else{
      pvals_full <- rep(1, length(all_genes))
      names(pvals_full) <- all_genes
      pvals_full <- data.frame(pval=pvals_full)
      pvals_full$gene=paste0(make.names(i),'_',rownames(pvals_full))
      return(pvals_full)
    }
  })
  list=list[!sapply(list,is.null)]
  if(length(list)>0){
    df.pval <- do.call(rbind,list)
    df.pval$method <- method
    return(df.pval)
  }else{
    return(NULL)
  }
}))

str(dat.pval)
dat.pval$method=ifelse(dat.pval$method=="CSIDE","C-SIDE",dat.pval$method)
alpha <- 0.05

dat.pval <- dat.pval %>%
  mutate(gene = sub("\\.", " ", gene)) %>%  
  separate(gene, into = c("cell_type", "gene_name"), sep = "_")  
df.sig <-  dat.pval %>%
  filter(pval < alpha)

# --- venn ---
sig_gene_list <- dat.pval %>%
  filter(pval < alpha) %>%
  group_by(method) %>%
  summarise(genes = list(unique(gene_name))) %>%
  deframe()  

# --- significant gene counts ---
all_combinations <- expand.grid(
  method = unique(dat.pval$method),
  cell_type = unique(dat.pval$cell_type),
  stringsAsFactors = FALSE
)

gene_counts <- dat.pval %>% 
  filter(pval < 0.05) %>%
  group_by(method, cell_type) %>%
  summarise(sig_gene_count = n_distinct(gene_name), .groups = "drop")

gene_counts_complete <- all_combinations %>%
  left_join(gene_counts, by = c("method", "cell_type")) %>%
  replace_na(list(sig_gene_count = 0))

intersection_info <- df.sig %>%
  group_by(cell_type, gene_name) %>%
  summarise(methods_in = list(unique(method)), method_count = n(), .groups = "drop")

shared2 <- intersection_info %>% filter(method_count == 2)
shared3 <- intersection_info %>% filter(method_count == 3)
shared4 <- intersection_info %>% filter(method_count == 4)



# --- Enrichment analysis ---
gene_sets_sig <- df.sig %>%
  group_by(method, cell_type) %>%
  filter(cell_type=='Tumor cells') %>% 
  summarise(genes = list(gene_name), .groups = "drop")

ego_results <- list()

for (i in seq_len(nrow(gene_sets_sig))) {
  method_i <- gene_sets_sig$method[i]
  cell_i <- gene_sets_sig$cell_type[i]
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
  ego_results[[paste(method_i, cell_i, sep = " ")]] <- ego_simple
}

all_terms <- unique(c(ego_results[[1]]$ID, ego_results[[2]]$ID, ego_results[[3]]$ID))

# --- Top genes ---
top_gene_list <- df.sig %>%
  filter(cell_type=='Tumor cells') %>% 
  arrange(method, pval) %>%
  group_by(method) %>%
  slice_head(n = 7) %>%
  summarise(top_genes = list(unique(gene_name))) %>%
  deframe()  

common_top_genes <- Reduce(intersect, top_gene_list[c("C-SIDE","STANCE","Celina")])
print(common_top_genes)

genes_to_plot <- c("PMEL", "VIM", "GAPDH", "ENO1")
expr_normalized <- scater::normalizeCounts(puck@counts)
expr_mat <- as.matrix(expr_normalized[genes_to_plot,rownames(pos) , drop = FALSE])

expr_df <- as.data.frame(t(expr_mat)) %>%
  rownames_to_column("spot") %>%
  pivot_longer(-spot, names_to = "gene", values_to = "expression")

coords <- puck@coords %>%
  as.data.frame() %>%
  rownames_to_column("spot")

plot_expr <- left_join(expr_df, coords, by = "spot")
plot_expr_log <- plot_expr %>%
  mutate(log_expr = log10(1 + expression))

# --- Cell proportion ---
cell_type_order <- rev(colnames(prop))
prop <- as.data.frame(prop)
colnames(prop) <- sub("\\.", " ", colnames(prop))

prop_long <- prop %>%
  rownames_to_column("coord") %>%
  pivot_longer(-coord, names_to = "cell_type", values_to = "value") %>%
  filter(value > 0.3) %>% 
  left_join(pos %>% rownames_to_column("coord"), by = "coord")
  
prop_long$cell_type <- factor(prop_long$cell_type, levels = cell_type_order)

# --- Save environment ---

my_theme <- 
  theme(
    axis.text   = element_text(size = 7),
    axis.title  = element_text(size = 7),
    legend.text = element_text(size = 7), 
    legend.title = element_blank(),
    strip.text  = element_text(size = 7)
  )

methods <- c("C-SIDE", "Celina", "spVC", "STANCE","CTSV","ctSVG")
colors <- c("#0073C2FF","#EFC000FF","#868686FF","#CD534CFF","#7AA6DCFF","#003C67FF")
method_colors <- setNames(colors, methods)


##
prop.use <- prop[, -which(colnames(prop) == "Low quality.cells")]
avg_props <- colMeans(prop.use, na.rm = TRUE)

avg_df <- data.frame(
  CellType = names(avg_props), 
  Proportion = avg_props      
)
p0 <- ggplot(avg_df, aes(x = '', y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity") +
  labs(y="",x = 'MBM13',fill="Cell type") +
  theme_void(base_size = 7) +
  theme(legend.position = 'right',
        axis.title = element_text(),
        plot.margin = margin(4, 2, 4, 0),
        legend.spacing.x = unit(0.05, "in"),
        legend.key.size   = unit(0.05, "in"), 
        legend.key.width = unit(0.05, "in") )+
  guides(fill= guide_legend(ncol = 1))  +
  ggtitle('Proportion')+
  scale_fill_manual(values = rev(RColorBrewer::brewer.pal(8,"Set2")),
                    guide = guide_colorbar(
                      title.position = "top",
                      label.hjust = 0.5
                    ))
p0 

common_cells <- intersect(rownames(prop.use), rownames(pos))
prop_filtered <- prop.use[common_cells, ]
pos_filtered <- pos[common_cells, ]

max_cell_types <- apply(prop_filtered, 1, function(x) {
  colnames(prop_filtered)[which.max(x)]
})

plot_data <- data.frame(
  cell_id = common_cells,
  x = pos_filtered$x,
  y = pos_filtered$y,
  cell_type = max_cell_types,
  max_proportion = apply(prop_filtered, 1, max)
)

p1 <- ggplot(plot_data, aes(x = x, y = y, color = cell_type)) +
  geom_point(size = 0.1, alpha = 0.8) +
  scale_color_manual(values = rev(RColorBrewer::brewer.pal(8,"Set2"))) +
  theme_void(base_size = 7)+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none",
        plot.margin = margin(0, 0, 0, 0)
  ) +
  # guides(color = guide_legend(override.aes = list(size = 3))) +
  coord_fixed()
# labs(color="cell type")+
p1
ggsave('cell.pdf', p1, width = 3.5, height = 3.5, units = "in")

p2 <- gene_counts_complete %>% 
  subset(cell_type=='Tumor cells') %>% 
  ggplot( aes(x = cell_type, y = sig_gene_count, fill = method)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "", y = "Gene count") +
  geom_text(aes(label = sig_gene_count),
            position = position_dodge(width = 0.9),  
            vjust = -0.1, 
            size = 2) +
  theme_minimal() +            
  my_theme +
  scale_fill_manual(values = method_colors)+
  theme(    plot.margin = margin(1, 0, 0, 0),
            legend.position = c(0.1,0.85))
p2
library(ggvenn)
# p3：Venn
p3 <- ggvenn(
  sig_gene_list,
  fill_color =  c("#0073C2FF","#7AA6DCFF","#EFC000FF" ,"#CD534CFF"),
  stroke_size = 0.1,
  show_percentage = FALSE,
  fill_alpha = 0.8,
  set_name_size = 2,   
  text_size = 2
  
)+theme( plot.margin = margin(0, 0, 0, 0))
p3
ggsave('Fig/p3.pdf')
# p4：spatial expr
p4 <- ggplot(plot_expr_log, aes(x = x, y = y, color = log_expr),alpha=0.5) +
  geom_point(size = 0.01) +
  scale_color_gradientn( colors = rev(brewer.pal(11, "RdYlBu"))) +  
  # scale_color_gradientn(
  #   colors = c("blue", "white", "red"),
  #   name = "log1p(Expression)"
  # ) + 
  coord_fixed() +
  facet_wrap(~gene) +
  theme_void()+
  my_theme+
  theme(
    legend.spacing.x = unit(0.05, "in"),
    # legend.key.size   = unit(0.05, "in"), 
    legend.key.width = unit(0.05, "in") ,
    plot.margin = margin(0, 2, 0, 0),
    strip.text  = element_text(size = 7, face = "plain",),
    axis.title = element_blank(),       
    axis.text = element_blank(),        
    axis.ticks = element_blank(),       
    axis.line = element_blank()        
  )
p4

library(dplyr)
results_unsim <- list()
library(clusterProfiler)
for (i in seq_len(nrow(gene_sets_sig))) {
  method_i <- gene_sets_sig$method[i]
  cell_i <- gene_sets_sig$cell_type[i]
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
  
  results_unsim[[paste(method_i, cell_i, sep = " ")]] <- ego
}
res <- do.call(rbind, lapply(names(results_unsim), function(method_name) {
  result_df <- results_unsim[[method_name]][, c('ID','Description', 'Count', 'p.adjust')]
  result_df$Method <- strsplit(method_name," ")[[1]][1]
  return(result_df)
}))

res_filtered <- res %>% filter(ID %in% all_terms)

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
pgo <- plot_grid(pgo_p+theme(legend.position = "none"),pgo_l,nrow=2,rel_heights = c(1,0.15))

library(cowplot)
library(ggplot2)

# -----------------------------
# 1. top
# -----------------------------
top_row <- plot_grid(
  p1, p0, p3, p4,
  ncol = 4,
  rel_widths = c(1.2,0.6,1.2,1.2), 
  labels = c("A","","B","C")
)

# -----------------------------
# 2. final
# -----------------------------
final_plot <- plot_grid(
  top_row, pgo,
  ncol = 1,
  rel_heights = c(1, 2), 
  labels = c('','D')
  
)

ggsave("./Fig/Fig7.pdf", final_plot, width = 6.67, height = 6.6, units = "in")


