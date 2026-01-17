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

ggsave('Fig/s/lung_stance.pdf',width = 6.69, height = 3)
