library(ggplot2)
library(ggrepel)
library(tidyr)
source('F:/ctSVGbench/my_theme.R')

# res.df$Datasets <- res.df$dataset
# res.df <- res.df %>%
#   separate(dataset, into = c("platform", "tissue"), sep = "_", extra = "merge") %>% 
#   select(platform, tissue, n_gene, n_spot,Datasets) %>% 
#   distinct()
res.df <- read.csv('./dataset.csv')
ggplot(res.df, aes(x = n_gene, y = n_spot, color = platform, label = tissue)) +
  geom_point(size = 1, alpha = 0.7) +
  geom_text_repel(size = 2, max.overlaps = 20) +  
  theme_bw() +
  my_theme+
  theme(legend.position = "bottom",
        plot.margin = margin(10, 10, 10, 10),
        legend.key.size   = unit(0.2, "in"), 
        legend.key.width = unit(0.2, "in"))+
  labs(x = "Number of genes", y = "Number of spots", color = "Platform")+
  scale_color_manual(values = brewer.pal(n = 8, name = "Set2"))+
  guides(color= guide_legend(
    nrow = 2))
ggsave("./Fig/s/datasets.pdf",width=6.69,height = 4)
