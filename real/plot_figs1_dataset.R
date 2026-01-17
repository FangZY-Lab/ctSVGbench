library(ggplot2)
library(ggrepel)
library(tidyr)
source('F:/ctSVGbench/my_theme.R')
datasets <- c(
  "MERFISH_hypothalamus",
  "SeqFish+_cortex",
  "SeqFish+_mouse_OB",
  "Slide-seq_tumor",
  "Slide-seqV2_hippocampus",
  "Slide-seqV2_melanoma_GSM6025935_MBM05_rep1",
  "Slide-seqV2_melanoma_GSM6025936_MBM05_rep2",
  "Slide-seqV2_melanoma_GSM6025937_MBM05_rep3",
  "Slide-seqV2_melanoma_GSM6025938_MBM06",
  "Slide-seqV2_melanoma_GSM6025939_MBM07",
  "Slide-seqV2_melanoma_GSM6025940_MBM08",
  "Slide-seqV2_melanoma_GSM6025944_MBM13",
  "Slide-seqV2_melanoma_GSM6025949_ECM08",
  "Slide-seqV2_melanoma_GSM6025950_ECM10",
  "Slide-seqV2_mouseOB",
  "ST_developmental_heart",
  "ST_PDAC",
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
  "StereoSeq_MDESTA",
  "Stereoseq_mosta_Dorsal_midbrain",
  "Stereoseq_mosta_E16.5_E1S3_whole_brain",
  "Visium_bladder",
  "Visium_intestine",
  "Visium_liver",
  "Visium_lymph_node",
  "Visium_mousebrain",
  "Visium_pancreas",
  "Visium_skin",
  "Visium_spleen",
  "Visium_tail",
  "VisiumHD_LUAD_2431",
  "VisiumHD_LUAD_6123",
  "VisiumHD_LUAD_6976",
  "VisiumHD_LUSC_5488",
  "VisiumHD_LUSC_7437",
  "VisiumHD_LUSC_7941"
)


data_info <- read.csv('data_info.csv')
ggplot(data_info, aes(x = n_gene, y = n_pos, color = tech_platform, label = tissue)) +
  geom_point(size = 1, alpha = 0.7) +
  geom_text_repel(size = 2, max.overlaps = 20) +  
  theme_bw() +
  my_theme+
  theme(legend.position = "bottom",
        plot.margin = margin(10, 10, 10, 10),
        legend.key.size   = unit(0.2, "in"), 
        legend.key.width = unit(0.2, "in"))+
  labs(x = "Number of genes", y = "Number of spots", color = "Platform")+
  scale_color_manual(values = brewer.pal(n = 10, name = "Set2"))+
  guides(color= guide_legend(
    nrow = 2))
ggsave("./Fig/s/datasets.pdf",width=6.69,height = 4)
