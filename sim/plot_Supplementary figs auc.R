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

# Define datasets and patterns
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
  "Stereoseq_mosta_E16.5_E1S3_whole_brain",
  "Stereoseq_mosta_Dorsal_midbrain"    
)#20

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
)#17

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

for(paramset in c('P2','P3') ){
  # ---- Calculate AUC for all datasets and patterns ----
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
      }))
      }))
  
  
  auc$methods <- factor(auc$methods,levels = c("spVC_1","spVC_2","C-SIDE", "Celina", "STANCE","CTSV","ctSVG"))
  auc$pattern <- stringr::str_to_sentence(auc$pattern)
  
  source('./sim/pre_fig3B.R')
  
  p1
  ggsave(sprintf('./Fig/s/sim-auc-%s.pdf',paramset),p1, width = 6.69, height = 4, units = "in")
}
