library(dplyr)
library(tidyr)
library(stringr)
library(reshape2)
library(ggplot2)
patterns <- c("pathology","hotspot", "stripe", "gradient",  "periodic", "neighbor")
paramset = 'P1'

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
)#17个

pval_long <- do.call(rbind, lapply(dts_sp, function(dt) {
  do.call(rbind, lapply(patterns, function(pt) {
  ds <- sprintf("sim_%s-%s-%s-rep1", dt, pt, paramset)
  file.orign <- sprintf('myRCTD_%s.rds', dt)
  prop.orign <- readRDS(here('sim','prop', file.orign))
  
  # Make column names syntactically valid
  colnames(prop.orign) <- make.names(colnames(prop.orign))
  
  # Select top cell types by abundance
  prop.use <- prop.orign[, names(tail(sort(colSums(prop.orign))))]
  
  # Ensure at least 6 cell types; add zero columns if needed
  if (ncol(prop.use) < 6) {
    cols.to.add <- 6 - ncol(prop.use)
    new.cols <- as.data.frame(matrix(0, nrow = nrow(prop.use), ncol = cols.to.add))
    colnames(new.cols) <- paste0("new", seq_len(cols.to.add))
    prop.use <- cbind(new.cols, prop.use)
  }
  
  # Read spatial coordinates and boundary
  pos.use <- readRDS(here('sim','pos', file.orign))
  colnames(pos.use) <- c("x","y")
  spot <- intersect(rownames(pos.use), rownames(prop.use))
  pos.use <- pos.use[spot, ]
  prop.use <- prop.use[spot, ]
  prop <- colSums(prop.use)[4]/sum(colSums(prop.use))
  res <- get_pvalue_wide(ds, svg_id)   # wide: gene x method
  
  # wide -> long
  tmp <- reshape2::melt(as.matrix(res$dat.pval.wide))
  colnames(tmp) <- c("gene", "method", "pval")
  tmp$pval <- as.numeric(tmp$pval)
  tmp$dataset <- ds
  tmp$pattern <- pt
  tmp$prop <- prop
  tmp

  }))
}))

pval_long <- pval_long %>%
  mutate(
    celltype = as.integer(str_extract(gene, "(?<=celltype)\\d+")),
    gene_id  = as.integer(str_extract(gene, "(?<=gene)\\d+"))
  )

pval_long <- pval_long %>%
  filter(celltype==4) %>% 
  mutate(
    neg_class = case_when(
      celltype == 4 & gene_id >= 1   & gene_id <= 200 ~ "non-target ctSVG",
      # celltype == 5 & gene_id >= 76  & gene_id <= 150 ~ "affected",
      # celltype == 6 & gene_id >= 1   & gene_id <= 75  ~ "affected",
      celltype %in% c(4, 5, 6) & gene_id >= 200 & gene_id <= 1200 ~ "non-ctSVG",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(neg_class))

fpr_df <- pval_long %>%
  group_by(dataset, neg_class,pattern, method, prop) %>%
  summarise(
    total_n = n(),
    sig_n = sum(pval < 0.05, na.rm = TRUE),
    fpr = sig_n / total_n,
    .groups = "drop"
  )

fpr_df %>% 
ggplot(aes(x = prop, y = fpr, color = method)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_smooth(method='lm')+
  facet_wrap(~neg_class)+
  # facet_wrap(~ pattern,ncol = 6) +
  labs(
    x = "Proportion (prop)",
    y = "False Positive Rate (FPR)",
    color = "Method"
  ) +
  theme_bw() +
  my_theme+
scale_color_manual(values = method_colors)
ggsave("Fig/s/ct4_prop.pdf", width = 6.69, height = 3)
