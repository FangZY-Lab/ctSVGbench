
calc_false_positive_rate <- function(dat.pval.wide, alpha = 0.05) {
  library(dplyr)
  library(tidyr)
  library(stringr)
  
  # 1. long format
  pval_long <- dat.pval.wide %>%
    as.matrix() %>%
    reshape2::melt()
  
  colnames(pval_long) <- c("gene", "method", "pval")
  pval_long$pval <- as.numeric(pval_long$pval)
  
  # 2. parse celltype and gene index
  pval_long <- pval_long %>%
    mutate(
      celltype = as.integer(str_extract(gene, "(?<=celltype)\\d+")),
      gene_id  = as.integer(str_extract(gene, "(?<=gene)\\d+"))
    )
  
  # 3. define negative gene classes
  pval_long <- pval_long %>%
    mutate(
      neg_class = case_when(
        # potentially affected negatives
        celltype == 4 & gene_id >= 1   & gene_id <= 200 ~ "affected",
        celltype == 5 & gene_id >= 76  & gene_id <= 150 ~ "affected",
        celltype == 6 & gene_id >= 1   & gene_id <= 75  ~ "affected",
        
        # pure negatives
        celltype %in% c(4, 5, 6) & gene_id >= 200 & gene_id <= 1200 ~ "non-affected",
        
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(neg_class))
  
  # 4. calculate false positive rate
  res_by_ct <- pval_long %>%
    group_by(method, celltype, neg_class) %>%
    summarise(
      n_total = n(),
      n_fp    = sum(pval < alpha, na.rm = TRUE),
      fpr     = n_fp / n_total,
      .groups = "drop"
    ) %>%
    mutate(level = "celltype")
  
  res_overall <- pval_long %>%
    group_by(method, neg_class) %>%
    summarise(
      n_total = n(),
      n_fp    = sum(pval < alpha, na.rm = TRUE),
      fpr     = n_fp / n_total,
      .groups = "drop"
    ) %>%
    mutate(
      celltype = "all",
      level = "overall"
    )
  
  res_by_ct <- res_by_ct %>%
    mutate(celltype = as.character(celltype))
  res <- bind_rows(res_by_ct, res_overall)
  
  return(res)
}
