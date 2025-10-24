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
source(here('real','utils',"rotate_bench.R"))
source('F:/ctSVGbench/my_theme.R')

datasets <- c(
  "Slide-seq_tumor",
  "Slide-seqV2_hippocampus",
  "Slide-seqV2_mouseOB",
  "ST_developmental_heart",
  "Visium_mousebrain"
)

for (dataset in datasets){
  methods <- c('C-SIDE','spVC','CELINA','STANCE')
  
  conc.res <- do.call(rbind,lapply(methods,function(method){
    dat.pval.wide <- get_wide_pval(dataset,method)
    data <- get_conc(dat.pval.wide,dataset=dataset,method=method)
    return(data)
  }))
  
  conc.res <- conc.res %>%
    mutate(method = recode(method,
                           "C-SIDE" = "C-SIDE",
                           "CELINA" = "Celina",
                           "STANCE" = "STANCE",
                           "spVC" = "spVC"))
  
  bin_size <- 5  # Change this to 50 if you want bins of 50 ranks instead of 5
  conc.res$rank_bin <- floor((conc.res$rank - 1) / bin_size) + 1
  
  # Automatically generate bin labels
  conc.res$rank_bin_label <- with(conc.res, {
    start <- (as.numeric(rank_bin) - 1) * bin_size + 1
    end <- as.numeric(rank_bin) * bin_size
    paste0(start, "-", end)
  })
  
  # Convert to factor for ordered bin display
  conc.res$rank_bin_label <- factor(conc.res$rank_bin_label, 
                                    levels = unique(conc.res$rank_bin_label))
  
  conc.res1 <- subset(conc.res,rank<30)                                  
  p1 <- ggplot(conc.res1, aes(x = rank_bin_label, y = conc, fill = method)) +
    stat_summary(fun = mean, geom = "bar", 
                 position = position_dodge(0.8), 
                 width = 0.7, color = "black") +
    stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), 
                 geom = "errorbar", 
                 position = position_dodge(0.8), 
                 width = 0.3) +
    facet_grid(~angle2) +
    theme_minimal() +
    my_theme +
    scale_fill_manual(values = method_colors) +
    labs(
      title = "",
      x = sprintf("Rank Bin (%d ranks per bin)", bin_size),
      y = "Top of Concordance",
      fill = ""
    ) 
  
  conc.res2 <- subset(conc.res,rank<30)                                  
  
  p2 <- ggplot(conc.res2, aes(x = rank, y = conc, color = method)) +
    geom_line() +
    facet_grid(~angle2) +
    theme_minimal() +
    my_theme +
    scale_color_manual(values = method_colors) +
    labs(
      title = "",
      x = "Rank",
      y = "Top of Concordance",
      color = ""
    ) 
  p2
  
  plot_grid(plotlist = list(p2,p1), nrow = 2, rel_heights = c(1,1), labels = c("A", "B"))
  ggsave(sprintf("./Fig/s/conc-%s.pdf", dataset), width = 6.9, height = 4)
}

library(purrr)
library(dplyr)
library(data.table)

conc.all <- lapply(datasets, function(dataset) {
  methods <- c('C-SIDE','spVC','CELINA','STANCE')
  
  conc.res <- do.call(rbind, lapply(methods, function(method) {
    dat.pval.wide <- get_wide_pval(dataset, method)
    data <- get_conc(dat.pval.wide, dataset = dataset, method = method)
    return(data)
  }))
  
  conc.res <- conc.res %>%
    mutate(method = recode(method,
                           "C-SIDE" = "C-SIDE",
                           "CELINA" = "Celina",
                           "STANCE" = "STANCE",
                           "spVC" = "spVC"))
  
  conc.res <- as.data.table(conc.res)
  
  # Compute mean concordance for ranks 1–30 grouped by method and angle
  mean_conc <- conc.res[rank >= 1 & rank <= 30, 
                        .(mean_conc = mean(conc, na.rm = TRUE)),
                        by = .(method, angle2)]
  mean_conc$datasets <- dataset
  return(mean_conc)
})

conc.df <- rbindlist(conc.all, idcol = NULL)

library(ggplot2)
library(viridis)
library(data.table)
conc.df$method <- factor(conc.df$method,levels = c("C-SIDE","spVC", "Celina", "STANCE"))
# Heatmap of mean concordance across datasets
p1 <- ggplot(conc.df, aes(x = method, y = datasets, fill = mean_conc)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = round(mean_conc, 3)), size = 2) +
  scale_fill_gradientn(colors = c("#F7FBFF", "#DEEBF7", "#C6DBEF", "#9ECAE1", "#6BAED6", "#4292C6"),
                       guide = guide_colorbar(
                         barwidth = unit(0.1, "in"),
                         barheight = unit(1, "in"),
                         label.hjust = 0.5)) +
  labs(fill="Concordance") +
  facet_wrap(~ angle2, 
             scales = "fixed",
             labeller = as_labeller(c(
               `rotate.30` = "0° vs. 30°",  
               `rotate.90` = "0° vs. 90°"   
             ))) +
  theme_minimal() +
  my_theme +
  labs(y="Datasets",x="Methods")+
  theme(
    strip.text = element_text(size = 7),
    panel.grid = element_blank()
  )
p1

summary_df_rob2 <- conc.df %>%
  group_by(method) %>%
  summarise(Concordance = mean(mean_conc, na.rm = TRUE)) %>% 
  column_to_rownames('method')

# Similarity calculation function
calc_similarity <- function(dat.pval.wide,
                            jaccard = FALSE,
                            jaccard_mode = "pval",  # "pval", "top", or "all"
                            p.threshold = 0.05,
                            top_n = 50,
                            cor = FALSE,
                            cor_method = "pearson") {
  
  jaccard_mode <- match.arg(jaccard_mode, c("pval", "top", "all"))
  cols <- colnames(dat.pval.wide)
  n <- length(cols)
  result <- list()
  
  # Calculate Jaccard index
  if (jaccard) {
    binary_mat <- switch(jaccard_mode,
                         pval = (dat.pval.wide < p.threshold) * 1,
                         top  = apply(dat.pval.wide, 2, function(x) as.integer(rank(x, ties.method = "average") <= top_n)),
                         all  = matrix(1, nrow = nrow(dat.pval.wide), ncol = ncol(dat.pval.wide),
                                       dimnames = list(rownames(dat.pval.wide), cols))
    )
    jaccard_mat <- matrix(NA, n, n, dimnames = list(cols, cols))
    for (i in 1:1) {
      for (j in 2:3) {
        if (i == j) {
          jaccard_mat[i, j] <- 1
        } else {
          intersect_count <- sum(binary_mat[, i] & binary_mat[, j])
          union_count <- sum(binary_mat[, i] | binary_mat[, j])
          jaccard_mat[i, j] <- ifelse(union_count == 0, NA, intersect_count / union_count)
        }
      }
    }
    result$jaccard <- jaccard_mat
  }
  
  # Calculate correlation matrix
  if (cor) {
    cor_mat <- matrix(NA, n, n, dimnames = list(cols, cols))
    for (i in 1:1) {
      for (j in 2:3) {
        cor_mat[i, j] <- if (i == j) 1 else cor(rank(dat.pval.wide[, i]),
                                                rank(dat.pval.wide[, j]),
                                                method = cor_method)
      }
    }
    result$cor <- cor_mat
  }
  
  return(result)
}

methods <- c('C-SIDE','spVC','CELINA','STANCE')

cor.res <- expand.grid(dataset = datasets, method = methods, stringsAsFactors = FALSE) %>%
  pmap_dfr(function(dataset, method) {
    dat.pval.wide <- get_wide_pval(dataset, method)
    data <- calc_similarity(dat.pval.wide,cor = TRUE)
    
    as.data.frame(as.table(data$cor)) %>%
      na.omit() %>%
      rename(angle1 = Var1, angle2 = Var2, correlation = Freq) %>%
      mutate(dataset = dataset, method = method)
  })

jac.alpha.res <- expand.grid(dataset = datasets, method = methods, stringsAsFactors = FALSE) %>%
  pmap_dfr(function(dataset, method) {
    dat.pval.wide <- get_wide_pval(dataset, method)
    data <- calc_similarity(dat.pval.wide,jaccard = TRUE,jaccard_mode = "pval")
    
    as.data.frame(as.table(data$jaccard)) %>%
      na.omit() %>% 
      rename(angle1 = Var1, angle2 = Var2, jaccard = Freq) %>%
      mutate(dataset = dataset,method = method)
  })

jac.top50.res <- expand.grid(dataset = datasets, method = methods, stringsAsFactors = FALSE) %>%
  pmap_dfr(function(dataset, method) {
    dat.pval.wide <- get_wide_pval(dataset, method)
    data <- calc_similarity(dat.pval.wide,jaccard = TRUE,jaccard_mode = "top")
    
    as.data.frame(as.table(data$jaccard)) %>%
      na.omit() %>% 
      rename(angle1 = Var1, angle2 = Var2, jaccard = Freq) %>%
      mutate(dataset = dataset,method = method)
  })

head(cor.res)

# Combine results for visualization
plot.df <- bind_rows(
  jac.alpha.res %>% rename(value = jaccard) %>% mutate(metric = "Jaccard_alpha"),
  jac.top50.res %>% rename(value = jaccard) %>% mutate(metric = "Jaccard_top50"),
  cor.res %>% rename(value = correlation) %>% mutate(metric = "Correlation")
)

plot.df <- plot.df %>%
  mutate(method = recode(method,
                         "C-SIDE" = "C-SIDE",
                         "CELINA" = "Celina",
                         "STANCE" = "STANCE",
                         "spVC" = "spVC")) %>% 
  mutate(angle=angle2)

# Rank distribution per metric
rank.df <- plot.df %>%
  group_by(metric, angle1, angle2, dataset) %>%
  mutate(rank = rank(-value, ties.method = "min")) %>%  
  ungroup()

rank.count <- rank.df %>%
  group_by(metric, method, rank) %>%
  summarise(count = n(), .groups = "drop")

p2 <- ggplot(rank.count, aes(x = rank, y = method, size = count, color = rank)) +
  geom_point(alpha = 0.85) +
  facet_wrap(
    ~ metric,
    labeller = labeller(
      metric = function(x) {
        recode(x,
               "Correlation" = "Correlation",
               "Jaccard_alpha" = "Jaccard (α = 0.05)",
               "Jaccard_top50" = "Jaccard (top 50 genes)"
        )
      }
    )
  )+
  scale_size_continuous(range = c(0.2, 6)) +
  scale_color_gradient(low = "darkblue", high = "lightblue", guide = "none") +
  theme_minimal() +
  my_theme +
  labs(x="Rank",y="Methods",size="Count")+
  theme(
    plot.margin = margin(0, 5, 0, 0),
    legend.key.size = unit(0.05, "in"), 
    legend.key.width = unit(0.05, "in"),   
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
p2

plot_grid(plotlist = list(p1,p2), nrow = 2, rel_heights = c(1,1), labels = c("A", "B"))
ggsave('./Fig/Fig4.pdf', width = 6.69, height = 4)

summary_df_rob1 <- plot.df %>%
  group_by(method,metric) %>%
  summarise(values = mean(value, na.rm = TRUE)) %>%
  pivot_wider(
    names_from = c(metric),    
    values_from = values,         
    names_glue = "{metric}") %>% 
  column_to_rownames('method')

summary_df <- read.csv("metrics_summary.csv", row.names = 1)
summary_df <- summary_df %>%
  mutate(
    Correlation   = summary_df_rob1[rownames(summary_df),"Correlation"],
    Jaccard_alpha = summary_df_rob1[rownames(summary_df),"Jaccard_alpha"],
    Jaccard_top50 = summary_df_rob1[rownames(summary_df),"Jaccard_top50"],
    Concordance   = summary_df_rob2[rownames(summary_df),"Concordance"]
  )

write.csv(summary_df, "metrics_summary.csv", row.names = TRUE)
