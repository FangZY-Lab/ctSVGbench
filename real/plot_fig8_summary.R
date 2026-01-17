library(funkyheatmap)
library(tidyverse)
library(purrr)

# ---- 1. input ----
summary_df <- read.csv("metrics_summary.csv", row.names = 1)
summary_df$Speed_hours <- NULL
summary_df <- summary_df[c("spVC", "Celina", "STANCE","CTSV","C-SIDE","ctSVG"),]

data <- summary_df %>%
  round(2) %>%
  rownames_to_column("id")
# data$minority <- c("yes","yes","no","no")
# data$majority <- c("yes","yes","yes","yes")

# ---- 2. group ----
metrics_list <- list(
  "AUC" = c("gradient_auc", "hotspot_auc", "neighbor_auc", "pathology_auc",
            "periodic_auc", "stripe_auc"),
  "Sensitivity" = c( "gradient_sensitivity", "hotspot_sensitivity",
                     "neighbor_sensitivity", "pathology_sensitivity",
                     "periodic_sensitivity", "stripe_sensitivity"),
  "Specificity"=c("specificity","affected_specificity_gradient", 
                  "affected_specificity_hotspot", 
                  "affected_specificity_neighbor", 
                  "affected_specificity_pathology", 
                  "affected_specificity_periodic", 
                  "affected_specificity_stripe"),
  "Robustness" = c("Jaccard_alpha", "Jaccard_top50","Correlation", "Concordance"),
  "Scalability" = c("CPU_time_5000", "RAM_GB_5000")
)

metrics_list <- lapply(metrics_list, function(cols) cols[cols %in% colnames(data)])

# ---- 3. column_info ----
cinfo <- tibble(
  id = c("id", unlist(metrics_list, use.names = FALSE)),
  group = c(NA, rep(names(metrics_list), lengths(metrics_list))),
  geom = c(
    "text",
    rep("rect", length(metrics_list$AUC)),
    rep("rect", length(metrics_list$Sensitivity)),
    rep("circle", 1),
    rep("circle",6),  
    rep("circle", length(metrics_list$Robustness)),
    rep("bar", length(metrics_list$Scalability))
  ),
  palette = c(
    NA,
    rep("pred_palette", length(metrics_list$AUC)),
    rep("sens_palette", length(metrics_list$Sensitivity)),
    rep("spec_palette", length(metrics_list$Specificity)),
    rep("rob_palette", length(metrics_list$Robustness)),
    rep("sca_palette", length(metrics_list$Scalability))
  ),
  options = vector("list", 1 + sum(lengths(metrics_list)))
)
cinfo$width <- rep(NA, ncol(data))
cinfo$width[c(25, 26)] <- 1.8
# ---- 4. palettes ----
palettes <- list(
  pred_palette = RColorBrewer::brewer.pal(9, "Blues")[4:2],
  spec_palette = RColorBrewer::brewer.pal(9, "Greens")[9:2],
  sens_palette = RColorBrewer::brewer.pal(9, "Purples")[5:2],
  rob_palette = RColorBrewer::brewer.pal(9, "YlOrBr")[9:2],
  sca_palette = RColorBrewer::brewer.pal(9, "Reds")[-8:-9]
)
palettes$mypalette <- RColorBrewer::brewer.pal(9, "Greys")[9:2]
palettes$black <- c(rep("black", 2))
# ---- 5. column_groups ----
column_groups <- tibble(
  MainCategory = c("Predictive performance", "Predictive performance","Predictive performance", "Robustness", "Scalability"),
  SubCategory  = c("AUC", "Specificity","Sensitivity", "Rotation invariance", "Compute cost"),
  group        = c("AUC", "Specificity","Sensitivity", "Robustness", "Scalability"),
  palette      = c("pred_palette", "pred_palette", "pred_palette","rob_palette", "sca_palette")
)

# ---- 6. legends ----
legends <-  list( 
  list(
    palette = "sens_palette",
    geom = "bar",
    enabled = FALSE
  ),  
  list(
    palette = "spec_palette",
    geom = "bar",
    enabled = FALSE
  ),
  list(
    palette = "pred_palette",
    geom = "bar",
    enabled = FALSE
  ),
  list(
    palette = "rob_palette",
    geom = "bar",
    enabled = FALSE
  ),
  list(
    palette = "mypalette",
    geom = "bar",
    title = "Scaled score",
    enabled = TRUE
  ),
  list(
    palette = "sca_palette",
    geom = "bar",
    title = "Scalability(scaled)",
    enabled = FALSE
  ))

label_map <- c(
  "gradient_auc" = "Gradient",
  "hotspot_auc" = "Hotspot",
  "neighbor_auc" = "Neighbor",
  "pathology_auc" = "Pathology",
  "periodic_auc" = "Periodic",
  "stripe_auc" = "Stripe",
  "specificity" = "permutation",
  "affected_specificity_gradient"="affected_gradient", 
  "affected_specificity_hotspot"="affected_hotspot",
  "affected_specificity_neighbor"="affected_neighbor",  
  "affected_specificity_pathology"="affected_pathology",  
  "affected_specificity_periodic"="affected_periodic",  
  "affected_specificity_stripe"="affected_stripe",
  "gradient_sensitivity" = "Gradient",
  "hotspot_sensitivity" = "Hotspot",
  "neighbor_sensitivity" = "Neighbor",
  "pathology_sensitivity" = "Pathology",
  "periodic_sensitivity" = "Periodic",
  "stripe_sensitivity" = "Stripe",
  "Correlation" = "Spearman",
  "Concordance" = "CAT",
  "Jaccard_alpha" = "Jaccard sig",
  "Jaccard_top50" = "Jaccard top",
  # "CPU_time_2000" = "CPU time_2k",
  # "RAM_GB_2000" = "Memory_2k",
  "CPU_time_5000" = "CPU time_5k",
  "RAM_GB_5000" = "Memory_5k"  
)
# cinfo[cinfo$id %in% c("minority","majority"), "directory"] <- "images"
# cinfo[cinfo$id %in% c("minority","majority"), "extension"] <- "png"
cinfo$name <- ifelse(cinfo$id %in% names(label_map),
                     label_map[cinfo$id],
                     cinfo$id)
# cinfo$name <- cinfo$id
# # ---- 7. scale ----

mm <- function(col) {
  col / max(col, na.rm = TRUE)
}

data$RAM_GB_5000 <- log2(data$RAM_GB_5000+1)
data$CPU_time_5000 <- log2(data$CPU_time_5000+1)
data[, 2:26] <- lapply(data[, 2:26], mm)

# ---- 8. funky heatmap ----
funky_heatmap(
  data,
  column_info = cinfo,
  column_groups = column_groups,
  palettes = palettes,
  legends = legends,
  scale_column = FALSE,
  position_args = 
    position_arguments(col_annot_offset = 2.5,
                       col_annot_angle =30))+theme_bw(base_size = 6)

ggsave('Fig/Fig8.2.pdf', width = 12, height = 5)

