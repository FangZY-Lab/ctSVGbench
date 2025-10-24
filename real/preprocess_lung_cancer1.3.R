library(Seurat)
library(ggplot2)
library(dplyr)
library(arrow)
library(here)
library(Seurat)
library(Matrix)
library(arrow)  
library(jsonlite)
library(ggplot2)
library(spacexr) #C-SIDE
library(SpatialExperiment)
library(tidyverse)
library(here)
library(data.table)
library(Matrix)
library(Seurat)
library(SummarizedExperiment)
library(data.table)
library(org.Hs.eg.db)   
library(AnnotationDbi)

process_sample <- function(sample_id, data_dir, row_range, col_range, cancer_type) {
  counts <- Read10X(file.path(data_dir, "filtered_feature_bc_matrix"))
  pos_file <- file.path(data_dir, "spatial/tissue_positions.parquet")
  tissue_pos <- read_parquet(pos_file)
  
  pos <- tissue_pos %>% filter(in_tissue == 1)
  p1 <- ggplot(pos, aes(x = array_row, y = array_col)) +
    geom_point()
  ggsave(paste0(sample_id, "_all_spots.pdf"), p1)
  
  selected_spots <- tissue_pos %>%
    filter(in_tissue == 1,
           array_row >= row_range[1] & array_row <= row_range[2],
           array_col >= col_range[1] & array_col <= col_range[2])
  message(sample_id, ": selected spots = ", nrow(selected_spots))
  
  p2 <- ggplot(selected_spots, aes(x = array_row, y = array_col)) +
    geom_point()
  ggsave(paste0(sample_id, "_selected_spots.pdf"), p2)
  
  pos_df <- data.frame(
    x = selected_spots$array_row,
    y = selected_spots$array_col,
    row.names = selected_spots$barcode
  )
  puck <- SpatialRNA(
    coords = pos_df,
    counts = counts[, intersect(colnames(counts), rownames(pos_df))]
  )
  
  saveRDS(puck, here("real", "puck",
                     paste0("myRCTD_VisiumHD_", cancer_type, "_", substr(sample_id,1,11), ".rds")))
  
}

process_sample("AXB-6976-A1",
               "/home/user/Fanglab1/yh/ctSVGbench/st_lusc+luad/AXB-6976-A1/square_008um",
               row_range = c(400, 580), col_range = c(150, 200), cancer_type = "LUAD")

process_sample("AXB-2431-A1",
               "/home/user/Fanglab1/yh/ctSVGbench/st_lusc+luad/AXB-2431-A1/square_008um",
               row_range = c(400, 500), col_range = c(400, 475), cancer_type = "LUAD")

process_sample("AXB-6123-A1",
               "/home/user/Fanglab1/yh/ctSVGbench/st_lusc+luad/AXB-6123-A1/square_008um",
               row_range = c(250, 430), col_range = c(650, 700), cancer_type = "LUAD")

process_sample("AXB-5488-D1",
               "/home/user/Fanglab1/yh/ctSVGbench/st_lusc+luad/AXB-5488-D1/square_008um",
               row_range = c(400, 500), col_range = c(380, 500), cancer_type = "LUSC")

process_sample("AXB-7437-D1",
               "/home/user/Fanglab1/yh/ctSVGbench/st_lusc+luad/AXB-7437-D1/square_008um",
               row_range = c(200, 520), col_range = c(100, 200), cancer_type = "LUSC")

process_sample("AXB-7941-D1",
               "/home/user/Fanglab1/yh/ctSVGbench/st_lusc+luad/AXB-7941-D1/square_008um",
               row_range = c(520, 800), col_range = c(475, 510), cancer_type = "LUSC")

