library(SpatialExperiment)
library(spacexr)
library(Matrix)
library(devtools)
library(ggplot2)
library(here)

# Function to extract a centered subset of spatial positions

# Function to extract a centered subset of spatial positions
get_centered_subset <- function(pos, 
                                         target_rows = 3000, 
                                         initial_radius_ratio = 0.3, 
                                         tol = 50, 
                                         max_attempts = 10000) {
  if (!all(c("x", "y") %in% colnames(pos))) {
    stop("colnames is not c('x','y')")
  }
  
  center_x <- median(pos$x)
  center_y <- median(pos$y)
  
  pos$distance_to_center <- sqrt((pos$x - center_x)^2 + (pos$y - center_y)^2)
  max_distance <- max(pos$distance_to_center)
  
  radius_ratio <- initial_radius_ratio
  step_size <- 0.00001
  current_rows <- 0
  attempts <- 0
  
  while (attempts < max_attempts) {
    current_radius <- radius_ratio * max_distance
    pos.subset <- pos[pos$distance_to_center <= current_radius, ]
    current_rows <- nrow(pos.subset)
    
    if (abs(current_rows - target_rows) <= tol) {
      break
    } else if (current_rows < target_rows) {
      radius_ratio <- radius_ratio + step_size
    } else {
      radius_ratio <- radius_ratio - step_size
    }
    
    if (radius_ratio < 0) {
      radius_ratio <- 0
      warning("fail")
      break
    }
    
    attempts <- attempts + 1
  }
  
  pos.subset$distance_to_center <- NULL
  original_rows <- nrow(pos)
  
  message(sprintf("original：%d，now：%d", original_rows, current_rows))
  
  if (attempts >= max_attempts) {
    warning(sprintf("exceed attempt times: %d", max_attempts))
  }
  
  return(pos.subset)
}
cell_sample_sizes <- c(100, 500, 1000, 2000, 5000, 8000, 10000, 20000, 50000, 100000)
for(target_rows in cell_sample_sizes){
  mylist <- readRDS('/home/user/Fanglab1/yh/ctSVGbench/real/input_sc/stereoseq_mosta_E16.5_E1S3.rds')
  pos <- mylist$pos
  pos.subset <- get_centered_subset(pos,target_rows = target_rows,initial_radius_ratio =round(target_rows/37363,4),  tol = 1, max_attempts = 10000)
}

# List of datasets
datasets <- c(
  "ST_PDAC",
  "Slide_seqV2_melanoma",
  "Visium_mousebrain",
  "ST_Developmentalheart",
  "Visium_liver",
  "StereoSeq_MDESTA",
  "StereoSeq_CBMSTA_Macaque",
  "SeqFish+_mouse_ob",
  "StereoSeq_CBMSTA_Marmoset",
  "Visium_spleen"
)

# process each dataset
lapply(datasets, function(dataset) {
  # read spatial coordinates and cell type proportions
  pos <- readRDS(file = here('real','pos', sprintf('myRCTD_%s.rds', dataset)))
  # prop <- readRDS(file = here('ctSVGbench','prop', sprintf('myRCTD_%s.rds', dataset)))
  
  # extract centered subset
  pos.subset <- get_centered_subset(pos)
  print(dataset)
  print(dim(pos))
  print(dim(pos.subset))
  # spot <- intersect(rownames(pos.subset), rownames(prop))
  
  # subset proportions
  # prop.subset <- prop[spot, ]
  
  # save subsetted data
  # saveRDS(prop.subset, file = here('ctSVGbench','prop', sprintf('myRCTD_%s.rds', dataset)))
  saveRDS(pos.subset, file = here('real','pos_subset', sprintf('myRCTD_%s.rds', dataset)))
})

