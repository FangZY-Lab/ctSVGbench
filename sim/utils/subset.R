library(SpatialExperiment)
library(spacexr)
library(Matrix)
library(devtools)
library(ggplot2)
library(here)

# Function to extract a centered subset of spatial positions
get_centered_subset <- function(pos, target_rows = 0.5 * nrow(pos), initial_range_ratio = 0.7, tol = 2) {
  # compute total range of x and y
  x_range <- max(pos$x) - min(pos$x)
  y_range <- max(pos$y) - min(pos$y)
  
  # compute center (median to reduce outlier effects)
  center_x <- median(pos$x)
  center_y <- median(pos$y)
  
  # dynamically adjust range_ratio to reach target_rows
  range_ratio <- initial_range_ratio
  step_size <- 0.005
  current_rows <- 0
  attempts <- 0
  max_attempts <- 10000
  
  while (attempts < max_attempts) {
    # compute current width and height
    subset_width <- range_ratio * x_range
    subset_height <- range_ratio * y_range
    
    # extract subset within the current window
    pos.subset <- pos[
      pos$x >= (center_x - subset_width / 2) & 
      pos$x <= (center_x + subset_width / 2) & 
      pos$y >= (center_y - subset_height / 2) & 
      pos$y <= (center_y + subset_height / 2), 
    ]
    
    current_rows <- nrow(pos.subset)
    
    # check if subset size is within tolerance
    if (abs(current_rows - target_rows) <= tol) {
      break
    } else if (current_rows < target_rows) {
      # expand range
      range_ratio <- range_ratio + step_size
    } else {
      # shrink range
      range_ratio <- range_ratio - step_size
    }
    
    attempts <- attempts + 1
  }
  
  # optional debug message
  message(sprintf("Final range_ratio: %.3f, Rows: %d", range_ratio, current_rows))
  
  return(pos.subset)
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
  pos <- readRDS(file = here('ctSVGbench','pos', sprintf('myRCTD_%s.rds', dataset)))
  prop <- readRDS(file = here('ctSVGbench','prop', sprintf('myRCTD_%s.rds', dataset)))
  
  # extract centered subset
  pos.subset <- get_centered_subset(pos)
  spot <- intersect(rownames(pos.subset), rownames(prop))
  
  # subset proportions
  prop.subset <- prop[spot, ]
  
  # save subsetted data
  saveRDS(prop.subset, file = here('ctSVGbench','prop', sprintf('myRCTD_%s_0.5s.rds', dataset)))
  saveRDS(pos.subset, file = here('ctSVGbench','pos', sprintf('myRCTD_%s_0.5s.rds', dataset)))
})
