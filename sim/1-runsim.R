# Set single-threading for BLAS/OMP to avoid multithreading conflicts
Sys.setenv(OPENBLAS_NUM_THREADS = "1")
Sys.setenv(OMP_NUM_THREADS = "1")

setwd('~/Fanglab1/yh/ctSVGbench')

library(RhpcBLASctl)
blas_set_num_threads(1)
library(here)
library(parallel)

# Create output directories if they don't exist
dir.create("sim/computation", recursive = TRUE, showWarnings = FALSE)
dir.create("sim/res", recursive = TRUE, showWarnings = FALSE)
dir.create("sim/plot", recursive = TRUE, showWarnings = FALSE)

# List of datasets to analyze
datasets <- c(
  "ST_PDAC",
  "Visium_liver",
  "Visium_mousebrain",
  "Visium_spleen",
  "StereoSeq_MDESTA",
  "StereoSeq_CBMSTA_Macaque",
  "Slide-seqV2_melanoma",
  "SeqFish+_mouse_ob",
  "StereoSeq_CBMSTA_Marmoset",
  "Slide-seq_tumor"
)

# Source the analysis function
source(here('sim','utils',"run_analysis_for_pattern.R"))

# Loop through each dataset
for (dataset in datasets) {
  # Read cell type proportions and spatial positions
  file.orign <- sprintf('myRCTD_%s.rds', dataset)
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
  spot <- intersect(rownames(pos.use), rownames(prop.use))
  pos.use <- pos.use[spot, ]
  prop.use <- prop.use[spot, ]
  boundary <- readRDS(here('sim','boundary', file.orign))
  
  # Define spatial patterns to simulate
  patterns <- c("hotspot", "stripe", "pathology", "neighbor", "periodic", "gradient")
  
  # Loop through patterns, replicates, and parameter sets
  for (pt in patterns) {
    for (rep in 1:1) {
      for (paramset in c("P1", "P2", "P3")) {
        # Reset single-threading before each run
        Sys.setenv(OPENBLAS_NUM_THREADS = "1")
        Sys.setenv(OMP_NUM_THREADS = "1")
        library(RhpcBLASctl)
        blas_set_num_threads(1)
        
        # Run the analysis
        run_analysis_for_pattern(
          pt, pos.use, prop.use, dt = dataset, boundary = boundary, 
          rep_id = rep, paramset = paramset, ncores = 70
        )
      }
    }
  }
}
