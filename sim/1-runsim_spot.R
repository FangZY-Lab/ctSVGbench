# Set single-threading for BLAS/OMP to avoid multithreading conflicts
Sys.setenv(OPENBLAS_NUM_THREADS = "1")
Sys.setenv(OMP_NUM_THREADS = "1")

setwd('~/Fanglab1/yh/ctSVGbench')

library(RhpcBLASctl)
blas_set_num_threads(1)
library(here)
library(parallel)

# Create output directories if they don't exist
dir.create("sim/res", recursive = TRUE, showWarnings = FALSE)

datasets <- c(
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
  colnames(pos.use) <- c("x","y")
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
          rep_id = rep, paramset = paramset, ncores = 80
        )
      }
    }
  }
}
