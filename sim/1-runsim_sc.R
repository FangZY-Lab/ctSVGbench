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

# List of datasets to analyze
datasets <- c(
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
  "stereoseq_mosta_E16.5_E1S3_whole_brain",
  "stereoseq_mosta_Dorsal_midbrain"
 )

# Source the analysis function
source(here('sim','utils',"run_analysis_for_pattern_sc.R"))

# Loop through each dataset
for (dataset in datasets) {
  # Read cell type proportions and spatial positions
  file.orign <- sprintf('myRCTD_%s.rds', dataset)
  prop.orign <- readRDS(here('sim','prop', file.orign))
  
  # Make column names syntactically valid
  colnames(prop.orign) <- make.names(colnames(prop.orign))
  
  # Select top cell types by abundance
  prop.use <- prop.orign[, names(tail(sort(colSums(prop.orign))))]
  prop.use <- prop.use[rowSums(prop.use)>0,]  
  
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
