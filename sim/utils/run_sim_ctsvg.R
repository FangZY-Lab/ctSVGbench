library(RhpcBLASctl)
blas_set_num_threads(1)
library(spacexr) #C-SIDE
library(spVC)
library(sp)
library(BPST)
library(Triangulation)
library(MGLM)
library(CELINA)
library(STANCE)
library(SpatialExperiment)
library(spacexr)
library(Matrix)
library(devtools)
library(reshape2)
library(dplyr)
library(here)
library(ctSVG)
library(CTSV)
library(SpatialExperiment)
library(here)
library(pscl)
library(qvalue)
source('/home/user/Fanglab1/yh/ctSVGbench/real/CTSV.R')
run_analysis_for_pattern <- function(pt,  
                                     pos.use, prop.use, dt = dataset, boundary, 
                                     rep_id = 1, paramset='P1',
                                     ncores = 10, sim_name=sim_name) {
  source(here('sim','utils','generate_sc.R'))
  st_code_path <- file.path("./sim/utils", paste0("generate_st_", paramset, ".R"))
  if (!file.exists(st_code_path)) {
    stop(paste("Script not found for paramset:", paramset))
  }
  source(st_code_path)
  message("Running analysis with ", st_code_path)    
  
  # normalize proportions
  prop.use <- sweep(prop.use, 1, rowSums(prop.use), '/')
  
  # generate spatial data and reference scRNA-seq
  stlist <- generate_spatial_data(pos = pos.use, cell_prop = prop.use, boundary = boundary, pattern = pt, seed = rep_id)
  refer.sc <- generate_sc(seed = rep_id)
  
  counts.sc <- refer.sc$expr_mat
  celltypes <- refer.sc$celltypes
  
  pos <- stlist$pos
  counts.orign <- stlist$counts
  mat <- as(counts.orign, "dgCMatrix")
  
  # filter low-expression genes
  low_genes <- which(Matrix::rowSums(mat > 0) < 20)
  remove_genes <- unique(low_genes)
  if(length(remove_genes) > 0){
    counts <- mat[-remove_genes,]
  } else {
    counts <- mat
  }
  
  library(RhpcBLASctl)
  blas_set_num_threads(1)
  
  
  # create SpatialRNA and reference objects
  puck <- SpatialRNA(coords = pos, counts = counts)
  reference <- Reference(counts = counts.sc, cell_types = as.factor(celltypes), min_UMI = 1) 
  myRCTD <- create.RCTD(puck, reference = reference, max_cores = ncores)
  
  # run RCTD
  if(grepl("Visium|ST", dt)){
    myRCTD <- run.RCTD(myRCTD, doublet_mode = "full") 
    doublet_mode = F
  } else {
    myRCTD <- run.RCTD(myRCTD, doublet_mode = "doublet") 
    doublet_mode = T
  }
  
  prop <- as.matrix(normalize_weights(myRCTD@results$weights))
  
  # ensure order matches
  pos <- pos[rownames(prop),]
  counts <- counts[, rownames(prop)]
  dims = as.data.frame(nrow(counts))
  write.csv(dims, file = here('sim','ngene_nofiltered', sprintf('%s.csv', sim_name)))
  
  # save RCTD proportions
  saveRDS(prop, here('sim','prop', paste0(sim_name, '.rds')))
  
  # C-SIDE
  
  CSIDE.results <- run.CSIDE.nonparam(myRCTD, 
                                      cell_type_threshold = 0,
                                      cell_types = names(tail(sort(colSums(prop)),3)), 
                                      fdr = 0.05, 
                                      doublet_mode = doublet_mode)
  
  res.cside <- CSIDE.results@de_results$all_gene_list
  saveRDS(res.cside, here('sim','res', sprintf('%s-C-SIDE.rds', sim_name)))
  
  
  # CELINA
  
  celltype_to_test <- names(tail(sort(colSums(prop)),3))
  normalized_counts <- scater::normalizeCounts(counts)
  Obj <- Create_Celina_Object(celltype_mat = t(prop), 
                              gene_expression_mat = as.matrix(normalized_counts), 
                              location = as.matrix(pos),
                              covariates = NULL)
  Obj@celltype_mat <- as.matrix(Obj@celltype_mat)
  
  Obj <- preprocess_input(Obj, 
                          cell_types_to_test = celltype_to_test,  
                          scRNA_count = as.matrix(counts.sc), 
                          sc_cell_type_labels = as.matrix(data.frame(celltypes, row.names = colnames(counts.sc))))
  
  Obj <- Calculate_Kernel(Obj)
  Obj <- Testing_interaction_all(Obj, num_cores = ncores)
  res.celina <- Obj@result
  saveRDS(res.celina, here('sim','res', sprintf('%s-CELINA.rds', sim_name)))
  
  
  # spVC
  Tr.cell <- TriMesh(boundary, n = 2)
  V <- as.matrix(Tr.cell$V) 
  Tr <- as.matrix(Tr.cell$Tr)  
  results <- test.spVC(Y = counts, X = prop, S = pos, V = V, Tr = Tr,
                       para.cores = ncores)
  
  saveRDS(results, here('sim','res', sprintf('%s-spVC.rds', sim_name)))
  
  
  # STANCE
  Obj.STANCE <- creatSTANCEobject(counts = counts,
                                  pos = pos,
                                  prop = prop,
                                  covariates = NULL)
  pred.STANCE <- data_preprocess(object = Obj.STANCE, normalized = FALSE)
  
  mySTANCE <- build_kernelMatrix(object = pred.STANCE)
  mySTANCE <- runTest1(object = mySTANCE, correction = F, pv.adjust = "BY")
  genes.list <- rownames(mySTANCE@gene_expression)
  utSVG.list <- genes.list[mySTANCE@Test_1$p_value < 0.05]
  
  mySTANCE <- runTest2(object = mySTANCE, 
                       Genes_to_test = utSVG.list, 
                       Cell_types_to_test = names(tail(sort(colSums(prop)),3)),
                       correction = F, ncores = ncores)
  
  res.stance <- mySTANCE@Test_2
  saveRDS(res.stance, here('sim','res', sprintf('%s-STANCE.rds', sim_name)))
  
  
  
  
  #fit the CTSV model
  spe <- SpatialExperiment(assay = counts[,which(colSums(counts) != 0)], colData = pos, spatialCoordsNames = c('x', 'y')) 
  CTSV.results <- CTSV(spe, W = prop, num_core = ncores)  
  top6_ct <- colnames(prop)
  res.ctsv.matrix <- CTSV.results$pval
  res.ctsv <- setNames(
    lapply(seq_along(top6_ct), function(i){
      data.frame(pval = pmin(res.ctsv.matrix[,i], res.ctsv.matrix[,i+length(top6_ct)]),
                 row.names = rownames(CTSV.results$qval))
    }
    ),
    top6_ct
  )
  top3_ct <- names(tail(sort(colSums(prop)),3))
  saveRDS(res.ctsv[top3_ct],here('sim','res',sprintf('%s-CTSV.rds',sim_name)))  
  cat(sprintf("Finished pattern %s", pt))
  
  
}
