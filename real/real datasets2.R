Sys.setenv(OPENBLAS_NUM_THREADS = "1")
Sys.setenv(OMP_NUM_THREADS = "1")
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
source('/home/user/Fanglab1/yh/ctSVGbench/real/CTSV.R')

ncores=20

datasets <- c(
  "Visium_mousebrain",
  "ST_Developmentalheart",
  "Visium_intestine",
  "Visium_pancreas",
  "StereoSeq_MDESTA",
  "Slide-seqV2_melanoma",
  "Visium_lymph_node",
  "Slide-seq_tumor"
  "Slide-seqV2_hippocampus",
  "StereoSeq_CBMSTA_Macaque",
  "Slide-seqV2_mouseOB",
  "Slide-seqV2_melanoma_GSM6025944_MBM13"
)

for (dataset in datasets){
  library(RhpcBLASctl)
  blas_set_num_threads(1)
  file <- sprintf('myRCTD_%s.rds',dataset)
  puck<- readRDS(here('real','puck',file))
  reference <- readRDS(here('real','reference',file))
  # reference <- readRDS("/home/user/Fanglab1/yh/ctSVGbench/real/reference/myRCTD_Slide-seqV2_melanoma_MBM.rds")  
  pos <- readRDS(here('real','pos_subset',file))
  boundary <- readRDS(here('real','boundary',file))
  counts.orign <- puck@counts[,rownames(pos)]
  counts.sc <- reference@counts
  celltypes <- reference@cell_types
  
  mito_genes = unique(c(grep("^MT-", rownames(counts.orign)), grep("^mt-", rownames(counts.orign))))
  mat <- as(counts.orign, "dgCMatrix")
  low_genes = which(Matrix::rowSums(mat > 0)<20)
  remove_genes = unique(c(mito_genes, low_genes))
  if(length(remove_genes)>0){
    counts = mat[-remove_genes,]
  }else{
    counts = mat
  }
  
  ncell <- data.frame(table(celltypes))
  if(any(ncell$Freq<25)){
    filter_cell <- ncell$celltypes[ncell$Freq<25]    
    filter_index <- which(celltypes %in% filter_cell )
    celltypes <- celltypes[-filter_index]
    new_level <- unique(as.character(celltypes))
    celltypes <- factor(celltypes, levels = new_level)
    counts.sc <- counts.sc[,-filter_index]
  }
  reference <- Reference(counts = counts.sc, cell_types = as.factor(celltypes), min_UMI = 1) 
  puck <- SpatialRNA(coords = pos, counts = counts)
  myRCTD <- create.RCTD(puck, reference = reference, max_cores = ncores)
  
  # run RCTD
  if(grepl("Visium|ST", dataset)){
    myRCTD <- run.RCTD(myRCTD, doublet_mode = "full") 
    doublet_mode = F
  } else {
    myRCTD <- run.RCTD(myRCTD, doublet_mode = "doublet") 
    doublet_mode = T
  }
  
  prop <- as.matrix(normalize_weights(myRCTD@results$weights))
  saveRDS(prop,here('real','prop',file))
  prop <- readRDS(here('real','prop',file))
  pos=pos[rownames(prop),]
  counts=counts[,rownames(prop)]
  gc1 <- gc(reset = TRUE)
  time = system.time({
    CSIDE.results<- run.CSIDE.nonparam(myRCTD, 
                                       cell_type_threshold = 0,
                                       cell_types = names(tail(sort(colSums(prop)),3)), fdr = 0.05, doublet_mode = doublet_mode)})
  
  gc2 <- gc()
  res.cside=CSIDE.results@de_results$all_gene_list
  saveRDS(res.cside,here('real','res',sprintf('%s-C-SIDE.rds',dataset)))
  
  computation = data.frame(
    time = time[['elapsed']],
    n_spot = ncol(counts),
    n_gene = nrow(counts),
    dataset = dataset,
    Peak_RAM_Used_MiB = sum(gc2[,6] - gc1[,6]),
    method = "C-SIDE"
  )
  
  write.table(computation, 
              file =here('real','computation',sprintf('%s-C-SIDE.csv',dataset)),
              sep = ',',
              row.names = FALSE)
  
  
  ### Create Celina object ---
  gc1 <- gc(reset = TRUE)
  time = system.time({
    celltype_to_test <- names(tail(sort(colSums(prop)),3))
    normalized_counts <- scater::normalizeCounts(counts)
    Obj = Create_Celina_Object(celltype_mat = t(prop), 
                               gene_expression_mat = as.matrix(normalized_counts), 
                               location = as.matrix(pos),
                               covariates = NULL)
    Obj@celltype_mat <- as.matrix(Obj@celltype_mat)
    
    Obj = preprocess_input(Obj, 
                           cell_types_to_test = celltype_to_test ,  
                           scRNA_count = as.matrix(counts.sc), 
                           sc_cell_type_labels = as.matrix(data.frame(celltypes,row.names = colnames(counts.sc))))
    
    Obj = Calculate_Kernel(Obj)
    Obj = Testing_interaction_all(Obj, num_cores=ncores)
  })
  gc2 <- gc()
  res.celina <- Obj@result
  saveRDS(res.celina,here('real','res',sprintf('%s-CELINA.rds',dataset)))
  
  computation = data.frame(
    time = time[['elapsed']],
    n_spot = ncol(counts),
    n_gene = nrow(counts),
    dataset = dataset,
    Peak_RAM_Used_MiB = sum(gc2[,6] - gc1[,6]),
    method = "CELINA"
  )
  
  write.table(computation, 
              file =here('real','computation',sprintf('%s-CELINA.csv',dataset)),
              sep = ',',
              row.names = FALSE)
  
  
  #fit the spVC model
  gc1 <- gc(reset = TRUE)
  time = system.time({
    Tr.cell <- TriMesh(boundary, n = 2) # n : triangulation fineness
    V <- as.matrix(Tr.cell$V) 
    Tr <- as.matrix(Tr.cell$Tr)  
    # Fit the spVC models
    results <- test.spVC(Y = counts, X = prop, S = pos, V = V, Tr = Tr,
                         para.cores = ncores)})
  
  gc2 <- gc()
  
  saveRDS(results,here('real','res',sprintf('%s-spVC.rds',dataset)))
  
  computation = data.frame(
    time = time[['elapsed']],
    n_spot = ncol(counts),
    n_gene = nrow(counts),
    dataset = dataset,
    Peak_RAM_Used_MiB = sum(gc2[,6] - gc1[,6]),
    method = "spVC"
  )
  
  write.table(computation, 
              file =here('real','computation',sprintf('%s-spVC.csv',dataset)),
              sep = ',',
              row.names = FALSE)
  
  ### Create STANCE object ---
  ncores=70

  gc1 <- gc(reset = TRUE)
  time = system.time({
    Obj.STANCE<- creatSTANCEobject(counts = counts,
                                   pos = pos,
                                   prop = prop,
                                   covariates = NULL) #The rownames of covariates matrix should match the colnames of counts matrix.
    
    pred.STANCE <- data_preprocess(object = Obj.STANCE, 
                                   normalized = FALSE)
    
    mySTANCE <- build_kernelMatrix(object = pred.STANCE)
    mySTANCE <- runTest1(object = mySTANCE,
                         correction = F, pv.adjust = "BY")
    head(mySTANCE@Test_1)
    genes.list <- rownames(mySTANCE@gene_expression)
    utSVG.list <- genes.list[mySTANCE@Test_1$p_value < 0.05]
    
    mySTANCE <- runTest2(object = mySTANCE, 
                         Genes_to_test = utSVG.list, 
                         Cell_types_to_test = names(tail(sort(colSums(prop)),3)),
                         correction = F, 
                         ncores = ncores)
  })
  gc2 <- gc()
  res.stance <- mySTANCE@Test_2
  saveRDS(res.stance,here('real','res',sprintf('%s-STANCE.rds',dataset)))
  
  computation = data.frame(
    time = time[['elapsed']],
    n_spot = ncol(counts),
    n_gene = nrow(counts),
    dataset = dataset,
    Peak_RAM_Used_MiB = sum(gc2[,6] - gc1[,6]),
    method = "STANCE"
  )
  
  write.table(computation, 
              file =here('real','computation',sprintf('%s-STANCE.csv',dataset)),
              sep = ',',
              row.names = FALSE)

#fit the CTSV model 
  spe <- SpatialExperiment(assay = counts[,rownames(prop)], colData = pos[rownames(prop),], spatialCoordsNames = c('x', 'y')) 

  CTSV.results <- CTSV(spe, W = as.matrix(prop), num_core = ncores) 
  res.ctsv.matrix <- CTSV.results$pval
  res.ctsv <- setNames(
    lapply(seq_along(colnames(prop)), function(i){
       data.frame(pval = pmin(res.ctsv.matrix[,i], res.ctsv.matrix[,i+ncol(prop)]),
                 row.names = rownames(CTSV.results$qval))
    }
    ),
    colnames(prop)
  )
  ct_total <- colSums(prop)
  top3_ct <- head(names(sort(ct_total, decreasing = TRUE)),3)

  saveRDS(res.ctsv[top3_ct],here('real','res',sprintf('%s-CTSV.rds',dataset)))  

}                




