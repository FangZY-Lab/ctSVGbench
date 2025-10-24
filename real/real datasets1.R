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

ncores=96
datasets <- c(
  "MERFISH_hypothalamus",
  "SeqFish+_mouse_ob",
  "ST_PDAC",
  "StereoSeq_CBMSTA_Marmoset",
  "Visium_liver",
  "Visium_skin",
  "Visium_spleen",
  "Slide-seqV2_melanoma",
  "Visium_bladder",
  "Visium_tail",
  "StereoSeq_mouseOB"
)


for (dataset in datasets){
  library(RhpcBLASctl)
  blas_set_num_threads(1)
  file <- sprintf('myRCTD_%s.rds',dataset)
  puck<- readRDS(here('real','puck',file))
  reference <- readRDS(here('real','reference',file))
  pos <- puck@coords[1:2]
  counts.orign <- puck@counts
  counts.sc <- reference@counts
  celltypes <- reference@cell_types
  
  ncell <- data.frame(table(celltypes))
  if(any(ncell$Freq<100)){
    filter_cell <- ncell$celltypes[ncell$Freq<100]    
    filter_index <- which(celltypes %in% filter_cell )
    celltypes <- celltypes[-filter_index]
    new_level <- unique(as.character(celltypes))
    celltypes <- factor(celltypes, levels = new_level)
    counts.sc <- counts.sc[,-filter_index]
  }
  
  mito_genes = unique(c(grep("^MT-", rownames(counts.orign)), grep("^mt-", rownames(counts.orign))))
  mat <- as(counts.orign, "dgCMatrix")
  low_genes = which(Matrix::rowSums(mat > 0)<20)
  remove_genes = unique(c(mito_genes, low_genes))
  if(length(remove_genes)>0){
    counts = mat[-remove_genes,]
  }else{
    counts = mat
  }
  
  puck <- SpatialRNA(coords = pos, counts = counts)
  reference <- Reference(counts = counts.sc, cell_types = as.factor(celltypes), min_UMI = 1) 
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
  # prop <- readRDS(here('real','prop',file))
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
    boundary <- readRDS(here('real','boundary',file))
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
}                




