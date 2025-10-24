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
setwd('~/Fanglab1/yh/ctSVGbench')
library(here)


dataset=  "Visium_skin"

for (i in c(1:100)){
  ncores=10
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
  
  mito_genes = unique(grep("^MT-", rownames(counts.orign),ignore.case = T))
  mat <- as(counts.orign, "dgCMatrix")
  low_genes = which(Matrix::rowSums(mat > 0)<20)
  remove_genes = unique(c(mito_genes, low_genes))
  if(length(remove_genes)>0){
    counts = mat[-remove_genes,]
  }else{
    counts = mat
  }
  set.seed(i)
  rownames(pos) <- sample(rownames(pos))
  puck <- SpatialRNA(coords=pos, counts=counts)

  reference <- Reference(counts = counts.sc, cell_types = as.factor(celltypes), min_UMI = 1) 
  myRCTD <- create.RCTD(puck, reference = reference, max_cores = ncores)
  myRCTD <- run.RCTD(myRCTD,  doublet_mode = 'full') 
  prop <- as.matrix(normalize_weights(myRCTD@results$weights))
  
  pos=pos[rownames(prop),]
  counts=counts[,rownames(prop)]
  
  CSIDE.results<- run.CSIDE.nonparam(myRCTD, 
                                     cell_type_threshold = 0,
                                     cell_types = names(tail(sort(colSums(prop)),3)), fdr = 0.05, doublet_mode = FALSE)
  
  res.cside=CSIDE.results@de_results$all_gene_list
  saveRDS(res.cside,here('real','res',sprintf('%s-C-SIDE-null%s.rds',dataset,i)))
  
  
  ### Create Celina object ---
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
  
  res.celina <- Obj@result
  saveRDS(res.celina,here('real','res',sprintf('%s-CELINA-null%s.rds',dataset,i)))
  
  #fit the spVC model
  boundary <- readRDS(here('real','boundary',file))
  Tr.cell <- TriMesh(boundary, n = 2) # n : triangulation fineness
  V <- as.matrix(Tr.cell$V) 
  Tr <- as.matrix(Tr.cell$Tr)  
  # Fit the spVC models
  results <- test.spVC(Y = counts, X = prop, S = pos, V = V, Tr = Tr,
                       para.cores = ncores)  
  saveRDS(results,here('real','res',sprintf('%s-spVC-null%s.rds',dataset,i)))
  
  ### Create STANCE object ---
  ncores=60
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
  
  res.stance <- mySTANCE@Test_2
  saveRDS(res.stance,here('real','res',sprintf('%s-STANCE-null%s.rds',dataset,i)))
}                




