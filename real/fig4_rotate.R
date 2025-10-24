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
setwd('~/yh')
library(here)

ncores=20

rotate_points <- function(original_points, angle_degrees) { 
  angle_radians <- angle_degrees * (pi / 180) # Convert degrees to radians 
  rotation_matrix <- matrix(c(cos(angle_radians), -sin(angle_radians), sin(angle_radians), cos(angle_radians)), nrow = 2, byrow = TRUE) 
  rotated_points <- t(rotation_matrix %*% t(original_points)) 
  output <- as.matrix(rotated_points) 
  colnames(output) <- c('x', 'y') 
  return(output) 
}

datasets <- c(
  "ST_PDAC",
  "Visium_spinal",
  "Slide-seqV2_melanoma",
  "Visium_skin",
  "SeqFish+_mouse_cortex_svz",
  "Visium_bladder",
  "Visium_tail",
  "Visium_liver",
  "SeqFish+_mouse_ob",
  "StereoSeq_mouseOB"
)
angles <- c(30,90)
for (dataset in datasets){
  for (angle_degrees in angles){
    
    file <- sprintf('myRCTD_%s.rds',dataset)
    puck<- readRDS(here('real','puck',file))
    reference <- readRDS(here('real','reference',file))
    pos.original <- puck@coords[1:2]
    pos <- rotate_points(pos.original, angle_degrees) 
    pos <- as.data.frame(pos)  
    counts.orign <- puck@counts
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
    saveRDS(prop,here('real','prop',file))
    # prop <- readRDS(here('real','prop',file))
    pos=pos[rownames(prop),]
    counts=counts[,rownames(prop)]
    
    gc1 <- gc(reset = TRUE)
    time = system.time({
      CSIDE.results<- run.CSIDE.nonparam(myRCTD, 
                                         cell_type_threshold = 1,
                                         cell_types = names(tail(sort(colSums(prop)),3)), fdr = 0.05, doublet_mode = doublet_mode)})
    
    gc2 <- gc()
    res.cside=CSIDE.results@de_results$all_gene_list
    saveRDS(res.cside,here('real','res',sprintf('%s-r%s-C-SIDE.rds',dataset,angle_degrees)))
    
    computation = data.frame(
      time = time[['elapsed']],
      n_spot = ncol(counts),
      n_gene = nrow(counts),
      dataset = dataset,
      Peak_RAM_Used_MiB = sum(gc2[,6] - gc1[,6]),
      method = "C-SIDE"
    )
    
    write.table(computation, 
                file =here('real','computation',sprintf('%s-r%s-C-SIDE.csv',dataset,angle_degrees)),
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
      
      for (each_celltype in celltype_to_test){
        Obj@genes_list[[each_celltype]] <- rownames(counts)
      }
      
      Obj = preprocess_input(Obj, 
                             cell_types_to_test = celltype_to_test ,  
                             scRNA_count = as.matrix(counts.sc), 
                             sc_cell_type_labels = as.matrix(data.frame(celltypes,row.names = colnames(counts.sc))))
      
      Obj = Calculate_Kernel(Obj)
      Obj = Testing_interaction_all(Obj, num_cores=ncores)
    })
    gc2 <- gc()
    res.celina <- Obj@result
    saveRDS(res.celina,here('real','res',sprintf('%s-r%s-CELINA.rds',dataset,angle_degrees)))
    
    computation = data.frame(
      time = time[['elapsed']],
      n_spot = ncol(counts),
      n_gene = nrow(counts),
      dataset = dataset,
      Peak_RAM_Used_MiB = sum(gc2[,6] - gc1[,6]),
      method = "CELINA"
    )
    
    write.table(computation, 
                file =here('real','computation',sprintf('%s-r%s-CELINA.csv',dataset,angle_degrees)),
                sep = ',',
                row.names = FALSE)
    
    
    #Sys.setenv(OPENBLAS_NUM_THREADS = "96")
    #fit the spVC model
    gc1 <- gc(reset = TRUE)
    time = system.time({
      boundary <- readRDS(here('real','boundary',file))
      boundary.r <- rotate_points(boundary, angle_degrees = angle_degrees)
      Tr.cell <- TriMesh(boundary.r, n = 2) # n : triangulation fineness
      V <- as.matrix(Tr.cell$V) 
      Tr <- as.matrix(Tr.cell$Tr)  
      # Fit the spVC models
      results <- test.spVC(Y = counts, X = prop, S = pos, V = V, Tr = Tr,
                           para.cores = ncores)})
    
    gc2 <- gc()
    
    saveRDS(results,here('real','res',sprintf('%s-r%s-spVC.rds',dataset,angle_degrees)))
    
    computation = data.frame(
      time = time[['elapsed']],
      n_spot = ncol(counts),
      n_gene = nrow(counts),
      dataset = dataset,
      Peak_RAM_Used_MiB = sum(gc2[,6] - gc1[,6]),
      method = "spVC"
    )
    
    write.table(computation, 
                file =here('real','computation',sprintf('%s-r%s-spVC.csv',dataset,angle_degrees)),
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
    saveRDS(res.stance,here('real','res',sprintf('%s-r%s-STANCE.rds',dataset,angle_degrees)))
    
    computation = data.frame(
      time = time[['elapsed']],
      n_spot = ncol(counts),
      n_gene = nrow(counts),
      dataset = dataset,
      Peak_RAM_Used_MiB = sum(gc2[,6] - gc1[,6]),
      method = "STANCE"
    )
    
    write.table(computation, 
                file =here('real','computation',sprintf('%s-r%s-STANCE.csv',dataset,angle_degrees)),
                sep = ',',
                row.names = FALSE)
  }                
  
}




