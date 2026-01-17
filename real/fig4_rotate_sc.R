Sys.setenv(OPENBLAS_NUM_THREADS = "1")
Sys.setenv(OMP_NUM_THREADS = "1")
library(RhpcBLASctl)
blas_set_num_threads(1)
setwd('/home/user/Fanglab1/yh/ctSVGbench/')
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
library(ctSVG)
library(Seurat)
library(Matrix)
library(arrow)  
library(jsonlite)
source('/home/user/Fanglab1/yh/ctSVGbench/real/CTSV.R')
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
  "StereoSeq_CBMSTA_Marmoset1_T514",
  "StereoSeq_CBMSTA_Mouse1_T189",
  "StereoSeq_CBMSTA_Mouse2_T349",
  "MERFISH_hypothalamus",
  "SeqFish+_cortex"  
)

angles <- c(0,30,90)
for (dataset in datasets){
  for (angle_degrees in angles){
    ncores=20    
    file <- sprintf('myRCTD_%s.rds',dataset)
    puck<- readRDS(here('real','puck',file))
    pos.original <- readRDS(here('real','pos_r',file))
    pos <- rotate_points(pos.original, angle_degrees) 
    pos <- as.data.frame(pos) 
    counts.orign <- puck@counts
    
    pos <- pos[intersect(rownames(pos),colnames(counts.orign)),] 
    counts.orign <- counts.orign[,intersect(rownames(pos),colnames(counts.orign))] 
    
    mito_genes = unique(c(grep("^MT-", rownames(counts.orign)), grep("^mt-", rownames(counts.orign))))
    mat <- as(counts.orign, "dgCMatrix")
    low_genes = which(Matrix::rowSums(mat > 0)<20)
    remove_genes = unique(c(mito_genes, low_genes))
    if(length(remove_genes)>0){
      counts = mat[-remove_genes,]
    }else{
      counts = mat
    }
    prop <- readRDS(here('real','prop',file))
    
    spot <- intersect(colnames(counts),rownames(prop))
    prop <- prop[spot,]
    ct_total <- colSums(prop)
    top3_ct <- names(sort(ct_total, decreasing = TRUE))[1:3]
    prop_top3 <- prop[, top3_ct, drop = FALSE]
    
    prop_top3 <- prop_top3[rowSums(prop_top3)>0,]
    counts <- counts[,rownames(prop_top3)]
    pos <- pos[rownames(prop_top3),]  
    
    counts <- as.matrix(counts)
    pos <- as.matrix(pos)
    prop_top3 <- as.matrix(prop_top3) 
    cell_type <- apply(prop_top3, 1, function(row) {
      return(colnames(prop_top3)[which.max(row)])
    })
    counts <- counts[Matrix::rowSums(counts != 0) > 0, ]  
    
    ### Create C-SIDE object ---
    puck <- SpatialRNA(coords = as.data.frame(pos), counts = counts)
    reference <- Reference(counts = counts, cell_types = factor(cell_type), min_UMI = -Inf)
    ct_tab <- table(reference@cell_types)
    keep_ct <- names(ct_tab[ct_tab >= 10])
    cell_type_filter <- cell_type[cell_type %in% keep_ct]
    reference <- Reference(counts = counts[,names(cell_type_filter)], cell_types = factor(cell_type_filter))
    
    
    myRCTD <- create.RCTD(spatialRNA = puck, reference = reference, max_cores = ncores,
                          gene_cutoff = -Inf, fc_cutoff = -Inf, gene_cutoff_reg = -Inf, fc_cutoff_reg = -Inf, UMI_min = -Inf,
                          UMI_max = Inf, counts_MIN = -Inf, UMI_min_sigma = -Inf, CELL_MIN_INSTANCE = -Inf)
    myRCTD@config[["MIN_OBS"]] <- -Inf
    myRCTD@config[["MIN_CHANGE_BULK"]] <- -Inf
    myRCTD@config[["MIN_CHANGE_REG"]] <- -Inf
    myRCTD@config[["CONFIDENCE_THRESHOLD"]] <- -Inf
    # myRCTD <- run.RCTD(RCTD = myRCTD)
    myRCTD@config$RCTDmode <- "full"
    myRCTD <- import_weights(myRCTD = myRCTD, weights = prop_top3)
    
    # res_cside <- run.CSIDE.nonparam(myRCTD = myRCTD, cell_type_threshold = -Inf, gene_threshold = -Inf, doublet_mode = FALSE, fdr = Inf)
    
    CSIDE.results <- tryCatch(
      {
        # Core logic: run CSIDE.nonparam while suppressing all warnings
        suppressWarnings(
          run.CSIDE.nonparam(
            myRCTD = myRCTD,
            cell_type_threshold = 0,
            cell_types = names(tail(sort(colSums(prop_top3)), 3)),
            fdr = 0.05,
            doublet_mode = FALSE
          )
        )
      },
      # Catch all errors and handle them
      error = function(e) {
        # Print error message for debugging (optional)
        cat("Error caught: ", e$message, "\n", sep = "")
        # Return a default value to avoid downstream errors
        NULL
      }
    )
    if (!is.null(CSIDE.results)) {
      res.cside <- CSIDE.results@de_results$all_gene_list
    } else {
      res.cside <- NULL
    }
    
    saveRDS(res.cside,here('real','res',sprintf('%s-r%s-C-SIDE.rds',dataset,angle_degrees)))
    
    ### Create Celina object ---
    celltype_to_test <- names(tail(sort(colSums(prop_top3)),3))
    normalized_counts <- scater::normalizeCounts(counts)
    Obj = Create_Celina_Object(celltype_mat = t(prop_top3), 
                               gene_expression_mat = as.matrix(normalized_counts), 
                               location = as.matrix(pos),
                               covariates = NULL)
    Obj@celltype_mat <- as.matrix(Obj@celltype_mat)
    
    Obj = preprocess_input(Obj, 
                           cell_types_to_test = celltype_to_test ,  
                           scRNA_count = as.matrix(counts), 
                           sc_cell_type_labels = as.matrix(data.frame(cell_type,row.names = colnames(counts))))
    
    Obj = Calculate_Kernel(Obj)
    Obj = Testing_interaction_all(Obj, num_cores=ncores)
    
    res.celina <- Obj@result
    saveRDS(res.celina,here('real','res',sprintf('%s-r%s-CELINA.rds',dataset,angle_degrees)))
    
    ### Create STANCE object ---
    ncores=110
    
    Obj.STANCE<- creatSTANCEobject(counts = counts,
                                   pos = pos,
                                   prop = prop_top3,
                                   covariates = NULL) #The rownames of covariates matrix should match the colnames of counts matrix.
    
    pred.STANCE <- data_preprocess(object = Obj.STANCE, 
                                   normalized = FALSE)
    
    mySTANCE <- build_kernelMatrix(object = pred.STANCE)
    mySTANCE <- runTest1(object = mySTANCE,
                         correction = F, pv.adjust = "BY")
    head(mySTANCE@Test_1)
    genes.list <- rownames(mySTANCE@gene_expression)
    utSVG.list <- genes.list[mySTANCE@Test_1$p_value < 0.05]
    
    res.stance = tryCatch(
      {
        suppressWarnings(
          mySTANCE <- runTest2(
            object = mySTANCE, 
            Genes_to_test = utSVG.list, 
            Cell_types_to_test = names(tail(sort(colSums(prop_top3)), 3)),
            correction = FALSE,
            ncores = ncores
          )  
        )  
        
        mySTANCE@Test_2
      },
      error = function(e) {
        message(e$message)  
        NULL
      }
    )  
    saveRDS(res.stance,here('real','res',sprintf('%s-r%s-STANCE.rds',dataset,angle_degrees)))
    
    #fit the ctsvg model
    coord <- pos
    colnames(coord) <- c("row","col")
    cell_types <- apply(prop_top3, 1, function(row) {
      return(colnames(prop_top3)[which.max(row)])
    })  
    d <- CreateSeuratObject(counts = counts) 
    d <- AddMetaData(object = d, metadata = cell_types, col.name = "celltype")
    Idents(d) <- "celltype"
    d <- NormalizeData(d)
    
    res.ctsvg <- tryCatch(
      {
        suppressWarnings(
          ctsvg_test(d = d, coord = coord, recluRes = NULL)
        )
      },
      error = function(e) {
        message("ctsvg_test failed: ", e$message)
        NULL
      }
    )
    
    
    saveRDS(res.ctsvg,here('real','res',sprintf('%s-r%s-ctsvg.rds',dataset,angle_degrees)))
    
    
    #fit the spVC model
    boundary <- readRDS(here('real','boundary_r',file))
    boundary.r <- rotate_points(boundary, angle_degrees = angle_degrees)
    Tr.cell <- TriMesh(boundary.r, n = 2) # n : triangulation fineness
    V <- as.matrix(Tr.cell$V) 
    Tr <- as.matrix(Tr.cell$Tr)  
    
    # Fit the spVC models safely
    res.spvc <- tryCatch(
      {
        suppressWarnings(
          test.spVC(
            Y = counts,
            X = prop_top3,
            S = pos,
            V = V,
            Tr = Tr,
            para.cores = ncores
          )
        )
        
      },
      error = function(e) {
        message("test.spVC failed: ", e$message)
        NULL
      }
    )
    
    saveRDS(res.spvc,here('real','res',sprintf('%s-r%s-spVC.rds',dataset,angle_degrees)))
    
    #fit the CTSV model 
    ncores=70
    spe <- SpatialExperiment(assay = counts[,rownames(prop_top3)], colData = pos[rownames(prop_top3),], spatialCoordsNames = c('x', 'y')) 
    CTSV.results <- CTSV(spe, W = as.matrix(prop_top3), num_core = ncores) 
    res.ctsv.matrix <- CTSV.results$pval
    res.ctsv <- setNames(
      lapply(seq_along(top3_ct), function(i){
        data.frame(pval = pmin(res.ctsv.matrix[,i], res.ctsv.matrix[,i+length(top3_ct)]),
                   row.names = rownames(CTSV.results$qval))
      }
      ),
      top3_ct
    )
    saveRDS(res.ctsv,here('real','res',sprintf('%s-r%s-CTSV.rds',dataset,angle_degrees)))
  }                
  
}




