crop_spatial_by_size <- function(pos, target_size, center = c("median")) {
  center <- match.arg(center)
  
  pos <- as.matrix(pos)
  stopifnot(ncol(pos) == 2)
  ctr <- apply(pos, 2, median)
  
  dist2ctr <- sqrt((pos[,1] - ctr[1])^2 + (pos[,2] - ctr[2])^2)
  
  ord <- order(dist2ctr)
  
  if (target_size > nrow(pos)) stop("target size exceeds total spots")
  idx <- ord[seq_len(target_size)]
  pos[idx, , drop = FALSE]
}
library(here)
cell_sample_sizes <- c(500, 1000, 2000, 5000, 8000, 10000, 20000, 50000, 80000,100000)

mylist <- readRDS('/home/user/Fanglab1/yh/ctSVGbench/real/input_sc/stereoseq_mosta_E16.5_E1S3.rds')

pos.orign <- mylist$pos
for (target_size in cell_sample_sizes){
  pos <- crop_spatial_by_size(pos.orign, target_size)
  mat <- as.matrix(mylist$counts)[,rownames(pos)] 
  prop <- mylist$prop[rownames(pos),]
  prop <- prop[,colSums(prop)>0]
  print(colSums(prop))
  mito_genes = unique(c(grep("^MT-", rownames(mat)), grep("^mt-", rownames(mat))))
  low_genes = which(Matrix::rowSums(mat > 0)<35)
  remove_genes = unique(c(mito_genes, low_genes))
  if(length(remove_genes)>0){
    counts = mat[-remove_genes,]
  }else{
    counts = mat
  }
  set.seed(123)
  gene_detect <- Matrix::rowSums(counts > 0)
  top_genes <- order(gene_detect, decreasing = TRUE)[1:100]
  counts <- counts[top_genes, , drop = FALSE]
  
  # counts <- counts[sample(1:nrow(counts),100),]
  obj <- list(counts=counts,prop=prop,pos=pos)
  saveRDS(obj,file=here('real','input_sc',sprintf("stereoseq_mosta_E16.5_E1S3-%s-cell.rds",as.integer(target_size))))
}



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
library(Seurat)
library(Matrix)
library(arrow)  
library(jsonlite)
source('/home/user/Fanglab1/yh/ctSVGbench/real/CTSV.R')
ncores=1


cell_sample_sizes <- c(500, 1000, 2000, 5000, 8000, 10000, 20000, 50000, 80000,100000)

for (target_size in cell_sample_sizes){
  file <- sprintf("stereoseq_mosta_E16.5_E1S3-%s-cell.rds",as.integer(target_size))
  mylist <- readRDS(file = here('real','input_sc',file))
  pos <- mylist$pos
  saveRDS(pos,file = here('real','pos',file))
}

for (target_size in cell_sample_sizes){
  ncores=1
  library(RhpcBLASctl)
  blas_set_num_threads(1)
  file <- sprintf("stereoseq_mosta_E16.5_E1S3-%s-cell.rds",as.integer(target_size))
  mylist <- readRDS(file = here('real','input_sc',file))
  
  counts <- mylist$counts
  pos <- mylist$pos
  prop <- mylist$prop
  colnames(pos) <- c("x","y")
  prop <- prop[rowSums(prop)>0,]
  colnames(prop) <- make.names(colnames(prop))
  counts <- counts[,rownames(prop)]
  pos <- pos[rownames(prop),]  
  
  counts <- as.matrix(counts)
  pos <- as.matrix(pos)
  prop <- as.matrix(prop) 
  cell_type <- apply(prop, 1, function(row) {
    return(colnames(prop)[which.max(row)])
  })
  
  ### Create C-SIDE object ---
  gc1 <- gc(reset = TRUE)
  time = system.time({
    puck <- SpatialRNA(coords = as.data.frame(pos), counts = counts)
    reference <- Reference(counts = counts, cell_types = factor(cell_type), min_UMI = -Inf)  
    ct_tab <- table(reference@cell_types)
    keep_ct <- names(ct_tab[ct_tab >= 10])
    cell_type_filter <- cell_type[cell_type %in% keep_ct]
    reference <- Reference(counts = counts[,names(cell_type_filter)], cell_types = factor(cell_type_filter), min_UMI = -Inf)
    
    
    myRCTD <- create.RCTD(spatialRNA = puck, reference = reference, max_cores = ncores,
                          gene_cutoff = -Inf, fc_cutoff = -Inf, gene_cutoff_reg = -Inf, fc_cutoff_reg = -Inf, UMI_min = -Inf,
                          UMI_max = Inf, counts_MIN = -Inf, UMI_min_sigma = -Inf, CELL_MIN_INSTANCE = -Inf)
    myRCTD@config[["MIN_OBS"]] <- -Inf
    myRCTD@config[["MIN_CHANGE_BULK"]] <- -Inf
    myRCTD@config[["MIN_CHANGE_REG"]] <- -Inf
    myRCTD@config[["CONFIDENCE_THRESHOLD"]] <- -Inf
    # myRCTD <- run.RCTD(RCTD = myRCTD)
    myRCTD@config$RCTDmode <- "full"
    myRCTD <- import_weights(myRCTD = myRCTD, weights = prop)
    
    # res_cside <- run.CSIDE.nonparam(myRCTD = myRCTD, cell_type_threshold = -Inf, gene_threshold = -Inf, doublet_mode = FALSE, fdr = Inf)
    
    CSIDE.results <- tryCatch(
      {
        # Core logic: run CSIDE.nonparam while suppressing all warnings
        suppressWarnings(
          run.CSIDE.nonparam(
            myRCTD = myRCTD,
            cell_type_threshold = 0,
            # cell_types = names(tail(sort(colSums(prop)), 3)),
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
    
    saveRDS(res.cside,here('real','res',sprintf('%s-C-SIDE.rds',dataset)))
  })
  
  gc2 <- gc()
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
  
  #fit the ctsvg model
  gc1 <- gc(reset = TRUE)
  time = system.time({
    coord <- pos
    colnames(coord) <- c("row","col")
    cell_types <- apply(prop, 1, function(row) {
      return(colnames(prop)[which.max(row)])
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
    
    
    saveRDS(res.ctsvg,here('real','res',sprintf('%s-ctsvg.rds',dataset)))
  })
  
  gc2 <- gc()
  computation = data.frame(
    time = time[['elapsed']],
    n_spot = ncol(counts),
    n_gene = nrow(counts),
    dataset = dataset,
    Peak_RAM_Used_MiB = sum(gc2[,6] - gc1[,6]),
    method = "ctSVG"
  )
  
  write.table(computation, 
              file =here('real','computation',sprintf('%s-ctsvg.csv',dataset)),
              sep = ',',
              row.names = FALSE)
  
  #fit the spVC model
  gc1 <- gc(reset = TRUE)
  time = system.time({
    boundary <- readRDS(here('real','boundary',file))
    Tr.cell <- TriMesh(boundary, n = 2) # n : triangulation fineness
    V <- as.matrix(Tr.cell$V) 
    Tr <- as.matrix(Tr.cell$Tr)  
    
    # Fit the spVC models safely
    res.spvc <- tryCatch(
      {
        suppressWarnings(
          test.spVC(
            Y = counts,
            X = prop,
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
    
    saveRDS(res.spvc,here('real','res',sprintf('%s-spVC.rds',dataset)))
  })
  
  gc2 <- gc()
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
  
  ### Create Celina object ---
  # celltype_to_test <- names(tail(sort(colSums(prop)),3))
  gc1 <- gc(reset = TRUE)
  time = system.time({
    normalized_counts <- scater::normalizeCounts(counts)
    Obj = Create_Celina_Object(celltype_mat = t(prop), 
                               gene_expression_mat = as.matrix(normalized_counts), 
                               location = as.matrix(pos),
                               covariates = NULL)
    Obj@celltype_mat <- as.matrix(Obj@celltype_mat)
    
    Obj = preprocess_input(Obj, 
                           cell_types_to_test = names(colSums(prop))[colSums(prop)>10] ,  
                           scRNA_count = as.matrix(counts), 
                           sc_cell_type_labels = as.matrix(data.frame(cell_type,row.names = colnames(counts))))
    
    Obj = Calculate_Kernel(Obj)
    Obj = Testing_interaction_all(Obj, num_cores=ncores)
    
    res.celina <- Obj@result
    saveRDS(res.celina,here('real','res',sprintf('%s-CELINA.rds',dataset)))
  })
  
  gc2 <- gc()
  computation = data.frame(
    time = time[['elapsed']],
    n_spot = ncol(counts),
    n_gene = nrow(counts),
    dataset = dataset,
    Peak_RAM_Used_MiB = sum(gc2[,6] - gc1[,6]),
    method = "Celina"
  )
  
  write.table(computation, 
              file =here('real','computation',sprintf('%s-CELINA.csv',dataset)),
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
    
    res.stance = tryCatch(
      {
        suppressWarnings(
          mySTANCE <- runTest2(
            object = mySTANCE, 
            Genes_to_test = utSVG.list, 
            # Cell_types_to_test = names(tail(sort(colSums(prop)), 3)),
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
    saveRDS(res.stance,here('real','res',sprintf('%s-STANCE.rds',dataset)))
  })
  
  gc2 <- gc()
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
  gc1 <- gc(reset = TRUE)
  time = system.time({
    spe <- SpatialExperiment(assay = counts[,rownames(prop)], colData = pos[rownames(prop),], spatialCoordsNames = c('x', 'y')) 
    CTSV.results <- CTSV(spe, W = as.matrix(prop), num_core = ncores) 
    res.ctsv.matrix <- CTSV.results$pval
    res.ctsv <- setNames(
      lapply(seq_along(colnames(prop)), function(i){
        data.frame(pval = pmin(res.ctsv.matrix[,i], res.ctsv.matrix[,i+length(colnames(prop))]),
                   row.names = rownames(CTSV.results$qval))
      }
      ),
      colnames(prop)
    )
    saveRDS(res.ctsv,here('real','res',sprintf('%s-CTSV.rds',dataset)))
  })
  gc2 <- gc()
  computation = data.frame(
    time = time[['elapsed']],
    n_spot = ncol(counts),
    n_gene = nrow(counts),
    dataset = dataset,
    Peak_RAM_Used_MiB = sum(gc2[,6] - gc1[,6]),
    method = "CTSV"
  )
  
  write.table(computation, 
              file =here('real','computation',sprintf('%s-CTSV.csv',dataset)),
              sep = ',',
              row.names = FALSE)
}