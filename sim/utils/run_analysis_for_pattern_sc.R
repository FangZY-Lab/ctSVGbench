library(RhpcBLASctl)
blas_set_num_threads(1)

library(spacexr) # C-SIDE
library(spVC)
library(sp)
library(BPST)
library(Triangulation)
library(MGLM)
library(CELINA)
library(STANCE)
library(SpatialExperiment)
library(Matrix)
library(devtools)
library(reshape2)
library(dplyr)
library(here)
library(ctSVG)
library(CTSV)
library(pscl)
library(qvalue)
library(Seurat)
library(arrow)
library(jsonlite)

run_analysis_for_pattern <- function(pt,  
                                     pos.use, prop.use, dt = dataset, boundary, 
                                     rep_id = 1, paramset='P1',
                                     ncores = 10) {
  source(here('sim','utils','generate_sc.R'))
  source(here('sim','utils','generate_st.R'))
  dropout_rate <- case_when(
    paramset == 'P1' ~ 0.1,
    paramset == 'P2' ~ 0.2,
    paramset == 'P3' ~ 0.3,
    TRUE ~ 0
  )
  prop.use <- sweep(prop.use, 1, rowSums(prop.use), '/')
  
  stlist <- generate_spatial_data(mu = 1, pos = pos.use, cell_prop = prop.use, boundary = boundary, pattern = pt, seed = rep_id, cell.level = TRUE, dropout_rate = dropout_rate)
  
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
  
  # dataset name with replicate ID
  dataset <- sprintf("sim_%s-%s-%s-rep%d", dt, pt, paramset, rep_id)
  prop <- stlist$prop
  counts <- counts[,rownames(prop)]
  pos <- pos[rownames(prop),]
  
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
  
  
  ## Create C-SIDE object ---
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
  
  saveRDS(res.cside,here('sim','res',sprintf('%s-C-SIDE.rds',dataset)))
  
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
  saveRDS(res.celina,here('sim','res',sprintf('%s-CELINA.rds',dataset)))
  
  ### Create STANCE object ---  
  ncores=120
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
  saveRDS(res.stance,here('sim','res',sprintf('%s-STANCE.rds',dataset)))
  
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
  saveRDS(res.ctsvg,here('sim','res',sprintf('%s-ctsvg.rds',dataset)))
  
  
  #fit the spVC model
  Tr.cell <- TriMesh(boundary, n = 2) # n : triangulation fineness
  V <- as.matrix(Tr.cell$V) 
  Tr <- as.matrix(Tr.cell$Tr)  
  
  # Fit the spVC models safely
  res.spvc <- tryCatch(
    {
      suppressWarnings(
        test.spVC(
          Y = counts,
          X = as.matrix(prop_top3),
          S = as.matrix(pos),
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
  
  saveRDS(res.spvc,here('sim','res',sprintf('%s-spVC_2.rds',dataset)))
  
  res.spvc_1step <- tryCatch(
    {
      suppressWarnings(
        test.spVC(
          Y = counts,
          X = as.matrix(prop_top3),
          S = as.matrix(pos),
          V = V,
          Tr = Tr,
          para.cores = ncores,
          twostep =F
        )
      )
      
    },
    error = function(e) {
      message("test.spvc_1step failed: ", e$message)
      NULL
    }
  )
  
  saveRDS(res.spvc_1step,here('sim','res',sprintf('%s-spVC_1.rds',dataset)))
  
  
  #fit the CTSV model 
  spe <- SpatialExperiment(assay = as.matrix(counts[,rownames(prop_top3)]), colData = pos[rownames(prop_top3),], spatialCoordsNames = c('x', 'y')) 
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
  saveRDS(res.ctsv,here('sim','res',sprintf('%s-CTSV.rds',dataset)))  
  
  
  cat(sprintf("Finished pattern %s", pt))
}
