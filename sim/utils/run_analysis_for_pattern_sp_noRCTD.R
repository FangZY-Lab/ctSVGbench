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
# library(ctSVG)
library(CTSV)
library(SpatialExperiment)
library(here)
library(pscl)
library(qvalue)

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
  # normalize proportions
  prop.use <- sweep(prop.use, 1, rowSums(prop.use), '/')
  
  # generate spatial data and reference scRNA-seq
  stlist <- generate_spatial_data(pos = pos.use, cell_prop = prop.use, boundary = boundary, pattern = pt, seed = rep_id, dropout_rate = dropout_rate)
  refer.sc <- generate_sc(seed = rep_id)
  
  counts.sc <- refer.sc$expr_mat
  celltypes <- refer.sc$celltypes
  prop <- stlist$prop
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
  
  cell_type <- apply(prop, 1, function(row) {
    return(colnames(prop)[which.max(row)])
  })
  counts <- as.matrix(counts)
  pos <- as.matrix(pos)
  prop <- as.matrix(prop)
  
  # C-SIDE
  puck <- SpatialRNA(coords = as.data.frame(pos), counts = counts)
  reference <- Reference(counts = counts, cell_types = factor(cell_type), min_UMI = -Inf)
  ct_tab <- table(reference@cell_types)
  keep_ct <- names(ct_tab[ct_tab >= 2])
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
  myRCTD <- import_weights(myRCTD = myRCTD, weights = prop)
  
  # res_cside <- run.CSIDE.nonparam(myRCTD = myRCTD, cell_type_threshold = -Inf, gene_threshold = -Inf, doublet_mode = FALSE, fdr = Inf)
  
  CSIDE.results <- tryCatch(
    {
      # Core logic: run CSIDE.nonparam while suppressing all warnings
      suppressWarnings(
        run.CSIDE.nonparam(
          myRCTD = myRCTD,
          cell_type_threshold = 0,
          cell_types = names(tail(sort(colSums(prop)), 3)),
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
  
  saveRDS(res.cside,here('sim','res',sprintf('%s-noRCTD-C-SIDE.rds',dataset)))
  
  
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
                          scRNA_count = as.matrix(counts), 
                          sc_cell_type_labels = as.matrix(data.frame(cell_type, row.names = colnames(counts))))
  
  Obj <- Calculate_Kernel(Obj)
  Obj <- Testing_interaction_all(Obj, num_cores = ncores)
  res.celina <- Obj@result
  saveRDS(res.celina, here('sim','res', sprintf('%s-noRCTD-CELINA.rds', dataset)))
  
  
  Tr.cell <- TriMesh(boundary, n = 2)
  V <- as.matrix(Tr.cell$V) 
  Tr <- as.matrix(Tr.cell$Tr)  
  
  results <- test.spVC(Y = counts, X = prop, S = pos, V = V, Tr = Tr,
                       para.cores = ncores, twostep =F)
  if(is.null(results)){
    res.spVC_1=results
  }else {
    if(all(rowSums(prop==1)>0)){
      ct_total <- colSums(prop)
      celltype_to_test <- names(sort(ct_total, decreasing = TRUE))[1:3]
      idx <- 1:3
      genes.v=names(results$results.full)
    }else{
      celltype_to_test <- names(tail(sort(colSums(prop)),3))
      idx=match(celltype_to_test,colnames(prop))
      genes.v=names(results$results.full)
    }
    if(length(genes.v)==0){
      res.spVC_1=NULL
    }
    res.spVC_1 <- lapply(idx,function(ct){
      
      key <- paste0("gamma_X", ct)
      
      pval <- vapply(genes.v, function(g) {
        x <- results$results.full[[g]]  
        if (!is.list(x)) return(NA_real_)
        pv <- x[["p.value"]]
        if (is.null(pv) || is.null(names(pv)) || !(key %in% names(pv))) return(NA_real_)
        as.numeric(pv[[key]])
      }, numeric(1))
      
      # names(pval)=sapply(strsplit(names(pval),"\\."),"[[",1)
      data.frame(pval = na.omit(pval))
    })    
    names(res.spVC_1) <- celltype_to_test
  }
  saveRDS(res.spVC_1, here('sim','res', sprintf('%s-noRCTD-spVC_1.rds', dataset)))
  
  results <- test.spVC(Y = counts, X = prop, S = pos, V = V, Tr = Tr,
                       para.cores = ncores)
  if(is.null(results)){
    res.spVC_2=results
  }else {
    if(all(rowSums(prop==1)>0)){
      ct_total <- colSums(prop)
      celltype_to_test <- names(sort(ct_total, decreasing = TRUE))[1:3]
      idx <- 1:3
      genes.v=names(results$results.varying)
    }else{
      celltype_to_test <- names(tail(sort(colSums(prop)),3))
      idx=match(celltype_to_test,colnames(prop))
      genes.v=names(results$results.varying)
    }
    if(length(genes.v)==0){
      res.spVC_2=NULL
    }
    res.spVC_2 <- lapply(idx,function(ct){
      pval=sapply(results$results.varying[genes.v],function(x){
        x$p.value[paste0("gamma_X", ct)]
      })
      names(pval)=sapply(strsplit(names(pval),"\\."),"[[",1)
      data.frame(pval = na.omit(pval))
    })
    names(res.spVC_2) <- celltype_to_test
  }
  saveRDS(res.spVC_2, here('sim','res', sprintf('%s-noRCTD-spVC_2.rds', dataset)))
  
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
  saveRDS(res.stance, here('sim','res', sprintf('%s-noRCTD-STANCE.rds', dataset)))
  
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
  saveRDS(res.ctsv[top3_ct],here('sim','res',sprintf('%s-noRCTD-CTSV.rds',dataset)))  
  cat(sprintf("Finished pattern %s", pt))  
}
