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
CTSV <- function(spe, W, num_core=1, BPPARAM = NULL){
  if (missing(spe) || !is(spe,"SpatialExperiment") || is.null(rownames(spe)) || is.null(colnames(spe))) {
    stop("Include SpatialExperient class object with rownames and colnames")
  }
  if (missing(W) || !is.matrix(W)) {
    stop("Include cell-type proportion matrix of the matrix type.")
  } 
  if(as.integer(num_core) != as.numeric(num_core)){
    stop("Input integer num of cores.")
  }
  Y <- as.matrix(t(assay(spe)))
  
  # Y <- t(assay(spe))
  loc <- spatialCoords(spe)
  if(sum(is.na(Y))>0 | sum(is.na(loc))>0 || sum(is.na(W)) > 0 || sum(rowSums(Y) == 0)>0 || sum(colSums(Y) == 0)>0 || sum(colSums(W) == 0)>0 || sum(rowSums(W) == 0)>0){
    stop("Remove NaNs, columns with all zeros and rows with all zeros in datasets.")
  }
  if(nrow(loc)!= nrow(W) || sum(rownames(W) != colnames(spe))>0){
    stop("Keep the number and names of spots consistent in gene expression matrix, location coordinate matrix and cell-type proportion matrix.")
  }
  if (is.null(BPPARAM)) {
    BPPARAM <- BiocParallel::MulticoreParam(workers = num_core)
  }
  # make sure the sum of cell type proportions is equal to 1 in each spot.    
  W <- W / rowSums(W)
  # number of genes
  G <- ncol(Y)
  # number of spots
  n <- nrow(loc)
  # number of cell types
  K <- ncol(W)
  # normalize cell-type proportion matrix W to ensure the summation across cell types in one spot is equal to one.
  W <- W / rowSums(W)
  # Center and normalize coordinates of spots to have mean zero and standard deviation one.
  S <- t(loc) - colMeans(loc)
  S <- t(S / apply(S, 1, sd))
  quan <- c(0.4,0.6)
  psi1 <- quantile(abs(S[,1]), quan)
  psi2 <- quantile(abs(S[,2]), quan)
  P_VAL <- array(NA, dim = c(G, 2*K, 5))
  pattern <- c("linear","gau1","gau2","cos1","cos2")
  for(fit_pat in pattern){
    if(fit_pat == "gau1"){
      h1 <- exp(-S[,1]^2 / 2 / psi1[1]^2)
      h2 <- exp(-S[,2]^2 / 2 / psi2[1]^2)
    }else if(fit_pat == "gau2"){
      h1 <- exp(-S[,1]^2 / 2 / psi1[2]^2)
      h2 <- exp(-S[,2]^2 / 2 / psi2[2]^2)
    }else if(fit_pat == "cos1"){
      h1 <- cos(2*pi*S[,1] / psi1[1])
      h2 <- cos(2*pi*S[,2] / psi2[1])
    } else if(fit_pat == "cos2"){
      h1 <- cos(2*pi*S[,1] / psi1[2])
      h2 <- cos(2*pi*S[,2] / psi2[2])
    }else{
      h1 <- S[,1]
      h2 <- S[,2]
    }
    # print(fit_pat)
    Tmp <- cbind(W * h1, W * h2, W)
    colnames(Tmp) <- seq_len(ncol(Tmp))
    res <- do.call(rbind,BiocParallel::bplapply(seq_len(G), CTSV:::.P_gene,BPPARAM = BPPARAM,Y=Y,Tmp = Tmp,h1=h1,h2=h2))
    P_VAL[,,match(fit_pat,pattern)] <- res
    rownames(P_VAL[,,match(fit_pat,pattern)]) <- colnames(Y)
  }
  P_VAL[which(is.na(P_VAL))] <- 1
  P_VAL[P_VAL == -1] <- 1
  P_VAL <- tan((0.5 - P_VAL)*pi)
  T_cau0 <- apply(P_VAL, c(1,2), mean)
  P_val <- 1-pcauchy(T_cau0)
  # convert q-values into q-values
  Q_val <- matrix(qvalue(c(P_val))$qvalue, G, 2*K)
  rownames(Q_val) <- colnames(Y)
  return(list("pval" = P_val,
              "qval" = Q_val))
}

run_analysis_for_pattern <- function(pt,  
                                     pos.use, prop.use, dt = dataset, boundary, 
                                     rep_id = 1, paramset='P1',
                                     ncores = 10) {
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
  
  
  # spVC
  Tr.cell <- TriMesh(boundary, n = 2)
  V <- as.matrix(Tr.cell$V) 
  Tr <- as.matrix(Tr.cell$Tr)  
  results <- test.spVC(Y = counts, X = prop, S = pos, V = V, Tr = Tr,
                       para.cores = ncores)
  
  saveRDS(results, here('sim','res', sprintf('%s-noRCTD-spVC.rds', dataset)))
  
  
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
    