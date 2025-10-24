generate_sc <- function(
    G = 1200, 
    mu = 1, 
    alpha = 1, 
    seed = 123 
){ 
  set.seed(seed)
  ct_num <- 6  # number of cell types
  cells_per_type <- 200  # cells per type
  N <- cells_per_type * ct_num
  
  # initialize expression matrix
  expr_mat <- matrix(0, nrow = G, ncol = N)
  rownames(expr_mat) <- paste0("gene", 1:G)
  colnames(expr_mat) <- paste0("cell", 1:N)

  for (ct in 1:ct_num) {
    # column index range for current cell type
    start_idx <- (ct - 1) * cells_per_type + 1
    end_idx <- ct * cells_per_type

    # marker gene range
    marker_genes <- ((ct - 1) * 25 + 201):((ct - 1) * 25 + 225)

    # generate counts
    for (i in start_idx:end_idx) {
      # non-marker genes
      expr_mat[-marker_genes, i] <- rnbinom(G - 25, mu = mu, size = 1 / alpha)
      # marker genes
      expr_mat[marker_genes, i] <- rnbinom(25, mu = mu * 10, size = 1 / alpha)
    }
  }

  # apply dropout
  dropout_mask <- matrix(rbinom(G * N, 1, 0.2), nrow = G, ncol = N)
  expr_mat[dropout_mask == 1] <- 0

  rownames(expr_mat) <- paste0("gene", 1:G)
  colnames(expr_mat) <- paste0("cell", 1:N)
  celltypes <- rep(paste0('celltype', 6:1), each = 200)
  names(celltypes) <- paste0("cell", 1:N)

  return(list(expr_mat = expr_mat, celltypes = celltypes))
}
