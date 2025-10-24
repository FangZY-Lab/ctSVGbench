generate_spatial_data <- function(
    G = 1200, 
    N = nrow(pos), 
    mu = 2, 
    alpha = 1, 
    pos, 
    cell_prop,
    boundary,
    pattern = c("hotspot", "stripe", "gradient", "neighbor", "pathology", "periodic"),
    shared_ct = c(5,6),  # shared cell types
    seed = 123
) {
  set.seed(seed)
  pattern <- match.arg(pattern)
  svg_num <- 200
  ct_num <- ncol(cell_prop)
  marker_num <- 25
  
  expr_mat <- matrix(NA, nrow = G, ncol = N)
  up_idx <- sample(1:svg_num, size = 10)
  down_idx <- setdiff(1:svg_num, up_idx)
  region_label <- rep(0, N)
  
  # true_label: 1=unique to shared_ct[1], 2=unique to shared_ct[2], 3=shared
  true_label <- numeric(G)
  true_label[1:75] <- 1
  true_label[76:150] <- 2
  true_label[151:200] <- 3
  
  # spatial patterns
  if (pattern == "hotspot") {
    center <- c(quantile(pos[,1], 0.55), quantile(pos[,2], 0.45))
    radius <- min(max(pos[,1]) - min(pos[,1]), max(pos[,2]) - min(pos[,2])) * 0.2
    region_label <- as.numeric((pos[,1] - center[1])^2 + (pos[,2] - center[2])^2 <= radius^2)
    
  } else if (pattern == "stripe") {
    x_mid <- quantile(pos[,1], 0.55)
    width <- (max(pos[,1]) - min(pos[,1])) / 10
    region_label <- as.numeric(abs(pos[,1] - x_mid) <= width)
    
  } else if (pattern == "gradient") {
    region_label <- (pos[,1] - min(pos[,1])) / (max(pos[,1]) - min(pos[,1]))
    region_label <- pmax(pmin(region_label, 1), 0.1)  
    
  } else if (pattern == "neighbor") {
    density <- cell_prop[, shared_ct[2]]
    region_label <- density / max(density)
    
  } else if (pattern == "pathology") {
    centers <- matrix(runif(4 * 2, min = min(pos), max = max(pos)), ncol = 2)
    region_label <- apply(pos, 1, function(p) {
      min(apply(centers, 1, function(c) sum((p - c)^2)))
    })
    region_label <- ifelse(region_label < quantile(region_label, 0.2), 1, 0)
    
  } else if (pattern == "periodic") {
    freq <- 6 * pi / (max(pos[,1]) - min(pos[,1]))  
    region_label <- sin(freq * pos[,1])  
    region_label <- (region_label - min(region_label)) / (max(region_label) - min(region_label)) 
    region_label <- pmax(region_label, 0.1)  
  }
  
  # fold-change matrix
  fold_change_matrix <- matrix(NA, nrow = svg_num, ncol = length(shared_ct))
  colnames(fold_change_matrix) <- paste0("celltype", shared_ct)
  rownames(fold_change_matrix) <- paste0("gene", 1:svg_num)
  
  # SVG gene generation
  for (i in 1:svg_num) {
    if (true_label[i] == 1) {
      base_mean.main <- mu * cell_prop[, shared_ct[1]]
      other_mean <- mu * rowSums(cell_prop[, -shared_ct[1]])
      fold_change <- ifelse(i %in% up_idx, sample(c(5, 10), 1),
                            sample(c(0.25, 0.5), 1))
      if (pattern %in% c("neighbor","gradient","periodic")) { 
        adj_mean.main <- base_mean.main * region_label * fold_change 
      } else { 
        adj_mean.main <- base_mean.main * ifelse(region_label == 1, fold_change, 1)
      }
      total_mean <- other_mean + adj_mean.main
      
    } else if (true_label[i] == 2) {
      base_mean.main <- mu * cell_prop[, shared_ct[2]]
      other_mean <- mu * rowSums(cell_prop[, -shared_ct[2]])
      fold_change <- ifelse(i %in% up_idx, sample(c(5, 10), 1),
                            sample(c(0.25, 0.5), 1))
      if (pattern %in% c("neighbor","gradient","periodic")) { 
        adj_mean.main <- base_mean.main * region_label * fold_change 
      } else { 
        adj_mean.main <- base_mean.main * ifelse(region_label == 1, fold_change, 1)
      }
      total_mean <- other_mean + adj_mean.main
      
    } else if (true_label[i] == 3) {
      # shared genes: independent fold-change for two cell types
      fc1 <- ifelse(i %in% up_idx, sample(c(5, 10), 1),
                    sample(c(0.25, 0.5), 1))
      fc2 <- ifelse(i %in% up_idx, sample(c(5, 10), 1),
                    sample(c(0.25, 0.5), 1))
      fold_change_matrix[i, 1] <- fc1
      fold_change_matrix[i, 2] <- fc2
      
      adj1 <- mu * cell_prop[, shared_ct[1]]
      adj2 <- mu * cell_prop[, shared_ct[2]]
      if (pattern %in% c("neighbor","gradient","periodic")) {
        adj1 <- adj1 * region_label * fc1
        adj2 <- adj2 * region_label * fc2
      } else {
        adj1 <- adj1 * ifelse(region_label == 1, fc1, 1)
        adj2 <- adj2 * ifelse(region_label == 1, fc2, 1)
      }
      base_mean.main <- adj1 + adj2
      other_mean <- mu * rowSums(cell_prop[, -shared_ct])
      total_mean <- other_mean + base_mean.main
    }
    
    expr_mat[i, ] <- rnbinom(N, mu = floor(total_mean), size = 1 / alpha)
  }
  
  # cell-type marker genes
  for (idx in (svg_num + 1):(svg_num + ct_num * marker_num)) {
    ct_rev <- ct_num - ((idx - svg_num - 1) %/% marker_num)
    total_mean <- mu * rowSums(cell_prop[, -ct_rev, drop = FALSE]) + mu * cell_prop[, ct_rev] * 10
    expr_mat[idx, ] <- rnbinom(N, mu = total_mean, size = 1 / alpha)
  }
  
  # background genes
  for (i in (svg_num + ct_num * marker_num + 1):G) {
    expr_mat[i, ] <- rnbinom(N, mu = mu, size = 1 / alpha)
  }
  
  # dropout
  dropout_rate <- 0.3
  zero_mask <- matrix(runif(length(expr_mat)) < dropout_rate, 
                      nrow = nrow(expr_mat), ncol = ncol(expr_mat))
  expr_mat[zero_mask] <- 0
  
  # output
  rownames(expr_mat) <- paste0("gene", 1:G)
  colnames(expr_mat) <- paste0("spot", 1:N)
  colnames(cell_prop) <- paste0("celltype", 1:ct_num)
  rownames(cell_prop) <- paste0("spot", 1:N)
  rownames(pos) <- paste0("spot", 1:N)
  
  return(list(
    prop = cell_prop, 
    counts = expr_mat, 
    pos = pos,
    boundary = boundary, 
    label = true_label,
    fold_change_matrix = fold_change_matrix
  ))
}
