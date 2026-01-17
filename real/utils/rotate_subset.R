get_subset_rect <- function(pos,
                            target_rows = 500,
                            aspect = 2,          
                            tol = 50,
                            max_attempts = 10000) {
  if (!all(c("x", "y") %in% colnames(pos))) stop("colnames is not c('x','y')")
  
  cx <- median(pos$x)
  cy <- median(pos$y)
  
  dx <- abs(pos$x - cx)
  dy <- abs(pos$y - cy)
  
  s_low <- 0
  s_high <- max(c(max(dx), max(dy))) * 2  
  
  attempts <- 0
  current_rows <- 0
  subset_idx <- rep(FALSE, nrow(pos))
  
  while (attempts < max_attempts) {
    s_mid <- (s_low + s_high) / 2
    half_w <- s_mid * aspect
    half_h <- s_mid
    
    subset_idx <- (dx <= half_w) & (dy <= half_h)
    current_rows <- sum(subset_idx)
    
    if (abs(current_rows - target_rows) <= tol) break
    
    if (current_rows < target_rows) {
      s_low <- s_mid
    } else {
      s_high <- s_mid
    }
    
    attempts <- attempts + 1
    
    if ((s_high - s_low) < 1e-12) break
  }
  
  pos_subset <- pos[subset_idx, , drop = FALSE]
  message(sprintf("original：%d，now：%d", nrow(pos), nrow(pos_subset)))
  if (attempts >= max_attempts) warning(sprintf("exceed attempt times: %d", max_attempts))
  
  pos_subset
}

datasets <- c(
  "StereoSeq_CBMSTA_Marmoset1_T514",
  "StereoSeq_CBMSTA_Mouse1_T189",
  "StereoSeq_CBMSTA_Mouse2_T349",
  "MERFISH_hypothalamus",
  "SeqFish+_cortex"  
)

for(dataset in datasets){
  file=sprintf('myRCTD_%s.rds',dataset)
  puck<- readRDS(here('real','puck',file))
  pos <- puck@coords
  colnames(pos) <- c("x","y")
  pos.subset <- get_subset_rect(pos)
  saveRDS(pos.subset,here('real','pos_r',file))
}
for (dataset in datasets){
  file=sprintf('myRCTD_%s.rds',dataset)
  pos <- readRDS(here('real','pos_r',file))
  ggplot(pos,aes(x=x,y=y))+
    geom_point()
  ggsave(here('real','pos_r',gsub('.rds','.pdf',file)))
  
}

