library(ggplot2)
library(scales)
library(here)
library(tidyr)
library(dplyr)
library(tibble)
library(tidyverse)
library(data.table)
source('./real/utils/get_wide_pval.R')
get_inters <- function(dataset){  
  dat.pval.wide <- get_wide_pval(dataset)
  dat.pval.binary <- as.data.frame((dat.pval.wide < 0.05) * 1)
  dat.pval.binary$Method_Count <- rowSums(dat.pval.binary)
  
  dat.plot <- dat.pval.binary %>%
    filter(Method_Count > 0) %>%
    mutate(Gene = rownames(.))
  
  
  dat.plot.long <- dat.plot %>%
    pivot_longer(cols = -c(Gene, Method_Count), names_to = "Method", values_to = "Detected") %>%
    filter(Detected == 1)
  
  inter_freq <- as.data.frame(table(dat.plot$Method_Count))
  colnames(inter_freq) <- c("intersect","freq")
  inter_freq$dataset <- dataset
  return(inter_freq)
}

get_conc<- function(dat,dataset=dataset) {
  nfeatures <- nrow(dat)
  k_max <- 200
  
  max_iters <- ncol(dat) * ncol(dat) * k_max
  df_res <- data.table(
    method1 = character(max_iters),
    method2 = character(max_iters),
    rank = integer(max_iters),
    conc = numeric(max_iters),
    dataset = character(max_iters),
    nfeatures = numeric(max_iters)
  )
  
  counter <- 1
  
  for (i in 1:ncol(dat)) {
    for (j in 1:ncol(dat)) {
      if (i != j) {      
        
        i.sorted <- names(sort(setNames(dat[, i], rownames(dat)), decreasing = FALSE))
        j.sorted <- names(sort(setNames(dat[, j], rownames(dat)), decreasing = FALSE))
        nfeatures <- nrow(dat)
        
        for (k in 1:k_max) {
          
          i.top <- i.sorted[1:k]
          j.top <- j.sorted[1:k]
          
          inters <- length(intersect(i.top, j.top))
          conc <- inters / k
          
          tmp <- data.table(
            method1 = colnames(dat)[i],
            method2 = colnames(dat)[j],
            rank = k,
            conc = conc,
            dataset = dataset,
            nfeatures = nfeatures
          )
          df_res[counter, ] <- tmp
          
          
          counter <- counter + 1
        }
      }
    }
  }
  df_res <- df_res[1:(counter - 1), ]
  
}
