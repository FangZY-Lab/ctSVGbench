library(ggplot2)
library(scales)
library(here)
library(tidyr)
library(dplyr)
library(tibble)
library(tidyverse)
library(data.table)

get_wide_pval <- function(dataset){
  prop <- readRDS(here('real','prop',sprintf('myRCTD_%s.rds',dataset)))  
  res.cside=readRDS(here('real','res',sprintf('%s-C-SIDE.rds',dataset)))
  res.celina=readRDS(here('real','res',sprintf('%s-CELINA.rds',dataset)))
  res.stance=readRDS(here('real','res',sprintf('%s-STANCE.rds',dataset)))
  # res.ctsv=readRDS(here('real','res',sprintf('%s-CTSV.rds',dataset)))

  spVC=readRDS(here('real','res',sprintf('%s-spVC.rds',dataset)))
  idx=match(names(res.celina),colnames(prop))
  genes.v=names(spVC$results.varying)
  res.spVC <- lapply(idx,function(ct){
    pval=sapply(spVC$results.varying[genes.v],function(x){
      x$p.value[paste0("gamma_X", ct)]
    })
    names(pval)=sapply(strsplit(names(pval),"\\."),"[[",1)
    data.frame(pval = na.omit(pval))
  })
  
  names(res.spVC) <- names(res.celina)
  
  all_lists <- list(
    CSIDE = res.cside, 
    spVC = res.spVC, 
    Celina = res.celina, 
    STANCE = res.stance
  )
  
  
  all_genes <- unique(unlist(lapply(all_lists, function(lst) {
    unlist(lapply(lst, rownames))
  })))
  
  dat.pval <- do.call(rbind,lapply(names(all_lists),function(method){
    
    lst <- all_lists[[method]]
    
    n <- min(3, length(names(lst)))
    
    list=lapply(names(lst)[1:n],function(i){  
      df <- lst[[i]]
      if(nrow(df)>0){
        
        pval_cols <- grep("p[_\\.]?val|pvalue|pval", colnames(df), ignore.case = TRUE) 
        
        if (length(pval_cols)!=1) {
          stop(sprintf("Data frame %s[[%d]] does not contain column pval ", method, i))
        }
        pvals <- p.adjust(df[,pval_cols],method = "BH")
        names(pvals) <- rownames(df)
        
        pvals_full <- rep(1, length(all_genes))
        names(pvals_full) <- all_genes
        
        pvals_full[names(pvals)] <- pvals
        pvals_full <- data.frame(pval=pvals_full)
        pvals_full$gene=paste0(make.names(i),"_",rownames(pvals_full))
        return(pvals_full)
        
      }else{
        pvals_full <- rep(1, length(all_genes))
        names(pvals_full) <- all_genes
        pvals_full <- data.frame(pval=pvals_full)
        pvals_full$gene=paste0(make.names(i),"_",rownames(pvals_full))
        return(pvals_full)
      }
      
    })
    list=list[!sapply(list,is.null)]
    if(length(list)>0){
      df.pval <- do.call(rbind,list)
      df.pval$method <- method
      
      return(df.pval)
    }else{
      return(NULL)
    }
    
  }))
  
  
  
  dat.pval$method <- recode(dat.pval$method,
                            "CSIDE" = "C-SIDE")
  str(dat.pval)
  dat.pval.wide <- dat.pval %>%
    pivot_wider(
      names_from = method,
      values_from = pval
    ) %>% 
    column_to_rownames( var = "gene") %>% 
    mutate(across(everything(), ~ replace_na(., 1)))
} 
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
