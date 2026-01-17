library(ggplot2)
library(ggpubr)
library(reshape2)
library(dplyr)
library(ggrepel)
library(here)
library(tidyr)
library(dplyr)
library(tibble)
library(tidyverse)

svg_id <- c(paste0("celltype5","gene",1:75),paste0("celltype6","gene",76:150),
            paste0("celltype5","gene",151:200),paste0("celltype6","gene",151:200))

### get_pvalue_wide ---
get_pvalue_wide <- function(dataset,svg_id,cell.level=F){
  
  res.cside=readRDS(here('sim','res',sprintf('%s-C-SIDE.rds',dataset)))
  res.celina=readRDS(here('sim','res',sprintf('%s-CELINA.rds',dataset)))
  res.stance=readRDS(here('sim','res',sprintf('%s-STANCE.rds',dataset)))
  res.ctsv <- readRDS(here("sim","res", sprintf("%s-CTSV.rds", dataset)))
  # prop <- readRDS(here('sim','prop',sprintf('myRCTD_%s.rds',dataset)))  
  
  spVC=readRDS(here('sim','res',sprintf('%s-spVC.rds',dataset)))
  if(is.null(spVC)){
    res.spVC=spVC
  }else {
    idx=match(names(res.celina),paste0("celltype",1:6))
    genes.v=names(spVC$results.varying)
    res.spVC <- lapply(idx,function(ct){
      pval=sapply(spVC$results.varying[genes.v],function(x){
        x$p.value[paste0("gamma_X", ct)]
      })
      names(pval)=sapply(strsplit(names(pval),"\\."),"[[",1)
      data.frame(pval = na.omit(pval))
    })    
    names(res.spVC) <- names(res.celina)

    # if(sum(sapply(res.spVC,\(i) nrow(i)))==0){
    #   res.spVC=NULL}
    
  }
  
  if(cell.level==TRUE){
    ctsvg=readRDS(here('sim','res',sprintf('%s-ctsvg.rds',dataset)))
    if(is.null(ctsvg)){
      res.ctsvg=ctsvg
    }else {
      res.ctsvg <- split(ctsvg, ctsvg$cluster)
      res.ctsvg <- lapply(res.ctsvg, \(df)
                          data.frame(
                            pval = df$pval,
                            row.names = df$gene
                          ))
    }
    
    all_lists <- list(
      CSIDE = res.cside, 
      spVC = res.spVC, 
      Celina = res.celina, 
      STANCE = res.stance,
      CTSV = res.ctsv, 
      ctSVG=res.ctsvg 
    )
  }else {
    all_lists <- list(
      CSIDE = res.cside, 
      spVC = res.spVC, 
      Celina = res.celina, 
      STANCE = res.stance,
      CTSV = res.ctsv  
    )
  }  
  
  all_genes <- unique(unlist(lapply(all_lists, function(lst) {
    unlist(lapply(lst, rownames))
  })))
  method <- c("CSIDE","spVC","Celina","STANCE","CTSV","ctSVG" )[5]
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
        pvals_full$gene=paste0(make.names(i),rownames(pvals_full))
        return(pvals_full)
        
      }else{
        pvals_full <- rep(1, length(all_genes))
        names(pvals_full) <- all_genes
        pvals_full <- data.frame(pval=pvals_full)
        pvals_full$gene=paste0(make.names(i),rownames(pvals_full))
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
  
  
  head(dat.pval.wide)
  label <- numeric(nrow(dat.pval.wide))
  names(label) <- rownames(dat.pval.wide)
  label[svg_id] <- 1
  
  return(list(dat.pval.wide=dat.pval.wide,label=label))
}


