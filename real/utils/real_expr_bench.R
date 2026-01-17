##real_expr_bench
get_pval_Stroma <- function(dataset){
  prop <- readRDS(here('real','prop',sprintf('myRCTD_%s.rds',dataset)))  
  res.cside=readRDS(here('real','res',sprintf('%s-C-SIDE.rds',dataset)))
  res.celina=readRDS(here('real','res',sprintf('%s-CELINA.rds',dataset)))
  res.stance=readRDS(here('real','res',sprintf('%s-STANCE.rds',dataset)))

ctsvg=readRDS(here('real','res',sprintf('%s-ctsvg.rds',dataset)))
  res.ctsvg <- split(ctsvg, ctsvg$cluster)
  res.ctsvg <- lapply(res.ctsvg, \(df)
                      data.frame(
                        pval = df$pval,
                        row.names = df$gene
                      ))  
  res.ctsv=readRDS(here('real','res',sprintf('%s-CTSV.rds',dataset)))                       
  
spVC=readRDS(here('real','res',sprintf('%s-spVC.rds',dataset)))
  if(is.null(spVC)){
    res.spVC=NULL
  }else {
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
  }
  
  all_lists <- list(CSIDE = res.cside,
                    spVC = if (is.null(spVC)) NULL else res.spVC,
                    Celina = res.celina, 
                    STANCE = res.stance,
                    ctSVG = res.ctsvg,
                    CTSV = res.ctsv)
  
  all_genes <- unique(unlist(lapply(all_lists, function(lst) {
    unlist(lapply(lst, rownames))
  })))
  
  dat.pval <- do.call(rbind,lapply(names(all_lists),function(method){
    
    lst <- all_lists[[method]]
    
    n <- min(1, length(names(lst)))
    
    list=lapply(c("Stroma"),function(i){  
      df <- lst[[i]]
      if(!is.null(df)){
    
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
        pvals_full$gene=rownames(pvals_full)
        return(pvals_full)
        
      }else{
        pvals_full <- rep(1, length(all_genes))
        names(pvals_full) <- all_genes
        pvals_full <- data.frame(pval=pvals_full)
        pvals_full$gene=rownames(pvals_full)
        return(pvals_full)
      }
      
      }else{
        pvals_full <- rep(1, length(all_genes))
        names(pvals_full) <- all_genes
        pvals_full <- data.frame(pval=pvals_full)
        pvals_full$gene=rownames(pvals_full)
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
  
  str(dat.pval)
    dat.pval$method=ifelse(dat.pval$method=="CSIDE","C-SIDE",dat.pval$method)

  
  dat.pval.wide <- dat.pval %>%
    pivot_wider(
      names_from = method,
      values_from = pval
    ) %>% 
    column_to_rownames( var = "gene") %>% 
    mutate(across(everything(), ~ replace_na(., 1)))
} 
