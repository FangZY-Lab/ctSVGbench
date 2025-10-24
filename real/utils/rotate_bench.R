
get_wide_pval <- function(dataset,method){
  if(method  %in% c('C-SIDE','CELINA','STANCE')){
    
    res.o <- readRDS(here('real','res',sprintf('%s-%s.rds',dataset,method)))
    res.r30 <- readRDS(here('real','res',sprintf('%s-r30-%s.rds',dataset,method)))
    # res.r60 <- readRDS(here('real','res',sprintf('%s-r60-%s.rds',dataset,method)))
    res.r90 <- readRDS(here('real','res',sprintf('%s-r90-%s.rds',dataset,method)))
    
  }else if(method=='spVC'){
    
    res.o <- pre.res(dataset,method,angle=0)
    res.r30 <- pre.res(dataset,method,angle=30)
    # res.r60 <- pre.res(dataset,method,angle=60)
    res.r90 <- pre.res(dataset,method,angle=90)
  }

  
  all_lists <- list(
    original = res.o, 
    rotate.30 = res.r30, 
    # rotate.60 = res.r60, 
    rotate.90 = res.r90
  )
  
  
  all_genes <- unique(unlist(lapply(all_lists, function(lst) {
    unlist(lapply(lst, rownames))
  })))
  dat.pval <- do.call(rbind,lapply(names(all_lists),function(angle){
    
    lst <- all_lists[[angle]]
    
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
        pvals_full$gene=paste0(make.names(i),'_',rownames(pvals_full))
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
      df.pval$angle <- angle
      return(df.pval)
    }else{
      return(NULL)
    }
    
  }))
  
  str(dat.pval)
  
  
  dat.pval.wide <- dat.pval %>%
    pivot_wider(
      names_from = angle,
      values_from = pval
    ) %>% 
    column_to_rownames( var = "gene") %>% 
    mutate(across(everything(), ~ replace_na(., 1)))
  
}

get_conc<- function(dat,dataset=dataset,method=method) {
  nfeatures <- nrow(dat)
  k_max <- 500
  
  max_iters <- ncol(dat) * ncol(dat) * k_max
  df_res <- data.table(
    angle1 = character(max_iters),
    angle2 = character(max_iters),
    method = character(max_iters),
    rank = integer(max_iters),
    conc = numeric(max_iters),
    dataset = character(max_iters),
    nfeatures = numeric(max_iters)
  )
  
  counter <- 1
  
  for (i in 1) {
    for (j in 2:ncol(dat)) {
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
            angle1 = colnames(dat)[i],
            angle2 = colnames(dat)[j],
            method = method,
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

pre.res <- function(dataset,method,angle=0){
  prop <- readRDS(here('real','prop',sprintf('myRCTD_%s.rds',dataset))) 
  res.celina <- readRDS(here('real','res',sprintf('%s-CELINA.rds',dataset)))
  if(angle==0){
    file=here('real','res',sprintf('%s-spVC.rds',dataset))
  }else if(angle==30){
    file=here('real','res',sprintf('%s-r30-%s.rds',dataset,method))
  }else if(angle==60){
    file=here('real','res',sprintf('%s-r60-%s.rds',dataset,method))
  }else if(angle==90){
    file=here('real','res',sprintf('%s-r90-%s.rds',dataset,method))
  }
  
  spVC=readRDS(file)
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
  return(res.spVC)
}