library(tidyverse)
library(readxl)
library(hrbrthemes)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(ggtext)
library(ggprism)
library(ggpmisc)
library(hrbrthemes)
library(pals)
library(viridis)
library(viridisLite)
library(scico)
library(reshape2)
library(ggrepel)
library(tibble)
library(pROC)
library(aplot)
library(here)
library(ggpubr)
library(dplyr)
library(tidyr)
library(dplyr)
library(tidyr)
library(ggpubr)
source('./my_theme.R')
source('./sim/utils/sim-bench.R')


dataset=  "SeqFish+_cortex"

methods <- c('C-SIDE','spVC','CELINA','STANCE',"CTSV","ctsvg")


reps <- c(1:100)
get_wide_pval <- function(dataset,i){
  prop <- readRDS(here('real','prop',sprintf('myRCTD_%s.rds',dataset)))  
  res.cside=readRDS(here('real','res',sprintf('%s-C-SIDE-null%s.rds',dataset,i)))
  res.celina=readRDS(here('real','res',sprintf('%s-CELINA-null%s.rds',dataset,i)))
  res.stance=readRDS(here('real','res',sprintf('%s-STANCE-null%s.rds',dataset,i)))
  # res.ctsv=readRDS(here('real','res',sprintf('%s-CTSV.rds',dataset)))
  res.ctsv=readRDS(here('real','res',sprintf('%s-CTSV-null%s.rds',dataset,i)))
  
  spVC=readRDS(here('real','res',sprintf('%s-spVC-null%s.rds',dataset,i)))
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
  
  ctsvg=readRDS(here('real','res',sprintf('%s-ctsvg-null%s.rds',dataset,i)))
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

datasets.all <- unlist(lapply(reps,function(i){
  unlist(lapply(methods,function(m){
    sprintf('%s-%s-null%s.rds',dataset,m,i)
  }))
}))

FPR_0.01 <- sapply(reps,function(i){
  pval=get_wide_pval(dataset = dataset,i)
  colSums(pval<0.01)/dim(pval)[1]
}) %>% 
  as.data.frame() %>%  
  tibble::rownames_to_column(var = "Method") %>% 
  gather(key = "sim", value = "FPR_0.01", -Method)

FPR_0.05 <- sapply(reps,function(i){
  pval=get_wide_pval(dataset = dataset,i)
  colSums(pval<0.05)/dim(pval)[1]
}) %>% 
  as.data.frame() %>%  
  tibble::rownames_to_column(var = "Method") %>% 
  gather(key = "sim", value = "FPR_0.05", -Method)

FPR_0.1 <- sapply(reps,function(i){
  pval=get_wide_pval(dataset = dataset,i)
  colSums(pval<0.1)/dim(pval)[1]
}) %>% 
  as.data.frame() %>%  
  tibble::rownames_to_column(var = "Method") %>% 
  gather(key = "sim", value = "FPR_0.1", -Method)


FPR <- cbind(FPR_0.01, FPR_0.05 = FPR_0.05$FPR, FPR_0.1 = FPR_0.1$FPR) %>%
  as.data.frame() %>%  
  gather(key = "cutoff", value = "FPR", -Method, -sim)

median.FPR <- plyr::ddply(FPR, ~cutoff + Method, function(x) {
  FPR <- median(x$FPR) 
  c("FPR" = FPR)
})

rank.FPR <- plyr::ddply(median.FPR, ~cutoff, function(x) {
  rank_value <- rank(-(x$FPR), ties.method = "average")  
  data.frame(x, rank = rank_value)
})

rank.FPR.mean <- plyr::ddply(rank.FPR, ~Method, function(x) {
  rank <- mean(x$rank)  
  c("rank" = rank)
})

ord0.01 <- FPR_0.01 %>% 
  group_by(Method) %>% 
  summarise(median = median(FPR_0.01)) %>% 
  mutate(order = order(median))

ord0.05 <- FPR_0.05 %>% 
  group_by(Method) %>% 
  summarise(median = median(FPR_0.05)) %>% 
  mutate(order = order(median))

ord0.1 <- FPR_0.1 %>% 
  group_by(Method) %>% 
  summarise(median = median(FPR_0.1)) %>% 
  mutate(order = order(median))

FPR_0.01$Method <- factor(FPR_0.01$Method, levels = unique(ord0.01$Method)[ord0.01$order], ordered = TRUE)
FPR_0.05$Method <- factor(FPR_0.05$Method, levels = unique(ord0.05$Method)[ord0.05$order], ordered = TRUE)
FPR_0.1$Method <- factor(FPR_0.1$Method, levels = unique(ord0.1$Method)[ord0.1$order], ordered = TRUE)

p5 <- ggplot(data = FPR, aes(color = Method)) +
  geom_boxplot(data = FPR_0.01, aes(x = "0.01", y = 1-FPR_0.01)) +
  geom_boxplot(data = FPR_0.05, aes(x = "0.05", y = 1-FPR_0.05)) +
  geom_boxplot(data = FPR_0.1, aes(x = "0.1", y = 1-FPR_0.1)) +
  # coord_trans(y = "sqrt") +
  geom_segment(aes(x = 1 - 0.5, xend = 1 + 0.5, y = 1-0.01, yend = 1-0.01), color = "red", linewidth = 0.2, lty = 2) +
  geom_segment(aes(x = 2 - 0.5, xend = 2 + 0.5, y = 1-0.05, yend = 1-0.05), color = "red", linewidth = 0.2, lty = 2) +
  geom_segment(aes(x = 3 - 0.5, xend = 3 + 0.5, y = 1-0.1, yend = 1-0.1), color = "red", linewidth = 0.2, lty = 2) +
  xlab("Nominal False Discovery Rate") + ylab("Specificity") +
  scale_color_manual(values = method_colors)+
  theme_minimal()+
  my_theme +
  theme(legend.position = "left",
  # legend.direction = "horizontal",
  #   legend.box = "horizontal",     
    legend.margin = margin(0, 0, 0, 0),
            legend.key.size = unit(0.05, "in"),
        legend.key.width = unit(0.05, "in"))+
  guides(
    color = guide_legend(
      ncol = 1,                
      byrow = TRUE,            
      title.hjust = 0.5    )
      
  )+labs(color="")
p5


summary_df_acc3 <- FPR %>% 
  mutate(cutoff=gsub('FPR_','FDR_',cutoff)) %>% 
  group_by(Method,cutoff) %>%
  summarise(
    Specificity = mean(1-FPR, na.rm = TRUE)  ) %>% 
  pivot_wider(
    names_from = cutoff,      
    values_from = Specificity,         
    names_glue = "{cutoff}_specificity"  
  ) %>% column_to_rownames('Method')

summary_df <- read.csv( "metrics_summary.csv", row.names = 1)
summary_df$specificity <- summary_df_acc3[rownames(summary_df),"FDR_0.05_specificity"]

write.csv(summary_df, "metrics_summary.csv", row.names = T)  
