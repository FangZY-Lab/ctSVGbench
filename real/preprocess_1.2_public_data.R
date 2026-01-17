library(spacexr) #C-SIDE
library(SpatialExperiment)
library(tidyverse)
library(here)
library(data.table)
library(Matrix)
library(Seurat)
library(SummarizedExperiment)

## Visium_lymphnode ---
counts.sc <- read.csv("~/yh/data/lymphnode_sc/counts.sc.csv")
celltypes <- read.csv("~/yh/data/lymphnode_sc/celltype.csv",row.names=1)
celltypes <- setNames(celltypes$new_celltype,rownames(celltypes))
celltypes <- factor(celltypes)
rownames(counts.sc) <- counts.sc$X
counts.sc$X <- NULL
counts.sc <- t(counts.sc)

reference <- Reference(counts = counts.sc, cell_types = celltypes, min_UMI = 1) 
saveRDS(reference,here('ctSVGbench','real','reference',"myRCTD_Visium_lymph_node.rds"))

## Visium_liver ---
obj=readRDS('~/yh/liver_mouseStSt_9celltypes.rds')
sc_idx <- obj$typeSample == "scRnaSeq"

counts.all <- obj@assays$RNA@counts
counts.sc <- counts.all[, sc_idx]
celltypes <- factor(obj$annot_cd45[sc_idx])
dim(counts.sc)
reference <- Reference(counts = counts.sc, cell_types = celltypes, min_UMI = 1) 
saveRDS(reference,here('ctSVGbench','real','reference',"myRCTD_Visium_liver.rds"))

obj=readRDS('~/yh/liver_mouseVisium_JB04.rds')
counts <- obj@assays$Spatial@counts
pos <- obj@images$image@coordinates[,c("imagerow","imagecol")]
colnames(pos) <- c('x','y')
puck <- SpatialRNA(coords=pos, counts=counts)
saveRDS(puck,here('ctSVGbench','real','puck',"myRCTD_Visium_liver.rds"))

## Visium_melanoma ---

obj=readRDS('~/yh/data/Visium_melanoma/melanoma_seurat_obj_filtered.rds') 
counts.sc <- obj@assays$RNA@counts
celltypes <- factor(obj$cell_type)
celltypes <- recode(celltypes,
                    "Monocyte/macrophage" = "Monocyte_macrophage")
celltypes <- recode(celltypes,
                    "T/NK cell" = "T_NK cell")                    
reference <- Reference(counts = counts.sc, cell_types = celltypes, min_UMI = 1) 
saveRDS(reference,here('ctSVGbench','real','reference',"myRCTD_Visium_melanoma.rds"))

obj=readRDS('~/yh/sc-ref/Visium_melanoma/melanoma_visium_sample02.rds')
counts <- obj@assays$Spatial@counts
table(obj@images$slice1@coordinates$tissue)
pos <- obj@images$slice1@coordinates[,c("imagerow","imagecol")]
colnames(pos) <- c('x','y')
puck <- SpatialRNA(coords=pos, counts=counts)
saveRDS(puck,here('ctSVGbench','real','puck',"myRCTD_Visium_melanoma.rds"))

## Visium_mousebrain ---
counts <- read.csv('~/yh/sc-ref/visium/st_count.csv', row.names=1)
pos <- read.csv('~/yh/sc-ref/visium/st_location.csv', row.names=1)
colnames(pos) <- c('x','y')
rownames(pos) <- gsub("-","\\.",rownames(pos))
puck <- SpatialRNA(coords=pos, counts=counts)
saveRDS(puck,here('ctSVGbench','real','puck',"myRCTD_Visium_mousebrain.rds"))

counts.sc <- read.csv('~/yh/sc-ref/visium/sc_mousebrain.csv',row.names=1)
celltypes <- read.csv('~/yh/sc-ref/visium/sc_celltype.csv',row.names=1)

Celltypes <- setNames(celltypes[,1],rownames(celltypes))
filter_index <- grep('LowQ',Celltypes)
Celltypes <- Celltypes[-filter_index]
Celltypes <- sapply(strsplit(Celltypes,'_'),"[[",1) 
Celltypes <- factor(Celltypes)

counts.sc <- t(counts.sc)
counts.sc <- counts.sc[,-filter_index]
dim(counts.sc)[2]==length(Celltypes)

reference <- Reference(counts = counts.sc, cell_types = Celltypes, min_UMI = 1)
saveRDS(reference,here('ctSVGbench','real','reference',"myRCTD_Visium_mousebrain.rds"))


## seqFish+_somatosensory ---
# counts <- t(read.csv('~/yh/sc-ref/seqFISH+/Out_gene_expressions_10000genes.csv',row.names=1))
# pos <- read.csv('~/yh/sc-ref/seqFISH+/Out_rect_locations.csv',row.names=1)
# pos <- pos[c("X","Y")]
# colnames(pos) <- c("x","y")
# puck <- SpatialRNA(coords=pos , counts=counts)
# saveRDS(puck,here('ctSVGbench','real','puck',"seqFish+_somatosensory.rds"))

# counts.sc <- read.table('~/yh/sc-ref/seqFISH+/raw_somatosensory_sc_exp.txt', header = TRUE)
# rownames(counts.sc) <- counts.sc$cell_id
# counts.sc$cell_id <- NULL
# celltypes <- read.table('~/yh/sc-ref/seqFISH+/somatosensory_sc_labels.txt')$V1
# celltypes <- setNames(celltypes,colnames(counts.sc))
# celltypes <- factor(celltypes)

# reference <- Reference(counts = counts.sc, cell_types = celltypes, min_UMI = 1) 
# saveRDS(reference,here('ctSVGbench','real','reference',"seqFish+_somatosensory.rds"))

### seqFish+_cortex --- 
counts <- t(read.csv('~/yh/sc-ref/seqFISH+_cortex+OB/cortex_svz_counts.csv'))
pos <- read.csv('~/yh/sc-ref/seqFISH+_cortex+OB/cortex_svz_cellcentroids.csv')
pos <- pos[c("X","Y")]
colnames(pos) <- c("x","y")
rownames(pos) <- colnames(counts) <- as.character(0:(nrow(pos)-1))
puck <- SpatialRNA(coords=pos , counts=counts)
saveRDS(puck,here('ctSVGbench','real','puck',"myRCTD_SeqFish+_mouse_cortex_svz.rds"))

prop <- read.csv('~/yh/sc-ref/seqFISH+_cortex+OB/cortex_svz_cell_type_annotations.csv',row.names=1)
prop$louvain <- as.character(prop$louvain)
Celltypes <- setNames(prop[,1],rownames(prop))

reference <- Reference(counts = counts, cell_types = factor(Celltypes), min_UMI = 1)
saveRDS(reference,here('ctSVGbench','real','reference',"myRCTD_SeqFish+_mouse_cortex_svz.rds"))



## seqFish+_OB ---
counts <- t(read.csv('~/yh/sc-ref/seqFISH+_cortex+OB/ob_counts.csv'))
pos <- read.csv('~/yh/sc-ref/seqFISH+_cortex+OB/ob_cellcentroids.csv')
pos <- pos[c("X","Y")]
colnames(pos) <- c("x","y")
rownames(pos) <- colnames(counts) <- as.character(0:(nrow(pos)-1))
puck <- SpatialRNA(coords=pos , counts=counts)
saveRDS(puck,here('ctSVGbench','real','puck',"myRCTD_SeqFish+_mouse_ob.rds"))

counts <- t(read.csv('~/yh/sc-ref/seqFISH+_cortex+OB/ob_counts.csv'))
prop <- read.csv('~/yh/sc-ref/seqFISH+_cortex+OB/OB_cell_type_annotations.csv',row.names=1)
prop$louvain <- as.character(prop$louvain)
Celltypes <- setNames(prop[,1],rownames(prop))
colnames(counts) <- as.character(0:(nrow(pos)-1))
reference <- Reference(counts = counts, cell_types = factor(Celltypes), min_UMI = 1)
saveRDS(reference,here('ctSVGbench','real','reference',"myRCTD_SeqFish+_mouse_ob.rds"))

#sc-ref=stereseq-mouse-OB

## ST_PDAC ---
obj=readRDS('~/yh/sc-ref/ST_PDAC/PDAC_GSM4100721.rds')
puck <- SpatialRNA(coords=obj$st_location , counts=obj$st_count)
saveRDS(puck,here('ctSVGbench','real','puck',"myRCTD_ST_PDAC.rds"))
celltypes <- setNames(obj$cell_type,colnames(obj$sc_count))
celltypes <- factor(celltypes)
celltypes <- recode(celltypes,
                    "T cells & NK cells " = "T cells and NK cells",
                    "Ductal - APOL1 high/hypoxic" = "Ductal APOL1 high_hypoxic",
                    "Ductal - CRISP3 high/centroacinar like" = "Ductal CRISP3 high_centroacinar like")
reference <- Reference(counts = obj$sc_count, cell_types = celltypes, min_UMI = 1) 
saveRDS(reference,here('ctSVGbench','real','reference',"myRCTD_ST_PDAC.rds"))



## ST_Developmentalheart ---
library(data.table)
file_path <- "~/yh/Developmentalheart/Filtered/share_files/all_cells_count_matrix_filtered.tsv.gz"
sc <- fread(file_path)
counts.sc <- data.frame(sc[,-1],row.names = sc[[1]])
meta_data <- read.delim("~/yh/Developmentalheart/Filtered/share_files/all_cells_meta_data_filtered.tsv.gz")
meta_data <- meta_data[match(colnames(counts.sc),meta_data$X),]
filter_doub <- which(meta_data$state=="Doublet")
meta_data <- meta_data[-filter_doub,]
counts.sc <- counts.sc[,-filter_doub]
celltypes <- setNames(meta_data$celltype,meta_data$X)
celltypes <- gsub("\\/"," or ",celltypes)

reference <- Reference(counts = counts.sc, cell_types = factor(celltypes), min_UMI = 1)
saveRDS(reference,here('ctSVGbench','real','reference',"myRCTD_ST_Developmentalheart.rds"))


st_path='~/yh/Developmentalheart/Filtered/filtered_ST_matrix_and_meta_data/filtered_matrix.tsv.gz'
st_meta='~/yh/Developmentalheart/Filtered/filtered_ST_matrix_and_meta_data/meta_data.tsv.gz'
st=fread(st_path)
meta=fread(st_meta)

rownames(st)=st$V1
st$V1 <- NULL

counts <- data.frame(st[,-1],row.names = st[[1]])
pos <- data.frame(x=meta$new_x,y=meta$new_y,row.names = meta$V1)
rownames(pos) <- paste0("X",rownames(pos))

library(biomaRt)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", version = 109)
ensembl_ids <- sub("\\.\\d+$", "", rownames(counts)) 
annot <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = ensembl_ids,
  mart = ensembl
)

rownames(counts) <- sub("\\.\\d+$", "", rownames(counts)) 
counts$ensembl_gene_id <- rownames(counts)
counts_merged <- merge(annot, counts, by = "ensembl_gene_id")
counts_merged <- counts_merged[!duplicated(counts_merged$hgnc_symbol),]
rownames(counts_merged) <- counts_merged$hgnc_symbol

counts_merged$ensembl_gene_id <- counts_merged$hgnc_symbol <- NULL
puck <- SpatialRNA(coords=pos, counts=counts_merged)
saveRDS(puck,here('ctSVGbench','real','puck',"myRCTD_ST_Developmentalheart.rds"))


## StereoSeq_mouse_olfactory bulb ---
obj=readRDS('~/yh/sc-ref/stereoseq_mouse/stereoseq_ob.RDS')
counts=t(obj$st_count)
puck <- SpatialRNA(coords=obj$st_location, counts=counts)
saveRDS(puck,here('ctSVGbench','real','puck',"myRCTD_StereoSeq_mouseOB.rds"))

obj=readRDS('~/yh/sc-ref/stereoseq_mouse/scRNA_seq_olfactory.RDS')
obj$cell_type_anno[,1] <- gsub("\\/"," or ",obj$cell_type_anno[,1])
celltypes <- setNames(obj$cell_type_anno[,1],rownames(obj$cell_type_anno))
celltypes <- factor(celltypes)

reference <- Reference(counts = obj$sc_count, cell_types = celltypes, min_UMI = 1) 
saveRDS(reference,here('ctSVGbench','real','reference',"myRCTD_StereoSeq_mouseOB.rds"))

# ## StereoSeq_zebrafish --- no raw counts
# counts <- t(read.csv('~/yh/sc-ref/stereoseq_zebrafish/ST_5.csv',row.names=1))
# pos <- read.csv('~/yh/sc-ref/stereoseq_zebrafish/ST_5_loc.csv',row.names=1)
# colnames(pos) <- c('x','y')
# puck <- SpatialRNA(coords=pos, counts=counts)
# saveRDS(puck,here('ctSVGbench','real','puck',"myRCTD_StereoSeq_zebrafish.rds"))


# counts.sc <- t(read.csv('~/yh/sc-ref/stereoseq_zebrafish/zf24_scRNA.csv',row.names=1))
# celltypes <- read.csv('~/yh/sc-ref/stereoseq_zebrafish/zf24_scRNA_type.csv',row.names=1)
# Celltypes <- setNames(celltypes[,1],rownames(celltypes))
# Celltypes <- sapply(strsplit(Celltypes,'_'),"[[",1) 

# Celltypes <- factor(Celltypes)
# reference <- Reference(counts = counts.sc, cell_types = Celltypes, min_UMI = 1)
# saveRDS(reference,here('ctSVGbench','real','reference',"myRCTD_StereoSeq_zebrafish.rds"))

## stereoseq_MBA ---
library(data.table)
library(Matrix)
st <- read.delim("~/yh/sc-ref/stereoseq_MBA/PFC_DT26.bin1.tsv.gz")
st <- data.table(st)

bin_size <- 50
st[, bin50_x := floor(x / bin_size)]
st[, bin50_y := floor(y / bin_size)]
st[, spot_bin50 := paste0(bin50_x, "_", bin50_y)]
st_sum <- st[, .(MIDCounts = sum(MIDCounts)), by = .(geneID, spot_bin50)]

gene_ids <- factor(st_sum$geneID)
spots <- factor(st_sum$spot_bin50)
expr_mat_bin50 <- sparseMatrix(
  i = as.integer(gene_ids),
  j = as.integer(spots),
  x = st_sum$MIDCounts,
  dimnames = list(levels(gene_ids), levels(spots))
)

coords <- unique(st[, .(spot_bin50, bin50_x, bin50_y)])
pos <- data.frame(x=coords$bin50_x,y=coords$bin50_y,row.names = coords$spot_bin50)
puck <- SpatialRNA(coords=pos, counts=expr_mat_bin50)
saveRDS(puck,here('ctSVGbench','real','puck',"myRCTD_StereoSeq_MBA.rds"))

sc <- readRDS('~/yh/sc-ref/stereoseq_MBA/snRNA-seq_macaque_brain_Seurat.RDS')
counts.sc <- Matrix::Matrix(round(sc@assays$SCT@counts), sparse = TRUE) 
celltypes <- sc$celltype
reference <- Reference(counts = counts.sc, cell_types = factor(celltypes), min_UMI = 1)
saveRDS(reference,here('ctSVGbench','real','reference',"myRCTD_StereoSeq_MBA.rds"))



## stereoseq_MDESTA ---
library(Seurat)
sc <- Read10X('~/yh/sc-ref/stereoseq_MDESTA/scRNAseq-1')
celltype <- read.csv('~/yh/sc-ref/stereoseq_MDESTA/scRNAseq-1/zm_celltype.csv',row.names=1)
barc=intersect(rownames(celltype),colnames(sc))
counts.sc <- sc[,barc]
celltype <- celltype[barc,]
celltypes <- setNames(celltype$cluster,rownames(celltype))
reference <- Reference(counts = counts.sc, cell_types = factor(celltypes), min_UMI = 1)
saveRDS(reference,here('ctSVGbench','real','reference',"myRCTD_StereoSeq_MDESTA.rds"))

library(data.table)
library(Matrix)
st <- fread('~/yh/sc-ref/stereoseq_MDESTA/MaizeEar_rep3_Raw_bin50.gem')
st[, spot := paste(x, y, sep = "_")]
gene_ids <- factor(st$geneID)
spots <- factor(st$spot)
counts <- sparseMatrix(
  i = as.integer(gene_ids),
  j = as.integer(spots),
  x = st$MIDCount,
  dimnames = list(levels(gene_ids), levels(spots))
) #long to wide
pos <- unique(st[, .(spot, x, y)])
pos <- as.data.frame(pos)
rownames(pos) <- pos$spot
pos <- pos[match(colnames(counts),rownames(pos)),]  
identical(rownames(pos),colnames(counts))
pos$spot <- NULL
puck <- SpatialRNA(coords=pos, counts=counts)
saveRDS(puck,here('ctSVGbench','real','puck',"myRCTD_StereoSeq_MDESTA.rds"))



## stereoseq_CBMSTA_Macaque ---
sc <- fread('~/yh/sc-ref/stereoseq_CBMSTA/Macaque_scRNA.csv')
ct <- read.csv('~/yh/sc-ref/stereoseq_CBMSTA/Macaque_scRNA_type.csv',row.names=1)
counts.sc <- t(data.frame(sc[,-1],row.names = sc[[1]]))
celltypes <- setNames(ct$ident,rownames(ct))
reference <- Reference(counts = counts.sc, cell_types = factor(celltypes), min_UMI = 1)
saveRDS(reference,here('ctSVGbench','real','reference',"myRCTD_StereoSeq_CBMSTA_Macaque.rds"))

st <- fread('~/yh/sc-ref/stereoseq_CBMSTA/Macaque_st.csv')
pos <- read.csv('~/yh/sc-ref/stereoseq_CBMSTA/Macaque_pos.csv',row.names=1)
counts <- data.frame(st[,-1],row.names = st[[1]])
rownames(pos) <- paste0("X",rownames(pos))
puck <- SpatialRNA(coords=pos, counts=counts)
saveRDS(puck,here('ctSVGbench','real','puck',"myRCTD_StereoSeq_CBMSTA_Macaque.rds"))

## stereoseq_CBMSTA_Marmoset ---
sc <- fread('~/yh/sc-ref/stereoseq_CBMSTA/Marmoset_scRNA.csv')
ct <- read.csv('~/yh/sc-ref/stereoseq_CBMSTA/Marmoset_scRNA_type.csv',row.names=1)
counts.sc <- t(data.frame(sc[,-1],row.names = sc[[1]]))
celltypes <- setNames(ct$ident,rownames(ct))
reference <- Reference(counts = counts.sc, cell_types = factor(celltypes), min_UMI = 1)
saveRDS(reference,here('ctSVGbench','real','reference',"myRCTD_StereoSeq_CBMSTA_Marmoset.rds"))

st <- fread('~/yh/sc-ref/stereoseq_CBMSTA/Marmoset_st.csv')
pos <- read.csv('~/yh/sc-ref/stereoseq_CBMSTA/Marmoset_pos.csv',row.names=1)
counts <- data.frame(st[,-1],row.names = st[[1]])
rownames(pos) <- paste0("X",rownames(pos))
puck <- SpatialRNA(coords=pos, counts=counts)
saveRDS(puck,here('ctSVGbench','real','puck',"myRCTD_StereoSeq_CBMSTA_Marmoset.rds"))

## stereoseq_CBMSTA_Mouse ---
sc <- fread('~/yh/sc-ref/stereoseq_CBMSTA/Mouse_scRNA.csv')
ct <- read.csv('~/yh/sc-ref/stereoseq_CBMSTA/Mouse_scRNA_type.csv',row.names=1)
counts.sc <- t(data.frame(sc[,-1],row.names = sc[[1]]))
celltypes <- setNames(ct$annotation,rownames(ct))
reference <- Reference(counts = counts.sc, cell_types = factor(celltypes), min_UMI = 1)
saveRDS(reference,here('ctSVGbench','real','reference',"myRCTD_StereoSeq_CBMSTA_Mouse.rds"))

st <- fread('~/yh/sc-ref/stereoseq_CBMSTA/Mouse_st.csv')
pos <- read.csv('~/yh/sc-ref/stereoseq_CBMSTA/Mouse_pos.csv',row.names=1)
counts <- data.frame(st[,-1],row.names = st[[1]])
rownames(pos) <- paste0("X",rownames(pos))
puck <- SpatialRNA(coords=pos, counts=counts)
saveRDS(puck,here('ctSVGbench','real','puck',"myRCTD_StereoSeq_CBMSTA_Mouse.rds"))

## stereoseq_CIRSTA ---
st=fread('~/yh/sc-ref/stereoseq_CIRSTA/st_counts.csv')
pos=fread('~/yh/sc-ref/stereoseq_CIRSTA/st_pos.csv')
counts <- data.frame(st[,-1],row.names = st[[1]])
coor <- data.frame(pos[,-1],row.names = pos[[1]])
rownames(coor) <- gsub("-","\\.",rownames(coor))
rownames(coor) <- paste0("X",rownames(coor))
index=grep("^X17D.M2.2_FS4", rownames(coor))
counts=counts[,index]
coor=coor[index,]
puck <- SpatialRNA(coords=coor, counts=counts)
saveRDS(puck,here('ctSVGbench','real','puck',"myRCTD_StereoSeq_CIRSTA.rds"))

sc=fread('~/yh/sc-ref/stereoseq_CIRSTA/scRNA.csv')
ct=fread('~/yh/sc-ref/stereoseq_CIRSTA/celltype.csv')
counts.sc <- t(data.frame(sc[,-1],row.names = sc[[1]]))
celltypes <- setNames(ct$annotation,ct$V1)
index=grep("D17.*M2", as.character(ct$V1),value = T)
counts.sc=counts.sc[,index]
celltypes=celltypes[index]
reference <- Reference(counts = counts.sc, cell_types = factor(celltypes), min_UMI = 1)
saveRDS(reference,here('ctSVGbench','real','reference',"myRCTD_StereoSeq_CIRSTA.rds"))

## mouse pancreatic_sc ---
library(Seurat)
obj2=Read10X('~/yh/sc-ref/panc_sc/GSE132188_RAW/GSM3852752')
obj3=Read10X('~/yh/sc-ref/panc_sc/GSE132188_RAW/GSM3852753')
obj4=Read10X('~/yh/sc-ref/panc_sc/GSE132188_RAW/GSM3852754')
obj5=Read10X('~/yh/sc-ref/panc_sc/GSE132188_RAW/GSM3852755')

ct=read.csv('~/yh/sc-ref/panc_sc/panc_celltype.csv',row.names=1)

colnames(obj2) <- paste0(colnames(obj2),"-0")
colnames(obj3) <- paste0(colnames(obj3),"-1")
colnames(obj4) <- paste0(colnames(obj4),"-2")
colnames(obj5) <- paste0(colnames(obj5),"-3")
gene=Reduce(intersect,list(rownames(obj2),rownames(obj3),rownames(obj4),rownames(obj5)))
obj=cbind(obj2[gene,],obj3[gene,],obj4[gene,],obj5[gene,])
umi=intersect(colnames(obj),rownames(ct))
ct=ct[umi,,drop=F]
obj=obj[,umi]

celltypes <- setNames(ct[,1],rownames(ct))
reference <- Reference(counts = obj, cell_types = factor(celltypes), min_UMI = 1)
saveRDS(reference,"~/yh/sc-ref/panc_sc/reference.rds")


## 
### Slide-seqV2_melanoma ---
library(data.table)
expr_mat <- fread("~/yh/sc-ref/SlideseqV2_melanoma/GSM6022264_MBM13_sn_counts.csv.gz")

gene_ids <- expr_mat[[1]]              
expr_numeric <- as.matrix(expr_mat[, -1])  
expr_numeric <- apply(expr_numeric, 2, function(x) as.numeric(trimws(x)))
rownames(expr_numeric) <- gene_ids
colnames(expr_numeric) <- colnames(expr_mat)[-1]
counts.sc <- expr_numeric

metadata=fread('~/yh/sc-ref/SlideseqV2_melanoma/GSE200218_sc_sn_metadata.csv.gz')
metadata_matched <- metadata[barcode_all %in% colnames(counts.sc), .(barcode_all, cell_type_main)]
celltypes <- metadata_matched[match(colnames(counts.sc), metadata_matched$barcode_all, nomatch = 0)]
celltypes <- setNames(make.names(celltypes$cell_type_main),celltypes$barcode_all)

reference <- Reference(counts = counts.sc, cell_types = factor(celltypes), min_UMI = 1)
saveRDS(reference,here('ctSVGbench','real','reference',"myRCTD_Slide-seqV2_melanoma.rds"))


counts <- fread("~/yh/sc-ref/SlideseqV2_melanoma/GSM6025944_MBM13_slide_raw_counts.csv.gz")
pos <- fread("~/yh/sc-ref/SlideseqV2_melanoma/GSM6025944_MBM13_slide_spatial_info.csv.gz")
counts <- as.data.frame(counts)
rownames(counts) <- counts[,1]
counts[,1] <- NULL

pos <- as.data.frame(pos[,2:4])
rownames(pos) <- pos[,1]
pos[,1] <- NULL
colnames(pos) <- c('x','y')
puck <- SpatialRNA(coords=pos, counts=counts)
saveRDS(puck,here('ctSVGbench','real','puck',"myRCTD_Slide-seqV2_melanoma.rds"))

### SlideseqV2_mouseOB ---
pos=fread("~/yh/sc-ref/SlideseqV2_mouseOB/GSM5173925_OB1_01_BeadLocationsForR.csv.gz")
pos=as.data.frame(pos)
rownames(pos) <- pos[,1]
pos[,1] <- NULL
colnames(pos) <- c('x','y')

obj=readRDS("~/yh/sc-ref/SlideseqV2_mouseOB/GSM5173925_OB1_Slide1.rds")
counts <- obj@assays$RNA@counts
puck <- SpatialRNA(coords=pos, counts=counts)
saveRDS(puck,here('ctSVGbench','real','puck',"myRCTD_Slide-seqV2_mouseOB.rds"))

##bash
cp myRCTD_SeqFish+_mouse_ob.rds myRCTD_Slide-seqV2_mouseOB.rds



### Visium_skin ---
data <- fread('~/yh/sc-ref/visium_skin/GSE144236_cSCC_counts.txt.gz')
setDF(data)  
rownames(data) <- data[[1]]
data <- data[ , -1]   


metadata <- fread('~/yh/sc-ref/visium_skin/GSE144236_patient_metadata_new.txt.gz')
setDF(metadata)
celltypes <- setNames(make.names(metadata$level1_celltype),metadata[,1])

reference <- Reference(counts = data, cell_types = factor(celltypes), min_UMI = 1)
saveRDS(reference,here('ctSVGbench','real','reference',"myRCTD_Visium_skin.rds"))


counts <- t(fread('~/yh/sc-ref/visium_skin/GSM4284316_P2_ST_rep1_stdata.tsv.gz'))
colnames(counts) <- as.character(unlist(counts[1, ]))
counts <- counts[-1, ]
counts <- as.matrix(counts)
storage.mode(counts) <- "numeric"

pos <- fread('~/yh/sc-ref/visium_skin/GSM4284316_spot_data-selection-P2_ST_rep1.tsv.gz')
pos[, spot_id := paste0(x, "x", y)]
pos <- as.data.frame(pos)
rownames(pos) <- pos$spot_id
counts <- counts[,match(pos$spot_id, colnames(counts))]

puck <- SpatialRNA(coords=pos[,1:2], counts=counts)
saveRDS(puck,here('ctSVGbench','real','puck',"myRCTD_Visium_skin.rds"))



### Visium_spleen ---
dir='~/yh/sc-ref/Visium_spleen/GSM6042726'
counts <- Read10X(data.dir = dir)
pos <- read.csv(gzfile("~/yh/sc-ref/Visium_spleen/GSM6042726_Round1_A1_Day7_Control_tissue_positions_list.csv.gz"),row.names=1)

pos1=pos[pos[,1]==1,2:3]
colnames(pos1) <- c('x','y')
spot <- intersect(rownames(pos1),colnames(counts))
pos1 <- pos1[spot,]
counts <- counts[,spot]
puck <- SpatialRNA(coords=pos1, counts=counts)
saveRDS(puck,here('ctSVGbench','real','puck',"myRCTD_Visium_spleen.rds"))

library(zellkonverter)
adata <- readH5AD("~/yh/sc-ref/Visium_spleen/tabula-muris-senis-facs-processed-official-annotations-Spleen.h5ad")
cell_annotations <- colData(adata)
celltypes <- data.frame(celltype=make.names(cell_annotations$cell_ontology_class),row.names =cell_annotations$cell )


library(tidyverse)
csv_files <- list.files("~/yh/sc-ref/Visium_spleen/Spleen", pattern = "\\.csv$", full.names = TRUE)
cell_data_list <- list()
for (file in csv_files) {
  cell_id <- tools::file_path_sans_ext(basename(file))
    df <- read.csv(file,header = F)
    cell_data_list[[cell_id]] <- setNames(df[,2], df[,1])
}

all_genes <- unique(unlist(lapply(cell_data_list, names)))

expr_matrix <- sapply(cell_data_list, function(x) {
  gene_counts <- rep(0, length(all_genes))
  names(gene_counts) <- all_genes
  gene_counts[names(x)] <- x
  gene_counts
})

counts.sc <- as.data.frame(expr_matrix)
rownames(counts.sc) <- all_genes
cells <- intersect(cell_annotations$cell,colnames(counts.sc))
celltypes <- celltypes[cells,,drop=F]
counts.sc <- counts.sc[,cells]
celltypes <- setNames(celltypes[,1],rownames(celltypes))
reference <- Reference(counts = counts.sc, cell_types = factor(celltypes), min_UMI = 1)
saveRDS(reference,here('ctSVGbench','real','reference',"myRCTD_Visium_spleen.rds"))

### Visium_bladder ---
library(zellkonverter)
adata <- readH5AD('~/yh/sc-ref/1/GSE169379_MIBC_snSeq.h5ad')
cell_annotations <- colData(adata)
celltypes <- setNames(cell_annotations$celltype,rownames(cell_annotations))
counts.sc <- SummarizedExperiment::assay(adata, "X")
reference <- Reference(counts = counts.sc, cell_types = factor(celltypes), min_UMI = 1)
saveRDS(reference,here('ctSVGbench','real','reference',"myRCTD_Visium_bladder.rds"))

adata <- readH5AD('~/yh/sc-ref/Visium_bladder/GSE171351_combined_visium.h5ad')
counts <- SummarizedExperiment::assay(adata, "X")
pos <- colData(adata)
pos1 <- pos[pos$Patient=='Bladder8',2:3]
puck <- SpatialRNA(coords=pos1, counts=counts)
saveRDS(puck,here('ctSVGbench','real','puck',"myRCTD_Visium_bladder.rds"))

### Visium_human_liver ---
cluster <- fread('~/yh/sc-ref/Visium_human_liver/E-MTAB-7407.clusters.tsv')
selected_K <- cluster[sel.K == TRUE, K]
cluster_assignments <- cluster[sel.K == TRUE, -c("sel.K", "K")]
clusters <- as.data.frame(t(cluster_assignments))
colnames(clusters) <- "cluster"
celltypes <- setNames(clusters$cluster,rownames(clusters))

path_prefix <- "~/yh/sc-ref/Visium_human_liver/E-MTAB-7407.aggregated_filtered_counts"
expr_mat <- readMM(paste0(path_prefix, ".mtx"))
gene_names <- fread(paste0(path_prefix, ".mtx_rows"), header = FALSE)[[1]]
spot_ids <- fread(paste0(path_prefix, ".mtx_cols"), header = FALSE)[[1]]

rownames(expr_mat) <- gene_names
colnames(expr_mat) <- spot_ids

reference <- Reference(counts = round(expr_mat), cell_types = factor(celltypes), min_UMI = 1)
saveRDS(reference,here('ctSVGbench','real','reference',"myRCTD_Visium_human_liver.rds"))

counts <- fread('~/yh/sc-ref/Visium_human_liver/GSM5093328_13-gan-gene_cell_exprs_table_symbol.txt.gz')
setDF(counts)
rownames(counts) <- counts$Symbol
counts$Symbol <- NULL
puck <- SpatialRNA(coords=pos, counts=counts)
saveRDS(puck,here('ctSVGbench','real','puck',"myRCTD_Visium_human_liver.rds"))


### Visium_spinal  --- 
metadata <- fread('~/yh/sc-ref/Visium_spinal /GSE184370_meta.txt.gz')
celltypes <- setNames(metadata$global_cell_type,metadata$barcode)

prefix <- "~/yh/sc-ref/Visium_spinal/GSE184370/"
expr_mat <- readMM(paste0(prefix, "GSE184370_courtine_lumbar-EES-rehab_UMI.mtx.gz"))
features <- fread(paste0(prefix, "GSE184370_features.txt.gz"), header = FALSE)
genes <- make.unique(features$V1) 
rownames(expr_mat) <- genes

barcodes <- fread(paste0(prefix, "GSE184370_barcodes.txt.gz"), header = FALSE)
colnames(expr_mat) <- barcodes$V1
reference <- Reference(counts = expr_mat, cell_types = factor(celltypes), min_UMI = 1)
saveRDS(reference,here('ctSVGbench','real','reference',"myRCTD_Visium_spinal.rds"))

pos <- as.data.frame(fread('~/yh/sc-ref/Visium_spinal/GSE184369_meta.txt.gz'))[,c('barcode','coord_x','coord_y','spot')]
barcodes_1_B1 <- pos$barcode[pos$spot=='Slide1_1_B1']
rownames(pos) <- pos$barcode
pos <- pos[barcodes_1_B1,2:3]

counts <- readMM('~/yh/sc-ref/Visium_spinal/GSE184369_courtine_lumbar-EES-rehab-spatial_UMI.mtx.gz')
features <- fread('~/yh/sc-ref/Visium_spinal/GSE184369_features.txt.gz',header=F)[[1]]
barcodes <- fread('~/yh/sc-ref/Visium_spinal/GSE184369_barcodes.txt.gz',header=F)[[1]]
colnames(counts) <- barcodes
rownames(counts) <- features
dim(counts)
counts <- counts[,barcodes_1_B1]
puck <- SpatialRNA(coords=pos, counts=counts)
saveRDS(puck,here('ctSVGbench','real','puck',"myRCTD_Visium_spinal.rds"))

### Visium_intestine --- 
anno <- fread('~/yh/sc-ref/Visium_intestine/GSE132465_GEO_processed_CRC_10X_cell_annotation.txt.gz')
celltypes <- setNames(make.names(anno$Cell_type) ,anno$Index)

counts.sc <- fread('~/yh/sc-ref/Visium_intestine/GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt.gz')
setDF(counts.sc)
rownames(counts.sc) <- counts.sc$Index
counts.sc$Index <- NULL
reference <- Reference(counts = counts.sc, cell_types = factor(celltypes), min_UMI = 1)
saveRDS(reference,here('ctSVGbench','real','reference',"myRCTD_Visium_intestine.rds"))

counts <- Read10X('~/yh/sc-ref/Visium_intestine/filtered_feature_bc_matrix')
pos <- read.csv('~/yh/sc-ref/Visium_intestine/spatial/tissue_positions_list.csv',header = F)
rownames(pos) <- pos[,1]
pos1=pos[pos$V2==1,3:4]
colnames(pos1) <- c('x','y')
spot <- intersect(rownames(pos1),colnames(counts))
pos1 <- pos1[spot,]
counts <- counts[,spot]
puck <- SpatialRNA(coords=pos1, counts=counts)
saveRDS(puck,here('ctSVGbench','real','puck',"myRCTD_Visium_intestine.rds"))

### Visium_tail ---
counts.sc <- Read10X('~/yh/sc-ref/Visium_tail/GSM7840723')
library(readxl)
metadata <- read_excel("~/yh/sc-ref/Visium_tail/GSM7840723_processed_metadata.xlsx")
metadata <- as.data.frame(metadata)
metadata[,1] <- sapply(strsplit(metadata[,1],"-"),"[[",1)
celltypes <- setNames(make.names(metadata$cell_type_annotation),metadata[,1])
reference <- Reference(counts = counts.sc, cell_types = factor(celltypes), min_UMI = 1)
saveRDS(reference,here('ctSVGbench','real','reference',"myRCTD_Visium_tail.rds"))

counts <- Read10X('~/yh/sc-ref/Visium_tail/GSM7840725')
pos <- fread('~/yh/sc-ref/Visium_tail/GSM7840725_B1_frog_6_24hpa_tissue_positions_list.csv.gz')
setDF(pos)
rownames(pos) <- pos[,1]
pos1=pos[pos$V2==1,3:4]
colnames(pos1) <- c('x','y')
spot <- intersect(rownames(pos1),colnames(counts))
pos1 <- pos1[spot,]
counts <- counts[,spot]
puck <- SpatialRNA(coords=pos1, counts=counts)
saveRDS(puck,here('ctSVGbench','real','puck',"myRCTD_Visium_tail.rds"))


### Visium_pancreas --- 
counts.sc <- Read10X('~/yh/sc-ref/Visium_pancreas/GSE243466')
metadata <- fread('~/yh/sc-ref/Visium_pancreas/GSE243466_Metadata_df.csv.gz')
metadata <- as.data.frame(metadata)
celltypes <- setNames(make.names(metadata$ManuscriptClustersByMaturation),metadata[,1])
reference <- Reference(counts = counts.sc, cell_types = factor(celltypes), min_UMI = 1)
saveRDS(reference,here('ctSVGbench','real','reference',"myRCTD_Visium_pancreas.rds"))

dir <- '~/yh/sc-ref/Visium_pancreas/XV0008C1'
seurat_obj <- Seurat::Load10X_Spatial(data.dir = dir)
list.files(dir)
list.files(dir, recursive = TRUE)

counts <- Read10X('~/yh/sc-ref/Visium_pancreas/XV0008C1/filtered_feature_bc_matrix')
pos <- read.csv('~/yh/sc-ref/Visium_pancreas/XV0008C1/spatial/XV0008C1_tissue_positions_list.csv',header = F)
rownames(pos) <- pos[,1]
pos1=pos[pos$V2==1,3:4]
colnames(pos1) <- c('x','y')
spot <- intersect(rownames(pos1),colnames(counts))
pos1 <- pos1[spot,]
counts <- counts[,spot]
puck <- SpatialRNA(coords=pos1, counts=counts)
saveRDS(puck,here('ctSVGbench','real','puck',"myRCTD_Visium_pancreas.rds"))


##slide_seq2_melanoma
expr_mat <- readMM("/home/user/Fanglab1/yh/sc-ref/slide_seq2_melanoma/GSE200218_sc_sn_counts.mtx.gz")

dim(expr_mat)
gene_names <- fread(
  input = "/home/user/Fanglab1/yh/sc-ref/slide_seq2_melanoma/GSE200218_sc_sn_gene_names.csv.gz",
  header = FALSE  
)
gene_vec <- as.character(gene_names[[1]])  
rownames(expr_mat) <- gene_vec[-1]
colnames(expr_mat) <- metadata[[1]]

metadata=fread('/home/user/Fanglab1/yh/sc-ref/slide_seq2_melanoma/GSE200218_sc_sn_metadata.csv.gz')
dim(metadata)
df_mpm <- metadata[patient %like% "MPM", c("barcode_all", "cell_type_main")]
df_mbm <- metadata[patient %like% "MBM", c("barcode_all", "cell_type_main")]

dim(df_mpm)
counts.sc.mpm <- expr_mat[,df_mpm$barcode_all]
counts.sc.mpm <- counts.sc.mpm[!duplicated(rownames(counts.sc.mpm)),]
celltypes_mpm <- setNames(make.names(df_mpm$cell_type_main),df_mpm$barcode_all)
reference.mpm <- Reference(counts = counts.sc.mpm, cell_types = factor(celltypes_mpm), min_UMI = 1)
saveRDS(reference.mpm,here('ctSVGbench','real','reference',"myRCTD_Slide-seqV2_melanoma_MPM.rds"))


dim(df_mbm)
counts.sc.mbm <- expr_mat[,df_mbm$barcode_all]
counts.sc.mbm <- counts.sc.mbm[!duplicated(rownames(counts.sc.mbm)),]
celltypes_mbm <- setNames(make.names(df_mbm$cell_type_main),df_mbm$barcode_all)
reference.mbm <- Reference(counts = counts.sc.mbm, cell_types = factor(celltypes_mbm), min_UMI = 1)
saveRDS(reference.mbm,here('ctSVGbench','real','reference',"myRCTD_Slide-seqV2_melanoma_MBM.rds"))
celltypes_mpm[1:3]

###
library(data.table)
library(here)

process_single_sample <- function(sample_id, data_dir, output_dir) {
  
  counts_file <- file.path(data_dir, paste0(sample_id, "_slide_raw_counts.csv.gz"))
  pos_file <- file.path(data_dir, paste0(sample_id, "_slide_raw_counts.csv.gz")) %>% 
    gsub("raw_counts", "spatial_info", .) 
  
  if (!file.exists(counts_file) || !file.exists(pos_file)) {
    return(FALSE)
  }
  
  counts <- fread(counts_file)
  counts <- as.data.frame(counts)
  rownames(counts) <- counts[, 1]  
  counts[, 1] <- NULL  
  
  pos <- fread(pos_file)
  pos <- as.data.frame(pos[, 2:4])  
  rownames(pos) <- pos[, 1] 
  pos[, 1] <- NULL  
  colnames(pos) <- c('x', 'y')  
  
  puck <- SpatialRNA(coords = pos, counts = counts)
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  output_rds <- file.path(output_dir, paste0("myRCTD_Slide-seqV2_melanoma_", sample_id, ".rds"))
  saveRDS(puck, output_rds)
  
  message(paste0("sample : ", sample_id, " done", output_rds))
  return(TRUE)
}

root_dir = "/home/user/Fanglab1/yh/sc-ref/slide_seq2_melanoma"
gz_files <- list.files(root_dir, pattern = "\\.gz$", full.names = FALSE)
sample_ids <- gz_files %>% 
  gsub("_slide.*\\.gz$", "", .) %>%  
  unique()  
output_dir <- here('ctSVGbench', 'real', 'puck')

lapply(sample_ids, function(sid) {
  process_single_sample(
    sample_id = sid,
    data_dir = root_dir,
    output_dir = output_dir
  )})


##### stereoseq_mosta
library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)

process_stereo_gem <- function(gem_gz_path) {
  file_name <- tools::file_path_sans_ext(tools::file_path_sans_ext(basename(gem_gz_path)))
  out_dir <- file.path(dirname(gem_gz_path), file_name)
  
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  gem_data <- read_tsv(
    file = gem_gz_path,
    show_col_types = FALSE
  )
  
  cell_col <- ifelse("cell_label" %in% colnames(gem_data), "cell_label", "cell_id")

umi_col <- ifelse("total_umi" %in% colnames(gem_data), "total_umi", "umi_count")

valid_data <- gem_data 
  filter(!!sym(cell_col) != 0)

counts <- valid_data %>%
  group_by(!!sym(cell_col), gene) %>%
  summarise(umi_value = sum(!!sym(umi_col)), .groups = "drop") %>%
  pivot_wider(names_from = gene, values_from = umi_value, values_fill = 0) %>%
  column_to_rownames(cell_col)
  
  pos <- valid_data %>%
    group_by(!!sym(cell_col)) %>%
    summarise(
      x = mean(x),
      y = mean(y),
      .groups = "drop"
    ) %>%
    column_to_rownames(cell_col)

  saveRDS(t(counts), file = file.path(out_dir, "counts.rds"))
  saveRDS(pos, file = file.path(out_dir, "pos.rds"))

  puck <- SpatialRNA(coords=pos, counts=counts)
  puck_save_root <- "/home/user/Fanglab1/yh/ctSVGbench/real/puck/"
  puck_save_file <- file.path(puck_save_root, paste0("myRCTD_", file_name, ".rds"))
  saveRDS(puck, file = puck_save_file)
}


files <- list.files("/home/user/Fanglab1/yh/sc-ref/stereoseq_cbmsta",pattern = ".gz",full.names = T)

for (f in files){
 process_stereo_gem(gem_gz_path = f)

}

###E16.5_E1S3
library(Matrix)   
library(readr)    
target_dir <- "/home/user/Fanglab1/yh/sc-ref/stereoseq_mosta/E16.5_E1S3"
counts_path <- file.path(target_dir, "counts.mtx")
annotation_path <- file.path(target_dir, "annotation.csv")
position_path <- file.path(target_dir, "position.csv")

counts_matrix <- readMM(counts_path)
gene_names <- readLines(file.path(target_dir, "gene_names.txt"))
cell_names <- readLines(file.path(target_dir, "cell_names.txt"))
dim(counts_matrix)
length(gene_names)
length(cell_names)

colnames(counts_matrix) <- gene_names
rownames(counts_matrix) <- cell_names

head(rownames(counts_matrix))
head(colnames(counts_matrix))

annotation_df <- read.csv(
  file = annotation_path,
  header = TRUE,  
  stringsAsFactors = FALSE  
)
head(annotation_df)

celltypes <- setNames(annotation_df$annotation, annotation_df$cell_id)
prop <- sapply(sort(unique(celltypes)), function(x) as.numeric(celltypes == x))
dimnames(prop) <- list(names(celltypes), sort(unique(celltypes)))
prop <- as.data.frame(prop)


position_df <- read.csv(
  file = position_path,
  header = TRUE,
  stringsAsFactors = FALSE
)
head(position_df)

pos=data.frame(x=position_df$x,y=position_df$y,row.names = position_df$cell_id)
counts <- t(counts_matrix)
dim(counts)
mylist <- list(counts=counts,pos=pos,prop=prop)
saveRDS(mylist,file = here('real','input_sc',"stereoseq_mosta_E16.5_E1S3.rds"))


###E16.5_E1S3_whole_brain    
target_dir <- "/home/user/Fanglab1/yh/sc-ref/stereoseq_mosta/E16.5_E1S3_whole_brain"
counts_path <- file.path(target_dir, "counts.mtx")
annotation_path <- file.path(target_dir, "annotation.csv")
position_path <- file.path(target_dir, "position.csv")

counts_matrix <- readMM(counts_path)
gene_names <- readLines(file.path(target_dir, "gene_names.txt"))
cell_names <- readLines(file.path(target_dir, "cell_names.txt"))
dim(counts_matrix)
length(gene_names)
length(cell_names)

colnames(counts_matrix) <- gene_names
rownames(counts_matrix) <- cell_names

head(rownames(counts_matrix))
head(colnames(counts_matrix))

annotation_df <- read.csv(
  file = annotation_path,
  header = TRUE,  
  stringsAsFactors = FALSE  
)
head(annotation_df)

celltypes <- setNames(annotation_df$annotation, annotation_df$cell_id)
prop <- sapply(sort(unique(celltypes)), function(x) as.numeric(celltypes == x))
dimnames(prop) <- list(names(celltypes), sort(unique(celltypes)))
prop <- as.data.frame(prop)


position_df <- read.csv(
  file = position_path,
  header = TRUE,
  stringsAsFactors = FALSE
)
head(position_df)

pos=data.frame(x=position_df$x,y=position_df$y,row.names = position_df$cell_id)
counts <- t(counts_matrix)
dim(counts)
mylist <- list(counts=counts,pos=pos,prop=prop)
saveRDS(mylist,file = here('real','input_sc',"stereoseq_mosta_E16.5_E1S3_whole_brain.rds"))


###Dorsal_midbrain
target_dir <- "/home/user/Fanglab1/yh/sc-ref/stereoseq_mosta/Dorsal_midbrain"
counts_path <- file.path(target_dir, "counts.mtx")
annotation_path <- file.path(target_dir, "annotation.csv")
position_path <- file.path(target_dir, "position.csv")

counts_matrix <- readMM(counts_path)
gene_names <- readLines(file.path(target_dir, "gene_names.txt"))
cell_names <- readLines(file.path(target_dir, "cell_names.txt"))
dim(counts_matrix)
length(gene_names)
length(cell_names)

colnames(counts_matrix) <- gene_names
rownames(counts_matrix) <- cell_names

head(rownames(counts_matrix))
head(colnames(counts_matrix))

annotation_df <- read.csv(
  file = annotation_path,
  header = TRUE,  
  stringsAsFactors = FALSE  
)
head(annotation_df)

celltypes <- setNames(annotation_df$annotation, annotation_df$cell_id)
prop <- sapply(sort(unique(celltypes)), function(x) as.numeric(celltypes == x))
dimnames(prop) <- list(names(celltypes), sort(unique(celltypes)))
prop <- as.data.frame(prop)


position_df <- read.csv(
  file = position_path,
  header = TRUE,
  stringsAsFactors = FALSE
)
head(position_df)

pos=data.frame(x=position_df$x,y=position_df$y,row.names = position_df$cell_id)
counts <- t(counts_matrix)
dim(counts)
mylist <- list(counts=counts,pos=pos,prop=prop)
saveRDS(mylist,file = here('real','input_sc',"stereoseq_mosta_Dorsal_midbrain.rds"))
