import pandas as pd
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)

#lymphnode_sc
sce = sc.read('~/yh/sc.h5ad')
sc.pp.filter_genes(sce, min_cells=3)
expr = sce.X.toarray() 
expr_df = pd.DataFrame(expr, index=sce.obs_names, columns=sce.var_names)
expr_df.to_csv("~/yh/data/lymphnode_sc/counts.sc.csv")
sce.obs[['new_celltype']].to_csv("~/yh/data/lymphnode_sc/celltype.csv")


#stereoseq_zebrafish
sce=sc.read('~/yh/sc-ref/stereoseq_zebrafish/zf24_scRNA.h5ad')
expr = sce.X.toarray() 
expr_df = pd.DataFrame(expr, index=sce.obs_names, columns=sce.var_names)
expr_df.to_csv("~/yh/sc-ref/stereoseq_zebrafish/zf24_scRNA.csv")
sce.obs[['celltype_new']].to_csv("~/yh/sc-ref/stereoseq_zebrafish/zf24_scRNA_type.csv")


## stereoseq_CIRSTA ---
sce=sc.read('~/yh/sc-ref/allStage.scRNAseq.expression.h5ad')
sce.obs.to_csv("~/yh/sc-ref/stereoseq_CIRSTA/celltype.csv")
expr = sce.X.toarray() 
expr_df = pd.DataFrame(expr, index=sce.obs_names, columns=sce.var_names)
expr_df.to_csv("~/yh/sc-ref/stereoseq_CIRSTA/scRNA.csv")

st=sc.read('~/yh/sc-ref/stereoseq_CIRSTA/allStage_display.expression.h5ad')
counts=st.X.toarray().astype(np.int32)
expr_df = pd.DataFrame(counts.T, index=st.var_names, columns=st.obs_names)
expr_df.to_csv("~/yh/sc-ref/stereoseq_CIRSTA/st_counts.csv")

pos=st.obsm['spatial']
spot_names = st.obs_names 
pos_df = pd.DataFrame(pos, columns=["x", "y"], index=spot_names)
pos_df.iloc[1:3,]
pos_df.to_csv("~/yh/sc-ref/stereoseq_CIRSTA/st_pos.csv")




## pancreatic endocrinogenesis ---
sce=sc.read('~/yh/GSE132188_adata.h5ad.h5')
sce.obs.clusters_fig2_final.to_csv("~/yh/sc-ref/panc_celltype.csv")
counts_sc=sce.X.toarray()


## #stereoseq_CBMSTA_Macaque ---
sce=sc.read('~/yh/sc-ref/stereoseq_CBMSTA/Macaque.sn.h5ad')
sce.obs[['ident']].to_csv("~/yh/sc-ref/stereoseq_CBMSTA/Macaque_scRNA_type.csv")
expr = sce.X.toarray() 
expr_int = expr.astype(np.int32)
expr_df = pd.DataFrame(expr_int, index=sce.obs_names, columns=sce.var_names)
expr_df.to_csv("~/yh/sc-ref/stereoseq_CBMSTA/Macaque_scRNA.csv")

st=sc.read('~/yh/sc-ref/stereoseq_CBMSTA/Macaque1_T40.h5ad')
counts_st=st.X.toarray()
counts = counts_st.astype(np.int32)
expr_df = pd.DataFrame(counts.T, index=st.var_names, columns=st.obs_names)
expr_df.to_csv("~/yh/sc-ref/stereoseq_CBMSTA/Macaque_st.csv")

pos=st.obsm['spatial']
spot_names = st.obs_names 
pos_df = pd.DataFrame(pos, columns=["x", "y"], index=spot_names)
pos_df.to_csv("~/yh/sc-ref/stereoseq_CBMSTA/Macaque_pos.csv")


## #stereoseq_CBMSTA_Marmoset ---
sce=sc.read('~/yh/sc-ref/stereoseq_CBMSTA/Marmoset.sn.h5ad')
sce.obs[['ident']].to_csv("~/yh/sc-ref/stereoseq_CBMSTA/Marmoset_scRNA_type.csv")
expr = sce.X.toarray() 
expr_int = expr.astype(np.int32)
expr_df = pd.DataFrame(expr_int, index=sce.obs_names, columns=sce.var_names)
expr_df.to_csv("~/yh/sc-ref/stereoseq_CBMSTA/Marmoset_scRNA.csv")

st=sc.read('~/yh/sc-ref/stereoseq_CBMSTA/Marmoset1_T478.h5ad')
counts_st=st.X.toarray()
counts = counts_st.astype(np.int32)
expr_df = pd.DataFrame(counts.T, index=st.var_names, columns=st.obs_names)
expr_df.to_csv("~/yh/sc-ref/stereoseq_CBMSTA/Marmoset_st.csv")

pos=st.obsm['spatial']
spot_names = st.obs_names 
pos_df = pd.DataFrame(pos, columns=["x", "y"], index=spot_names)
pos_df.to_csv("~/yh/sc-ref/stereoseq_CBMSTA/Marmoset_pos.csv")

# ===== Lung LUSC =====
import scanpy as sc
import pandas as pd
import numpy as np
import csv
from scipy import sparse

file_path = "./a4dbbb30-3d3d-4760-b6d3-bc899f748cf7.h5ad"

# 1. Load data in backed mode (low memory usage)
adata = sc.read_h5ad(file_path, backed="r")

# 2. Filter cells belonging to the target disease
disease_name = "squamous cell lung carcinoma"
sq_cells = adata.obs_names[adata.obs['disease'] == disease_name]

# 3. Create a subset and save to disk (memory-efficient)
adata_sq = adata[sq_cells, :].copy(filename="sq_cells.h5ad")

# 4. Load the subset into memory
adata_sq = sc.read_h5ad("sq_cells.h5ad").to_memory()

# 5. Export cell annotations
annotations_df = adata_sq.obs
annotations_df.to_csv("sq_cells_annotations.csv", index=True)

# 6. Export count matrix
if hasattr(adata_sq.X, "toarray"):
    counts_df = pd.DataFrame(
        adata_sq.raw.X.toarray(),
        index=adata_sq.raw.obs_names,
        columns=adata_sq.raw.var_names
    )
else:
    counts_df = pd.DataFrame(
        adata_sq.raw.X,
        index=adata_sq.raw.obs_names,
        columns=adata_sq.raw.var_names
    )

counts_df.to_csv("sq_cells_counts.csv", index=True)


# ===== Lung LUAD =====
import scanpy as sc
import pandas as pd
import numpy as np
from anndata import AnnData
import csv
from scipy import sparse

file_path = "./a4dbbb30-3d3d-4760-b6d3-bc899f748cf7.h5ad"

# 1. Load data in backed mode (low memory usage)
adata = sc.read_h5ad(file_path, backed="r")

# 2. Filter cells belonging to the target disease
disease_name = "lung adenocarcinoma"
ad_cells = adata.obs_names[adata.obs['disease'] == disease_name]

# 3. Create a subset and save to disk (memory-efficient)
adata_ad = adata[ad_cells, :].copy(filename="ad_cells.h5ad")

# 4. Load the subset into memory
adata_ad = sc.read_h5ad("ad_cells.h5ad").to_memory()

# 5. Export cell annotations
annotations_df = adata_ad.obs
annotations_df.to_csv("ad_cells_annotations.csv", index=True)

# 6. Export count matrix (with row and column names)
counts_with_names = pd.DataFrame(
    data=adata_ad.raw.X.toarray(),
    index=adata_ad.raw.obs_names,   # cell names
    columns=adata_ad.raw.var_names  # gene names
)

# Save to CSV
counts_with_names.to_csv("counts_matrix_with_names_ad.csv", index=True)
