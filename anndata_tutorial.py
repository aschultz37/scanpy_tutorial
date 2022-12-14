import numpy as np
import pandas as pd
import anndata as ad
from scipy.sparse import csr_matrix

#AnnData tutorial
print(ad.__version__)

#Create basic AD object w sparse count information (e.g. gene expression
#counts)
counts = csr_matrix(np.random.poisson(1, size=(100,2000)), dtype=np.float32)
adata = ad.AnnData(counts)
print(adata)

##initial data we passed, stored as a sparse matrix in adata.X
print(adata.X)

#Provide index for the obs and var axes
adata.obs_names = [f"Cell_{i:d}" for i in range(adata.n_obs)]
adata.var_names = [f"Gene_{i:d}" for i in range(adata.n_vars)]
print(adata.obs_names[:10])

#Add metadata at obs and var levels
##N.B. adata.obs and adata.var are Pandas DataFrames
ct = np.random.choice(["B", "T", "Monocyte"], size=(adata.n_obs,))
adata.obs["cell_type"] = pd.Categorical(ct) #Categoricals are more efficient
print(adata.obs)

##see that AnnData object has been updated
print(adata)

#Can subset the AnnData using the randomly selected cell types (ct)
bdata = adata[adata.obs.cell_type == "B"]
print(bdata)

#Can have multi-dimensional metadata (e.g. UMAP embedding)
#Use .obsm or .varm for multi-dimensional metadata
#.obsm must == .n_obs, and .varm must == .n_vars

##random matrix that we can interpret as UMAP embedding
adata.obsm["X_umap"] = np.random.normal(0, 1, size=(adata.n_obs, 2))
##random gene-level metadata
adata.varm["gene_stuff"] = np.random.normal(0, 1, size=(adata.n_vars, 5))
print(adata.obsm)

##see that AnnData object has been updated
print(adata)

#Notes on .obsm/.varm
#The array-like metadata can originate from pd.dataframe, scipy.sparse_matrix,
#or np.dense_array
#In scanpy, their values (columns) are not easily plotted, whereas items from
#.obs are easily plotted on things such as UMAP plots

#Unstructured metadata (.uns) can be anything (e.g. list, dict)
adata.uns["random"] = [1,2,3]
print(adata.uns)

#Can store multiple forms of the core data in Layers
adata.layers["log_transformed"] = np.log1p(adata.X)
print(adata)

##can return a DataFrame of a layer
##.obs_names --> row, .var_names --> column
layer_df = adata.to_df(layer="log_transformed")
print(layer_df)

#Writing results to disk
adata.write("tutorial_results.h5ad", compression="gzip")

#