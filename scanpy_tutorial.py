import numpy as np
import pandas as pd
import scanpy as sc

def unique_clusters(zdata):
    '''Find all unique clusters from AnnData object
        and return a list.'''
    clust_list = []
    for cat in zdata.obs['leiden']:
        if cat not in clust_list:
            clust_list.append(cat)
    return clust_list

def parse_clusters(zdata, cluster_list):
    '''Based on a list of all clusters in the AnnData object,
        make new objects which are the subset of cells corresponding
        to each cluster. Returns a dictionary of objs keyed by cluster ID.'''
    clust_dict = {}
    for cluster in cluster_list:
        tmpdata = zdata[zdata.obs.leiden == cluster]
        clust_dict[cluster] = tmpdata
    return clust_dict

# Set up general settings
sc.settings.verbosity = 3
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')

results_file = 'output/tutorial_WT-0-7.h5ad'  # output of preprocessing
                                              # clustering

# Read count matrix into AnnData object
adata = sc.read_10x_mtx(
    'raw/',                    # directory with the .mtx file
    var_names='gene_symbols',  # use gene symbols for var names
    cache=True)                # cache file for faster reading

adata.var_names_make_unique()

print(adata)

# PREPROCESSING
# show genes that yield highest fraction of counts in each cell, for all
sc.pl.highest_expr_genes(adata, n_top=20, )

# basic filtering
sc.pp.filter_cells(adata, min_genes=200) # filter cells w < 200 genes
sc.pp.filter_genes(adata, min_cells=3)   # filter genes found in < 3 cells

# compute mitochondrial gene metrics
adata.var['mt'] = adata.var_names.str.startswith('mt-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None,
                           log1p=False, inplace=True)

# violin plots of some computer quality measures
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)

# remove cells with too many mt genes express or too many total counts
sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt')
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')

adata = adata[adata.obs.n_genes_by_counts < 2500, :]
adata = adata[adata.obs.pct_counts_mt < 5, :]

# total-count normalize data matrix X to 10,000 reads per cell
sc.pp.normalize_total(adata, target_sum=1e4) # make counts comparable b/w cells

# log transform the data
sc.pp.log1p(adata)

# identify highly-variable genes
sc.pp.highly_variable_genes(adata, min_mean=0.0125,
                            max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(adata)

# set .raw attribute of AnnData obj to the normalized & logarithmized data
# for later use in differential testing/visualizations; freezes state of obj
# can get back an Anndata of the object in .raw by calling .raw.to_adata()
adata.raw = adata

# filter based on highly variable genes
adata = adata[:, adata.var.highly_variable]

# regress out the effects of total counts per cell and %mt genes expressed
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])

# scale data to unit variance; clip values > std dev 10
sc.pp.scale(adata, max_value=10)

# PRINCIPAL COMPONENT ANALYSIS (PCA)
# reduces dimensionality of the data; reveals main axes of variation & denoises
sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca(adata, color='Cxcr6')

# inspect contribution of PCs to total variance in data
sc.pl.pca_variance_ratio(adata, log=True)

#save result
adata.write(results_file)
print(adata)

# COMPUTING THE NEIGHBORHOOD GRAPH
# these values based on Seurat tutorial
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=20)

# EMBEDDING THE NEIGHBORHOOD GRAPH
# embed the graph in two dimensions using UMAP
sc.tl.umap(adata)
sc.pl.umap(adata, color=['Cxcr6', 'Cxcr4', 'Ccr7'])
#sc.tl.paga(adata)
#sc.pl.paga(adata, plot=False)
#sc.tl.umap(adata, init_pos='paga')

# CLUSTERING THE NEIGHBORHOOD GRAPH
# directly clusters neighborhood graph of cells from prev. section
sc.tl.leiden(adata, resolution=0.7)
sc.pl.umap(adata, color=['leiden', 'Cxcr6', 'Cxcr4'])

# HEATMAP OF MARKER GENES
# 10 most variable genes per cluster
# find all clusters
clusters = unique_clusters(adata)

# create AnnData obj of each cluster & store as dict keyed by cluster
data_clusters = parse_clusters(adata, clusters)

# for each cluster, find 10 highest-variance genes


adata.write(results_file)