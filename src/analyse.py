from numpy.random import default_rng
import matplotlib.pyplot as plt
import scanpy as sc
import squidpy as sq
from anndata import AnnData
import argparse

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Load AnnData, run spatial analysis, and save results")
parser.add_argument("--input", type=str, required=True, help="Path to the input AnnData (.h5ad) file")
parser.add_argument("--output", type=str, required=True, help="Path to save the processed AnnData (.h5ad) file")
parser.add_argument("--prefix", type=str, default="spatial_", help="Prefix for output plots")
args = parser.parse_args()

# Load AnnData
adata = sc.read(args.input)
print(adata)

# FIX GENE NAMES
if adata.var_names.dtype in ['int64', 'int32']:
    if 'gene_symbols' in adata.var.columns:
        adata.var_names = adata.var['gene_symbols'].values
    else:
        adata.var_names = [f"gene_{i}" for i in range(adata.n_vars)]

adata.var_names_make_unique()

# -------------------------------
# QC (RAW COUNTS - CORRECT PLACE)
# -------------------------------
adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

# -------------------------------
# PREPROCESSING
# -------------------------------
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)

sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000)
sc.pp.scale(adata)

# -------------------------------
# DIMENSIONALITY REDUCTION + CLUSTERING
# -------------------------------
sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.leiden(adata)

print(adata)

# -------------------------------
# SPATIAL PLOTS
# -------------------------------
fig, ax = plt.subplots(figsize=(8, 8))
sq.pl.spatial_scatter(
    adata,
    shape=None,
    color="leiden",
    size=50,
    ax=ax
)
plt.savefig(f"figures/{args.prefix}_leiden.png")
plt.close()

sq.gr.spatial_neighbors(adata, radius=3.0)

fig, ax = plt.subplots(figsize=(8, 8))
sq.pl.spatial_scatter(
    adata,
    color="leiden",
    connectivity_key="spatial_connectivities",
    size=50,
    ax=ax
)
plt.savefig(f"figures/{args.prefix}_neighbors.png")
plt.close()

# -------------------------------
# ONE QC FIGURE (PER CLUSTER, AFTER LEIDEN)
# -------------------------------
fig = sc.pl.violin(
    adata,
    keys=["total_counts", "n_genes_by_counts", "pct_counts_mt"],
    groupby="leiden",
    multi_panel=True,
    stripplot=False,
    show=False
)
plt.savefig(f"figures/{args.prefix}_qc_per_cluster.png")
plt.close()

# -------------------------------
# SAVE
# -------------------------------
adata.write(args.output)
print(f"Processed AnnData saved to {args.output}")
