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

# --- Standard preprocessing ---
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)

# --- Highly Variable Genes ---
sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000)

# --- Scale the data ---
sc.pp.scale(adata)

# --- PCA / neighbors / UMAP / Leiden ---
sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.leiden(adata)
print(adata)

# --- Spatial scatter plot ---
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

# --- Compute spatial neighbors ---
sq.gr.spatial_neighbors(adata, radius=3.0)

# --- Spatial scatter with edges ---
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

# Save processed AnnData
adata.write(args.output)
print(f"Processed AnnData saved to {args.output}")
