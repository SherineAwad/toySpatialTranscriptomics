from numpy.random import default_rng
import scanpy as sc
from anndata import AnnData
import numpy as np
import argparse

# --- Parse command-line arguments ---
parser = argparse.ArgumentParser(description="Generate structured toy AnnData with spatial signal and multiple clusters/tissues")
parser.add_argument("--output", type=str, default="toy_spatial.h5ad", help="Path to save the AnnData (.h5ad) file")
args = parser.parse_args()

# --- Setup ---
rng = default_rng(42)
n_cells = 100
n_genes = 100
cells_per_subcluster = n_cells // 4  # 4 clusters

# --- Generate structured counts ---
counts = rng.poisson(lam=2, size=(n_cells, n_genes))

# Inject cluster-specific signals
# Cluster 0: liver subcluster 1
counts[0*cells_per_subcluster:1*cells_per_subcluster, 0:10] += rng.poisson(lam=5, size=(cells_per_subcluster, 10))
# Cluster 1: liver subcluster 2
counts[1*cells_per_subcluster:2*cells_per_subcluster, 10:20] += rng.poisson(lam=5, size=(cells_per_subcluster, 10))
# Cluster 2: kidney subcluster 1
counts[2*cells_per_subcluster:3*cells_per_subcluster, 20:30] += rng.poisson(lam=5, size=(cells_per_subcluster, 10))
# Cluster 3: kidney subcluster 2
counts[3*cells_per_subcluster:4*cells_per_subcluster, 30:40] += rng.poisson(lam=5, size=(cells_per_subcluster, 10))

# --- Spatial coordinates (separate regions per tissue) ---
coords_cluster0 = rng.normal(loc=2, scale=0.5, size=(cells_per_subcluster, 2))
coords_cluster1 = rng.normal(loc=3, scale=0.5, size=(cells_per_subcluster, 2))
coords_cluster2 = rng.normal(loc=8, scale=0.5, size=(cells_per_subcluster, 2))
coords_cluster3 = rng.normal(loc=9, scale=0.5, size=(cells_per_subcluster, 2))
coordinates = np.vstack([coords_cluster0, coords_cluster1, coords_cluster2, coords_cluster3])

# --- Create AnnData ---
adata = AnnData(counts, obsm={"spatial": coordinates})

# Names
adata.var_names = [f"Gene{i}" for i in range(n_genes)]
adata.obs_names = [f"Cell{i}" for i in range(n_cells)]

# --- Assign clusters (for reference / annotation) ---
adata.obs["leiden"] = [0]*cells_per_subcluster + [1]*cells_per_subcluster + [2]*cells_per_subcluster + [3]*cells_per_subcluster

# --- Assign slides (aligned with tissues) ---
adata.obs["slide"] = ["slide1"]*(2*cells_per_subcluster) + ["slide2"]*(2*cells_per_subcluster)

# --- Map to tissue ---
tissue_map = {"slide1": "liver", "slide2": "kidney"}
adata.obs["tissue"] = [tissue_map[s] for s in adata.obs["slide"]]

# --- Print ---
print(adata)
print("\nMetadata (obs):")
print(adata.obs["leiden"].value_counts())
print(adata.obs["tissue"].value_counts())

# --- Save ---
adata.write(args.output)
print(f"\nSaved AnnData to {args.output}")
