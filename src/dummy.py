from numpy.random import default_rng
import scanpy as sc
from anndata import AnnData
import argparse

# --- Parse command-line arguments ---
parser = argparse.ArgumentParser(description="Generate toy AnnData with spatial coordinates and slide/tissue info")
parser.add_argument("--output", type=str, default="toy_spatial.h5ad", help="Path to save the AnnData (.h5ad) file")
args = parser.parse_args()

# --- Generate toy data ---
rng = default_rng(42)
n_cells = 10
n_genes = 100

# Technical slides (dummy)
slides_available = ["slide1", "slide2"]  

# Mapping slides to dummy tissue (for testing only)
tissue_map = {"slide1": "liver", "slide2": "kidney"}  

# Gene counts
counts = rng.integers(0, 15, size=(n_cells, n_genes))  # 10 cells × 100 genes

# Spatial coordinates
coordinates = rng.uniform(0, 10, size=(n_cells, 2))  # x, y coordinates

# --- Create AnnData object ---
adata = AnnData(counts, obsm={"spatial": coordinates})

# --- Add gene and cell names ---
adata.var_names = [f"Gene{i}" for i in range(n_genes)]
adata.obs_names = [f"Cell{i}" for i in range(n_cells)]

# --- Randomized slide assignment (technical origin) ---
adata.obs["slide"] = rng.choice(slides_available, size=n_cells)

# --- Optional tissue assignment based on slide (for pipeline testing) ---
adata.obs["tissue"] = [tissue_map[s] for s in adata.obs["slide"]]

# --- Print summary ---
print(adata)
print("\nMetadata (obs):")
print(adata.obs)

# --- Save to .h5ad file ---
adata.write(args.output)
print(f"\nSaved AnnData to {args.output}")
