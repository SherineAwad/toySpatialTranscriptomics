from numpy.random import default_rng
import scanpy as sc
from anndata import AnnData
import numpy as np
import argparse

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Generate toy AnnData with spatial coordinates")
parser.add_argument("--output", type=str, default="toy_spatial.h5ad", help="Path to save the AnnData (.h5ad) file")
args = parser.parse_args()

# Generate toy data
rng = default_rng(42)
counts = rng.integers(0, 15, size=(10, 100))  # feature matrix (10 cells x 100 genes)
coordinates = rng.uniform(0, 10, size=(10, 2))  # spatial coordinates
image = rng.uniform(0, 1, size=(10, 10, 3))  # image (optional)

# Create AnnData object with spatial coordinates
adata = AnnData(counts, obsm={"spatial": coordinates})

# Optional: add variable names and cell names
adata.var_names = [f"Gene{i}" for i in range(counts.shape[1])]
adata.obs_names = [f"Cell{i}" for i in range(counts.shape[0])]

# Save to .h5ad file
adata.write(args.output)
print(f"Saved AnnData to {args.output}")
