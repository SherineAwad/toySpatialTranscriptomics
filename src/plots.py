#!/usr/bin/env python

import argparse
import scanpy as sc
import squidpy as sq
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# -------------------------
# Args
# -------------------------
parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True)
parser.add_argument("--prefix", default="fig")
args = parser.parse_args()

prefix = args.prefix

# -------------------------
# Load data
# -------------------------
adata = sc.read(args.input)

if "cell_type" not in adata.obs:
    raise ValueError("cell_type missing")

# Ensure spatial coordinates exist
if "spatial" not in adata.obsm:
    if "X_spatial" in adata.obsm:
        adata.obsm["spatial"] = adata.obsm["X_spatial"]
    else:
        raise ValueError("No spatial coordinates found")

# =========================================================
# 1. Composition
# =========================================================
plt.figure(figsize=(10, 6))
adata.obs["cell_type"].value_counts().plot(kind="bar")
plt.title("Cell Type Composition")
plt.tight_layout()
plt.savefig(f"figures/{prefix}_celltype_composition.png", dpi=150, bbox_inches="tight")
plt.close()
print("✓ Figure 1 saved")

# =========================================================
# 2. Spatial graph
# =========================================================
if "spatial_connectivities" not in adata.obsp:
    sq.gr.spatial_neighbors(adata, coord_type="generic", n_neighbors=6)

fig, ax = plt.subplots(figsize=(10, 10))
sc.pl.embedding(adata, basis="spatial", color="cell_type", size=30, ax=ax, show=False)

if "spatial_connectivities" in adata.obsp:
    from scipy.sparse import find
    coords = adata.obsm["spatial"]
    edges = find(adata.obsp["spatial_connectivities"])
    step = max(1, len(edges[0]) // 2000)
    for idx in range(0, len(edges[0]), step):
        i, j = edges[0][idx], edges[1][idx]
        if i < j:
            ax.plot([coords[i,0], coords[j,0]], [coords[i,1], coords[j,1]], 
                   'gray', alpha=0.2, linewidth=0.3)

plt.title("Spatial Graph")
plt.tight_layout()
plt.savefig(f"figures/{prefix}_spatial_graph.png", dpi=150, bbox_inches="tight")
plt.close()
print("✓ Figure 2 saved")

# =========================================================
# 3. Neighborhood enrichment
# =========================================================
sq.gr.nhood_enrichment(adata, cluster_key="cell_type")

fig, ax = plt.subplots(figsize=(12, 10))
sq.pl.nhood_enrichment(adata, cluster_key="cell_type", cmap="viridis", ax=ax)
plt.title("Neighborhood Enrichment")
plt.tight_layout()
plt.savefig(f"figures/{prefix}_nhood_enrichment.png", dpi=150, bbox_inches="tight")
plt.close()
print("✓ Figure 3 saved")

