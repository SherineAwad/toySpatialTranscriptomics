#!/usr/bin/env python

import argparse
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

# -------------------------
# Args
# -------------------------
parser = argparse.ArgumentParser(description="Marker-based cell type annotation")
parser.add_argument("--input", required=True, help="Input AnnData (.h5ad)")
parser.add_argument("--output", required=True, help="Output AnnData (.h5ad)")
parser.add_argument("--markers", required=True, help="CSV: cell_type,markers")
parser.add_argument("--prefix", default="annotated", help="Output prefix for plots/files")
parser.add_argument("--min_genes", type=int, default=2)
parser.add_argument("--unknown_thresh", type=float, default=0.1)
args = parser.parse_args()

# -------------------------
# Load data
# -------------------------
adata = sc.read(args.input)
print(adata)

# Fix: Make variable names unique and ensure strings
adata.var_names_make_unique()
adata.var_names = adata.var_names.astype(str)

# -------------------------
# Load markers
# -------------------------
marker_df = pd.read_csv(args.markers)

if "cell_type" not in marker_df.columns or "markers" not in marker_df.columns:
    raise ValueError("Marker file must contain: cell_type, markers")

marker_dict = {
    row["cell_type"]: [g.strip() for g in row["markers"].split(";")]
    for _, row in marker_df.iterrows()
}

print(f"Loaded {len(marker_dict)} cell types")

# -------------------------
# Score genes
# -------------------------
score_df = pd.DataFrame(index=adata.obs_names)

for cell_type, genes in marker_dict.items():
    # Filter to genes that exist in adata
    genes = [g for g in genes if g in adata.var_names]

    if len(genes) < args.min_genes:
        print(f"Skipping {cell_type}: only {len(genes)} genes found")
        continue

    sc.tl.score_genes(
        adata,
        gene_list=genes,
        score_name=f"score_{cell_type}"
    )

    score_df[cell_type] = adata.obs[f"score_{cell_type}"].values

# -------------------------
# Normalize scores
# -------------------------
score_df = (score_df - score_df.mean()) / (score_df.std() + 1e-9)

# -------------------------
# Assign labels
# -------------------------
adata.obs["cell_type"] = score_df.idxmax(axis=1)
adata.obs["cell_type_score"] = score_df.max(axis=1)

adata.obs.loc[
    adata.obs["cell_type_score"] < args.unknown_thresh,
    "cell_type"
] = "Unknown"

# -------------------------
# Save AnnData
# -------------------------
adata.write(args.output)
print(f"Saved AnnData → {args.output}")

# -------------------------
# PLOTS - saved to figures/
# -------------------------

# UMAP
if "X_umap" in adata.obsm:
    fig, ax = plt.subplots(figsize=(10, 8))
    sc.pl.umap(
        adata,
        color="cell_type",
        size=40,
        ax=ax,
        show=False
    )
    plt.tight_layout()
    plt.savefig(f"figures/{args.prefix}_umap_celltype.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("✓ UMAP plot saved to figures/")

# Spatial
if "spatial" in adata.obsm:
    import squidpy as sq
    fig, ax = plt.subplots(figsize=(10, 10))
    sq.pl.spatial_scatter(
        adata,
        color="cell_type",
        size=40,
        shape=None,
        ax=ax
    )
    plt.tight_layout()
    plt.savefig(f"figures/{args.prefix}_spatial_celltype.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("✓ Spatial plot saved to figures/")

print("DONE ✔")
