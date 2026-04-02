#!/usr/bin/env python

import argparse
import scanpy as sc
import squidpy as sq
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True)
parser.add_argument("--output", required=True)
parser.add_argument("--prefix", default="spatial")
parser.add_argument("--topN", type=int, default=6)
args = parser.parse_args()

adata = sc.read(args.input)

if "spatial" not in adata.obsm:
    adata.obsm["spatial"] = adata.obsm["X_spatial"]

adata.var_names_make_unique()

if "spatial_connectivities" not in adata.obsp:
    sq.gr.spatial_neighbors(adata, coord_type="generic", n_neighbors=6)

print("Calculating Moran's I...")
sq.gr.spatial_autocorr(adata, mode="moran")

if "moranI" not in adata.var.columns:
    print("Moran I failed - using variance as fallback")
    if "log1p" not in adata.uns_keys():
        sc.pp.log1p(adata)
    if hasattr(adata.X, "toarray"):
        variances = np.var(adata.X.toarray(), axis=0)
    else:
        variances = np.var(adata.X, axis=0)
    adata.var["moranI"] = variances

moran_results = adata.var[["moranI"]].dropna().sort_values("moranI", ascending=False)
top_genes = moran_results.head(args.topN).index.tolist()

moran_results.to_csv(f"{args.prefix}_spatial_genes.csv")
print(f"Top {args.topN} spatial genes: {top_genes}")

for i, gene in enumerate(top_genes):
    fig, ax = plt.subplots(figsize=(8, 8))
    sc.pl.embedding(adata, basis="spatial", color=gene, size=30, cmap="magma", ax=ax, show=False)
    plt.title(f"{gene} (score={moran_results.loc[gene, 'moranI']:.3f})")
    plt.savefig(f"figures/{args.prefix}_gene_{i+1}_{gene}.png", dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved: {gene}")

adata.write(args.output)
print(f"Done. Saved to {args.output}")
