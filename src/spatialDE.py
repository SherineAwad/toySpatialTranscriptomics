#!/usr/bin/env python

import argparse
import scanpy as sc
import squidpy as sq
import pandas as pd
import matplotlib.pyplot as plt

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

sq.gr.spatial_autocorr(adata, mode="moran")

moran_df = adata.uns["moranI"]

# genes are index, scores are in column "I"
moran_results = moran_df[["I"]].copy()
moran_results = moran_results.rename(columns={"I": "moranI"})
moran_results = moran_results.sort_values("moranI", ascending=False)

top_genes = moran_results.head(args.topN).index.tolist()

moran_results.to_csv(f"{args.prefix}_spatial_genes.csv")

for i, gene in enumerate(top_genes):
    fig, ax = plt.subplots(figsize=(8, 8))
    sc.pl.embedding(
        adata,
        basis="spatial",
        color=gene,
        size=30,
        cmap="magma",
        ax=ax,
        show=False
    )
    plt.savefig(f"figures/{args.prefix}_gene_{i+1}_{gene}.png", dpi=150, bbox_inches="tight")
    plt.close()

adata.write(args.output)
