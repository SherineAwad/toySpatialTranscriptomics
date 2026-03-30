#!/usr/bin/env python

import scanpy as sc
import squidpy as sq
import argparse
import os
import pandas as pd

# -------------------------
# Arguments
# -------------------------
parser = argparse.ArgumentParser(description="Spatial differential expression by tissue")
parser.add_argument("--input", type=str, required=True, help="Preprocessed AnnData (.h5ad)")
parser.add_argument("--output", type=str, required=True, help="Output file for DE results (.csv) and AnnData (.h5ad)")
parser.add_argument("--prefix", type=str, default="spatial_de", help="Prefix for heatmap plot")
parser.add_argument("--top", type=int, default=5, help="Top N genes per tissue to include in heatmap")
args = parser.parse_args()

# -------------------------
# Load AnnData
# -------------------------
adata = sc.read(args.input)
print(adata)

# -------------------------
# Validate tissue column
# -------------------------
if "tissue" not in adata.obs:
    raise ValueError("AnnData.obs must contain a 'tissue' column for spatial DE.")

# Remove NaN tissues
adata = adata[~adata.obs["tissue"].isna()].copy()

if adata.obs["tissue"].nunique() < 2:
    raise ValueError("Need at least 2 tissue groups for DE analysis.")

print("Tissues found:", adata.obs["tissue"].unique())

# -------------------------
# Create figures directory
# -------------------------
fig_dir = "figures"
os.makedirs(fig_dir, exist_ok=True)
sc.settings.figdir = fig_dir

# -------------------------
# Run DE per tissue
# -------------------------
sc.tl.rank_genes_groups(
    adata,
    groupby="tissue",
    method="wilcoxon",
    pts=True
)

# -------------------------
# Save DE results CSV
# -------------------------
# CSV is always separate, ensure .csv extension
csv_path = args.output
if not csv_path.endswith(".csv"):
    csv_path = csv_path + ".csv"

de_results = sc.get.rank_genes_groups_df(adata, group=None)
de_results.to_csv(csv_path, index=False)
print(f"Spatial DE results saved to {csv_path}")

# -------------------------
# Save AnnData exactly as --output
# -------------------------
adata.write(args.output)
print(f"AnnData object saved exactly to {args.output}")

# -------------------------
# Prepare top genes for heatmap
# -------------------------
top_genes_all = []
for tissue in adata.obs["tissue"].unique():
    df_tissue = sc.get.rank_genes_groups_df(adata, group=tissue)
    top_genes = df_tissue["names"].head(args.top).tolist()
    top_genes_all.extend(top_genes)

# Remove duplicates
top_genes_all = list(dict.fromkeys(top_genes_all))

# -------------------------
# Heatmap
# -------------------------
sc.pl.heatmap(
    adata,
    var_names=top_genes_all,
    groupby="tissue",
    show=False,
    save=f"_{args.prefix}_top_genes.png"
)

print(f"\nHeatmap of top {args.top} genes per tissue saved in '{fig_dir}' as '{args.prefix}_top_genes.png'")
