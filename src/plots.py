import scanpy as sc
import squidpy as sq
import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Arguments
parser = argparse.ArgumentParser(description="Plot preprocessed AnnData with optional marker genes")
parser.add_argument("--input", type=str, required=True, help="Path to input AnnData (.h5ad)")
parser.add_argument("--prefix", type=str, default="spatial_", help="Prefix for output plots")
parser.add_argument("--markers", type=str, default="markers.txt", help="Path to marker genes file (one gene per line)")
args = parser.parse_args()

# Load AnnData
adata = sc.read(args.input)
print(adata)

# --- UMAP plot ---
if "X_umap" in adata.obsm:
    sc.pl.umap(
        adata,
        color="leiden" if "leiden" in adata.obs else None,
        size=500,  # smaller dot size
        save=f"{args.prefix}_umap.png"
    )


# --- Moran's I for spatial autocorrelation ---
try:
    sq.gr.spatial_autocorr(adata, mode="moran")
    # Safe way to get Moran's I values
    if "moranI" in adata.var.columns:
        moran_series = adata.var['moranI']
    elif "moranI" in adata.uns:
        moran_series = pd.Series(adata.uns["moranI"]["I"], index=adata.uns["moranI"].index)
    else:
        moran_series = None

    if moran_series is not None and not moran_series.empty:
        # Top 5 genes for quick overview
        top_genes = list(moran_series.sort_values(ascending=False).head(5).index)
        print("Top Moran's I genes:", top_genes)
        sq.pl.spatial_scatter(
            adata,
            color=top_genes,
            shape=None,
            size=50,
            save=f"{args.prefix}_top_moranI.png"
        )

        # Gene-by-gene plots for marker genes
        try:
            markers = []
            with open(args.markers, "r") as f:
                markers = [line.strip() for line in f if line.strip() in adata.var_names]
            if markers:
                for g in markers:
                    sq.pl.spatial_scatter(
                        adata,
                        color=g,
                        shape=None,
                        size=50,
                        save=f"{args.prefix}_gene_{g}.png"
                    )
            else:
                print("No markers found in AnnData.var_names; skipping marker plots.")
        except FileNotFoundError:
            print(f"Markers file {args.markers} not found; skipping gene-by-gene plots.")

    else:
        print("Moran's I not available; skipping top genes and marker plots.")

except Exception as e:
    print("Spatial autocorrelation (Moran's I) could not be computed:", e)

print("All plots saved with prefix:", args.prefix)
