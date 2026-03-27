#!/usr/bin/env python

import argparse
import scanpy as sc
import squidpy as sq
import pandas as pd

# -------------------------
# Arguments
# -------------------------
parser = argparse.ArgumentParser(description="Annotate AnnData and plot using cell types")
parser.add_argument("--input", type=str, required=True, help="Input AnnData (.h5ad)")
parser.add_argument("--output", type=str, required=True, help="Output annotated AnnData (.h5ad)")
parser.add_argument("--annotations", type=str, required=True, help="Annotation file (leiden,cell_type)")
parser.add_argument("--prefix", type=str, default="annotated", help="Prefix for plots")
args = parser.parse_args()

# -------------------------
# Load AnnData
# -------------------------
adata = sc.read(args.input)
print(adata)

# -------------------------
# Load annotations
# -------------------------
try:
    ann = pd.read_csv(args.annotations, sep=None, engine="python")
    if not {"leiden", "cell_type"}.issubset(ann.columns):
        raise ValueError("Annotation file must contain columns: leiden, cell_type")
except Exception as e:
    raise RuntimeError(f"Failed to read annotations file: {e}")

# -------------------------
# Apply annotations
# -------------------------
adata.obs["leiden"] = adata.obs["leiden"].astype(str)
ann["leiden"] = ann["leiden"].astype(str)

mapping = dict(zip(ann["leiden"], ann["cell_type"]))
adata.obs["cell_type"] = adata.obs["leiden"].map(mapping)

# Report missing clusters
missing = sorted(set(adata.obs["leiden"]) - set(mapping.keys()))
if len(missing) > 0:
    print(f"Warning: clusters {missing} have no annotation (set to NaN)")

# -------------------------
# Decide coloring
# -------------------------
if "cell_type" in adata.obs and not adata.obs["cell_type"].isna().all():
    color_by = "cell_type"
    print("Using cell_type for plotting")
else:
    color_by = "leiden"
    print("cell_type missing/empty → falling back to leiden")

# -------------------------
# Save annotated AnnData
# -------------------------
adata.write(args.output)
print(f"Annotated AnnData saved to {args.output}")

# -------------------------
# Plots
# -------------------------

# UMAP
if "X_umap" in adata.obsm:
    sc.pl.umap(
        adata,
        color=color_by,
        size=50,
        save=f"{args.prefix}_umap.png"
    )

# Spatial scatter
if "spatial" in adata.obsm:
    sq.pl.spatial_scatter(
        adata,
        color=color_by,
        shape=None,
        size=50,
        save=f"{args.prefix}_spatial.png"
    )

print("Annotation plots completed with:", color_by)
