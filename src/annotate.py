#!/usr/bin/env python

import argparse
import scanpy as sc
import squidpy as sq
import pandas as pd

# -------------------------
# Arguments
# -------------------------
parser = argparse.ArgumentParser(description="Annotate AnnData and plot using cell types/tissue")
parser.add_argument("--input", type=str, required=True, help="Input AnnData (.h5ad)")
parser.add_argument("--output", type=str, required=True, help="Output annotated AnnData (.h5ad)")
parser.add_argument("--annotations", type=str, required=True, help="Annotation file (leiden, cell_type, tissue)")
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
        raise ValueError("Annotation file must contain columns: leiden, cell_type (tissue optional)")
except Exception as e:
    raise RuntimeError(f"Failed to read annotations file: {e}")

# -------------------------
# Apply annotations
# -------------------------
adata.obs["leiden"] = adata.obs["leiden"].astype(str)
ann["leiden"] = ann["leiden"].astype(str)

# --- Map cell_type ---
mapping_cell_type = dict(zip(ann["leiden"], ann["cell_type"]))
adata.obs["cell_type"] = adata.obs["leiden"].map(mapping_cell_type)

# --- Map tissue if present ---
if "tissue" in ann.columns:
    mapping_tissue = dict(zip(ann["leiden"], ann["tissue"]))
    adata.obs["tissue"] = adata.obs["leiden"].map(mapping_tissue)
    print("Tissue annotation applied")
else:
    print("No tissue column in annotation CSV → skipping tissue annotation")

# --- Report missing clusters ---
missing_cell_type = sorted(set(adata.obs["leiden"]) - set(mapping_cell_type.keys()))
if missing_cell_type:
    print(f"Warning: clusters {missing_cell_type} have no cell_type annotation (set to NaN)")

if "tissue" in ann.columns:
    missing_tissue = sorted(set(adata.obs["leiden"]) - set(mapping_tissue.keys()))
    if missing_tissue:
        print(f"Warning: clusters {missing_tissue} have no tissue annotation (set to NaN)")

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
    if "cell_type" in adata.obs:
        sc.pl.umap(
            adata,
            color="cell_type",
            size=50,
            save=f"{args.prefix}_umap_celltype.png"
        )
    if "tissue" in adata.obs:
        sc.pl.umap(
            adata,
            color="tissue",
            size=50,
            save=f"{args.prefix}_umap_tissue.png"
        )

# Spatial scatter
if "spatial" in adata.obsm:
    if "cell_type" in adata.obs:
        sq.pl.spatial_scatter(
            adata,
            color="cell_type",
            shape=None,
            size=50,
            save=f"{args.prefix}_spatial_celltype.png"
        )
    if "tissue" in adata.obs:
        sq.pl.spatial_scatter(
            adata,
            color="tissue",
            shape=None,
            size=50,
            save=f"{args.prefix}_spatial_tissue.png"
        )

print("Annotation plots completed for cell_type and tissue")
