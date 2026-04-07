import scanpy as sc
import squidpy as sq
import argparse
import matplotlib.pyplot as plt
import numpy as np
import scipy.sparse as sp

# Arguments
parser = argparse.ArgumentParser(description="Filter AnnData with QC plots and optional MT genes")
parser.add_argument("--input", type=str, required=True, help="Path to input AnnData (.h5ad) file")
parser.add_argument("--output", type=str, required=True, help="Path to save filtered AnnData (.h5ad) file")
parser.add_argument("--prefix", type=str, default="spatial_", help="Prefix for QC plot files")
parser.add_argument("--min_genes", type=int, default=50, help="Minimum genes per cell")
parser.add_argument("--min_cells", type=int, default=3, help="Minimum cells per gene")
parser.add_argument("--mt", type=str, default="MT", help="Prefix for mitochondrial genes; pass NULL to ignore")
args = parser.parse_args()

# Load AnnData
adata = sc.read(args.input)

# FIX: make gene names unique (required for safe slicing)
adata.var_names_make_unique()

# Compute QC metrics
def compute_qc(adata, mt_prefix):
    adata.obs['n_counts'] = adata.X.sum(axis=1)
    adata.obs['n_genes'] = (adata.X > 0).sum(axis=1)

    # Only calculate percent_mt if mt_prefix is not NULL
    if mt_prefix.upper() != "NULL":
        mt_cols = adata.var_names[adata.var_names.str.startswith(mt_prefix)]

        if len(mt_cols) > 0:
            if sp.issparse(adata.X):
                mt_counts = adata[:, mt_cols].X.toarray().sum(axis=1)
            else:
                mt_counts = adata[:, mt_cols].X.sum(axis=1)

            adata.obs['percent_mt'] = mt_counts / adata.obs['n_counts']
        else:
            adata.obs['percent_mt'] = 0.0
            print(f"No genes found with prefix '{mt_prefix}'. percent_mt set to 0 for QC plots.")

        return ['n_genes', 'n_counts', 'percent_mt']

    else:
        print("Mitochondrial gene QC ignored.")
        return ['n_genes', 'n_counts']


# Prefilter QC
qc_metrics = compute_qc(adata, args.mt)

sc.pl.violin(
    adata,
    qc_metrics,
    jitter=0.4,
    multi_panel=True,
    save=f"{args.prefix}pre_filter.png"
)

# Filter
sc.pp.filter_cells(adata, min_genes=args.min_genes)
sc.pp.filter_genes(adata, min_cells=args.min_cells)

# Postfilter QC
qc_metrics = compute_qc(adata, args.mt)

sc.pl.violin(
    adata,
    qc_metrics,
    jitter=0.4,
    multi_panel=True,
    save=f"{args.prefix}post_filter.png"
)

# Save filtered AnnData
adata.write(args.output)
print(f"Filtered AnnData saved to {args.output}")
