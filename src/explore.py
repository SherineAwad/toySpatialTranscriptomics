import scanpy as sc
import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True)
args = parser.parse_args()

adata = sc.read_h5ad(args.input)

print(adata)

gene_names = adata.var_names

print("\nTop 4 gene names:")
print(gene_names[:4])

mt_mask = gene_names.str.upper().str.startswith("MT-")
mt_genes = gene_names[mt_mask]

print("\nTop 4 mitochondrial genes:")
print(mt_genes[:4])
