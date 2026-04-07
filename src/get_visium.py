import scanpy as sc

adata = sc.datasets.visium_sge(sample_id="V1_Breast_Cancer_Block_A_Section_1")
adata.write("visium_bc.h5ad")
