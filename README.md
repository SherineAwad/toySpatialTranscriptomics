# Toy Spatial Transcriptomics Pipeline

This repository contains a minimal set of scripts for a toy spatial transcriptomics workflow. The scripts are organized in the `src/` directory.


## Create a dummy spatial scRNA dataset using `dummy.py` script 

- **Number of cells (observations):** 10  
- **Number of genes (variables):** 100  
- **Data type:** raw counts (integers from 0 to 15)  

## Cell Metadata (`adata.obs`)

- **Cell IDs:** `Cell0`, `Cell1`, … `Cell9`  
- No additional annotations  

## Gene Metadata (`adata.var`)

- **Gene IDs:** `Gene0`, `Gene1`, … `Gene99`  
- No extra annotations  

## Expression Matrix (`adata.X`)

- Shape: `(10, 100)`  
- Each row corresponds to a cell, each column to a gene  
- Values: integer counts between 0 and 15  

## Spatial Coordinates (`adata.obsm["spatial"]`)

- Shape: `(10, 2)`  
- Floating-point coordinates representing a 2D spatial layout  
- Each row corresponds to a cell’s (x, y) location  


## Visualization: Filtering  QC

To remove low-quality cells and features, we filtered out cells with too few detected genes and genes expressed in too few cells.  
A QC violin plot was generated to visualize these metrics and check that the filtering worked as expected.


### Before filter 

![Pre-filter QC Violin Plot](figures/violintoy_spatialpre_filter.png?v=3)

### Post filter 

![Post-filter QC Violin Plot](figures/violintoy_spatialpost_filter.png?v=3)


## Analysis

In this step, we process the filtered data to prepare it for downstream analysis:


1. **Standard preprocessing**:  
   - The total counts for each cell are normalized to make cells comparable.  
   - Logarithmic transformation (`log1p`) is applied to stabilize variance across genes.  

2. **Highly Variable Genes (HVGs)**:  
   - We identify the top 2000 most variable genes using the Seurat method.  
   - Optionally, the dataset can be subsetted to include only these highly variable genes for downstream analysis.  

3.. **Scaling**: Each gene is scaled to zero mean and unit variance, ensuring that highly expressed genes do not dominate the analysis.  

4. **Dimensionality reduction and clustering**:  
   - **PCA** reduces the data to key components capturing most of the variation.  
   - **Neighbors graph** is computed to identify similar cells.  
   - **UMAP** projects cells into 2D for visualization.  
   - **Leiden clustering** groups cells into clusters based on similarity.  



### Spatial scatter plot

   - Shows the location of each cell in the tissue.  
   - Cells are colored by their **Leiden cluster**, making it easy to see which clusters occupy which regions.  
   - Helps identify spatial patterns or distinct domains of cell populations.

![Leiden clusters](figures/toy_spatial_leiden.png?v=3)  



### Spatial scatter with edges

   - Computes spatial neighbors for each cell based on a defined radius.  
   - Plots cells with **edges connecting neighboring cells** to show local spatial relationships.  
   - Cell colors still reflect Leiden clusters, while edges reveal how clusters interact spatially.  
   - Useful for examining spatial connectivity and potential interactions between cell groups.

![Spatial neighbors](figures/toy_spatial_neighbors.png?v=3)


### UMAP Plot

A 2D representation of cells based on gene expression.  
- Each dot is a cell, colored by its Leiden cluster.  
- Cells close together have similar gene expression.  
- This plot shows **expression similarity**, not tissue location.

![UMAP Plot](figures/umaptoy_spatial_umap.png?v=3)



### Top Spatially Variable Genes (Moran's I)

This plot highlights the **genes with the strongest spatial patterns** in the tissue.  

- **Moran's I** is a statistic that measures how gene expression is spatially autocorrelated — in other words, whether high or low expression tends to cluster together in space.  
- The top 5 genes with the highest Moran's I values are selected.  
- Cells are plotted in their tissue positions, colored by the expression of these top genes.  
- Purpose: Quickly visualize genes whose expression shows **strong spatial organization**, which can reveal tissue structure or spatially distinct cell populations.  

![Top Moran's I genes](figures/toy_spatial_top_moranI.png?v=3)




### Marker Genes Spatial Plots

Cells are plotted in their tissue positions and colored by the expression of selected marker genes.  
Each gene has its own plot, showing where it is highly or lowly expressed in the tissue.  

##### As an example we plot marker genes Gene92 and Gene77
![Example Marker Gene](figures/toy_spatial_gene_Gene92.png?v=3)

![Example Marker Gene](figures/toy_spatial_gene_Gene77.png?v=3)


### Now annotate 


We annotating clusters where Leiden cluster labels are replaced with biological **cell type annotations**.



![](figures/umapannotated_umap_celltype.png?v=3)
![](figures/umapannotated_umap_tissue.png?v=3)
![](figures/annotated_spatial_celltype.png?v=3)
![](figures/annotated_spatial_tissue.png?v=3)


## References

- **Scanpy** – single-cell analysis in Python:  
  [Wolf et al., Genome Biology 2018](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1382-0)  
  Python package: [https://scanpy.readthedocs.io](https://scanpy.readthedocs.io)

- **Squidpy** – spatial single-cell analysis in Python:  
  [Palla et al., Nature Methods 2022](https://www.nature.com/articles/s41592-021-01322-8)  
  Python package: [https://squidpy.readthedocs.io](https://squidpy.readthedocs.io)
