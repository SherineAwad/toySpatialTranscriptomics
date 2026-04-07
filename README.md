# Spatial Transcriptomics Pipeline

## Data
Visium breast cancer spatial transcriptomics: 3,798 spots, 36,601 genes. Each spot contains spatial coordinates (x,y) and gene expression counts. The tissue section includes tumor regions, stroma, immune infiltrates, and vascular elements.


#### Before Filtering  
![]( figures/violinvisium_bcpre_filter.png?v=3)

#### After Filtering (We almost did not filter as the data is already filtered) 
![](figures/violinvisium_bcpost_filter.png?v=3)

## Results

**Leiden clustering (spatial projection)** – Spots are colored by Leiden clusters derived from transcriptional similarity (PCA → neighbors → UMAP → Leiden). The plot shows how these expression-based clusters are distributed across tissue space, highlighting spatial organization of transcriptionally similar regions.
![](figures/visium_bc_leiden.png?v=3)

**Spatial neighbor graph** – Spatial scatter plot with connectivity edges showing physically adjacent spots within a fixed radius. This graph defines local spatial relationships used for downstream spatial statistics and neighborhood analysis.
![](figures/visium_bc_neighbors.png?v=3)

#### Look at QC per cluster leiden 

![](figures/visium_bc_qc_per_cluster.png?v=2)

**UMAP of annotated cell types** - Dimensionality reduction showing how cell types cluster in expression space. Well-separated groups indicate successful marker-based annotation with distinct transcriptional identities.

![](figures/visium_bc_umap_celltype.png?v=2) 

**Spatial cell type map** – Spatial scatter plot of spots colored by annotated cell types, showing how different cell populations are distributed across the tissue architecture. This reveals spatial organization and potential regional enrichment of specific cell types.

![](figures/visium_bc_spatial_celltype.png?v=3)

### Cell Type Composition
This bar plot shows the overall abundance of each annotated cell type in the dataset. It provides a quick overview of cellular heterogeneity and highlights dominant or rare populations.
![](figures/visium_bc_celltype_composition.png?v=3)

### Spatial Graph
This visualization shows the spatial organization of spots/cells and their neighborhood connections. Nodes are colored by cell type, and edges represent spatial adjacency inferred from nearest-neighbor relationships, illustrating local tissue structure.
![](figures/visium_bc_spatial_graph.png?v=3)

### Neighborhood Enrichment
This heatmap shows pairwise enrichment or depletion of cell type neighborhoods. Positive values indicate that two cell types co-occur spatially more often than expected by chance, while negative values indicate avoidance.

## Spatially Variable Genes (SVG)

Spatially variable genes are genes whose expression levels are not randomly distributed across the tissue but instead show structured spatial patterns. These genes can highlight functional tissue regions, spatial domains, or microenvironment niches.

**Moran’s I** measures spatial autocorrelation — whether nearby spatial spots have similar expression values for a given gene. It is computed per gene using spatial neighborhood graphs (e.g., from Squidpy). Values typically range from -1 to +1, where positive values indicate spatial clustering and values near zero indicate weak or no spatial structure.

In this analysis, genes are ranked by their Moran’s I score (column **"I"** from Squidpy output), and the top genes represent the strongest spatially structured expression patterns in the tissue.

- **High Moran’s I (> 0.5)**: Gene shows strong spatial clustering in specific regions (e.g., tissue compartments, tumor core, or stromal zones)
- **Moderate Moran’s I (0.2–0.5)**: Gene shows moderate spatial structure or gradients across tissue
- **Low Moran’s I (< 0.2)**: Gene shows weak or no spatial organization and is more uniformly or randomly distributed


![](figures/visium_bc_CRISP3.png?v=2)
![](figures/visium_bc_CPB1.png?v=2)
![](figures/visium_bc_IGKC.png?v=2)
![](figures/visium_bc_IGLC2.png?v=2)
![](figures/visium_bc_IGHG3.png?v=2)
![](figures/visium_bc_MALAT1.png?v=2)

[View Spatially Variable Genes (Moran’s I)](https://docs.google.com/spreadsheets/d/18yUKNfJoli7Hmwz5dQvJJVYpO5k_upeTbhsZso9qO24/edit?usp=sharing) 

## References

- **Scanpy** – single-cell analysis in Python: [Wolf et al., Genome Biology 2018](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1382-0)
- **Squidpy** – spatial single-cell analysis in Python: [Palla et al., Nature Methods 2022](https://www.nature.com/articles/s41592-021-01322-8)
