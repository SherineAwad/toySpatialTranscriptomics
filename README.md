# Spatial Transcriptomics Pipeline

## Data
Visium breast cancer spatial transcriptomics: 3,798 spots, 36,601 genes. Each spot contains spatial coordinates (x,y) and gene expression counts. The tissue section includes tumor regions, stroma, immune infiltrates, and vascular elements.

## Results

![](figures/visium_BC_leiden.png?v=2)

**Leiden clustering** - Unsupervised clustering identifies distinct transcriptional regions across the tissue. Each color represents a cluster of spots with similar expression patterns, revealing heterogeneous tissue architecture and spatial compartments.

![](figures/figures/visium_BC_neighbors.png?v=2)

**Spatial neighbor graph** - Connectivity network showing physically adjacent spots. This defines spatial relationships used for neighborhood analysis and autocorrelation statistics.

![](figures/figures/visium_BC_umap_celltype.png?v=2)

**UMAP of annotated cell types** - Dimensionality reduction showing how cell types cluster in expression space. Well-separated groups indicate successful marker-based annotation with distinct transcriptional identities.

![](figures/figures/visium_BC_2_spatial_celltype.png?v=2)

**Spatial cell type distribution** - Maps each annotated cell type to its physical location. Reveals tumor heterogeneity across regions, immune infiltration patterns at invasive margins, and stromal tissue organization.

![](figures/visium_BC_spatial_celltype.png?v=2)

**Spatial cell types (alternate view)** - Same annotations with different visualization parameters, confirming spatial patterns are robust.

![](figures/visium_BC_3_spatial_graph.png?v=2)

**Spatial graph with cell types** - Neighbor network overlaid with cell type annotations, showing connectivity patterns between different cell populations.

![](figures/visium_BC_5_gene_3_SERF2.png?v=2)

**SERF2 spatial expression** - Distribution of SERF2 gene expression across the tissue. Highlights regions where this gene is active.

![](figures/visium_BC_umap_celltype.png?v=2)

**UMAP cell types (reference)** - Reference UMAP showing the relationship between all annotated cell types.

![](figures/visium_BC_5_gene_1_RPL41.png?v=2)

**RPL41 spatial expression** - Ribosomal protein RPL41 expression pattern. High expression marks regions of active protein synthesis and metabolic activity.

![](figures/visium_BC_5_gene_4_RPL13A.png?v=2)

**RPL13A spatial expression** - Ribosomal protein RPL13A distribution. Shows spatial heterogeneity in ribosomal gene expression across the tissue.

![](figures/visium_BC_5_gene_2_RPS27.png?v=2)

**RPS27 spatial expression** - Ribosomal protein RPS27 pattern. Co-localizes with other ribosomal genes indicating regional metabolic activity.

![](figures/visium_BC_4_nhood_enrichment.png?v=2)

**Neighborhood enrichment heatmap** - Cell type pair co-occurrence analysis. Red (positive values) indicates pairs found together more often than expected by chance, suggesting spatial attraction and potential functional interactions. Blue (negative values) indicates spatial exclusion or mutual avoidance.

![](figures/visium_BC_1_celltype_composition.png?v=2)

**Cell type composition** - Bar plot quantifying frequency of each annotated cell type across all spatial spots. Reveals dominant populations (tumor, immune, stromal) and rare cell types in the breast cancer microenvironment.

## References

- **Scanpy** – single-cell analysis in Python: [Wolf et al., Genome Biology 2018](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1382-0)
- **Squidpy** – spatial single-cell analysis in Python: [Palla et al., Nature Methods 2022](https://www.nature.com/articles/s41592-021-01322-8)
