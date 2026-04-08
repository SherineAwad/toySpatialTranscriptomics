# Spatial Transcriptomics Pipeline

## Data
Visium breast cancer spatial transcriptomics: 3,798 spots, 36,601 genes. Each spot contains spatial coordinates (x,y) and gene expression counts. The tissue section includes tumor regions, stroma, immune infiltrates, and vascular elements.

## Preprocessing 

### Before Filtering  
![]( figures/violinvisium_bcpre_filter.png?v=3)

### After Filtering (We almost did not filter as the data is already filtered) 
![](figures/violinvisium_bcpost_filter.png?v=3)

## Analysis  

**Leiden clustering (spatial projection)** – Spots are colored by Leiden clusters derived from transcriptional similarity (PCA → neighbors → UMAP → Leiden). The plot shows how these expression-based clusters are distributed across tissue space, highlighting spatial organization of transcriptionally similar regions.
![](figures/visium_bc_leiden.png?v=3)

**Spatial neighbor graph** – This plot visualizes how spatial spots are connected to one another based on a defined distance radius.

**What each dot represents:** A spatial capture spot, colored by its Leiden cluster assignment.

**What the connections show:** Lines between spots indicate that two spots are within the specified radius (here, 3.0 distance units) and are therefore considered spatial neighbors.

**What this plot tells you:**

- Whether the chosen radius appropriately captures local spatial relationships (not too sparse, not too dense)
- Whether spots from the same Leiden cluster form contiguous connected regions or are broken into isolated pieces
- How different clusters are spatially arranged relative to each other

![](figures/visium_bc_neighbors.png?v=4)

### Look at QC per cluster leiden 

![](figures/visium_bc_qc_per_cluster.png?v=2)

**UMAP of annotated cell types** - Dimensionality reduction showing how cell types cluster in expression space. Well-separated groups indicate successful marker-based annotation with distinct transcriptional identities.

![](figures/visium_bc_umap_celltype.png?v=2) 

**Spatial cell type map**  - This plot shows the spatial arrangement of cell types across the tissue section. Each dot represents a spatial capture spot, colored by its assigned cell type.

The position of each dot reflects its original tissue coordinates, allowing direct visualization of where each cell population resides and how different cell types are distributed relative to one another.

![](figures/visium_bc_spatial_celltype.png?v=3)

### Cell Type Composition
This bar plot shows the overall abundance of each annotated cell type in the dataset. It provides a quick overview of cellular heterogeneity and highlights dominant or rare populations.
![](figures/visium_bc_celltype_composition.png?v=3)

## Spatial Graph with Neighbor Connections

This plot shows spatial spots colored by cell type, with gray lines connecting neighboring spots.

The lines define which spots are considered neighbors based on the chosen distance or number of nearest neighbors. 

Use this plot to quickly verify that no spot is isolated and to see which cell types sit adjacent to each other. The lines are dense and not meant to be traced individually — this is primarily a quality control visualization.

![](figures/visium_bc_spatial_graph.png?v=3)

### Neighborhood Enrichment

This heatmap reveals which cell types tend to be spatial neighbors and which tend to avoid each other.

**How it works:** For each pair of cell types (row vs column), the method counts how often cells of the column type appear as direct spatial neighbors of cells of the row type. This is compared to a random distribution to calculate an enrichment score.

**Reading the heatmap:**

- **Positive scores (lighter colors)** = the two cell types are found next to each other more often than expected by chance → spatial attraction
- **Scores near zero (dark purple)** = the two cell types are neighbors as often as expected by chance → no spatial preference
- **Negative scores (intermediate colors)** = the two cell types are found next to each other less often than expected → spatial avoidance

**The diagonal (same cell type vs itself):** This shows whether a cell type tends to cluster with its own kind. High positive values mean cells of that type form dense same-type clusters. Low or negative values mean they are scattered or intermingled with other types.

**Why same cell types are not necessarily close neighbors:** In many tissues, cells of the same type are dispersed rather than clustered. For example, immune cells may be spaced apart with other cell types sitting between them. Two B cells can exist in the same tissue region without being direct spatial neighbors. Therefore, a low diagonal value does not mean the cell type is absent — it means that when you find one cell, the cells immediately surrounding it are often of different types.

Use this plot to identify which cell populations co-localize in tissue architecture and which occupy distinct spatial niches.

![](figures/visium_bc_nhood_enrichment.png?v=2) 


## Gene expression spatial plots for few genes

For each selected marker gene, the plot below shows its spatial expression pattern across the tissue section. 

- **Each dot** represents a spatial capture spot
- **Dot color** indicates expression level (darker = lower expression, lighter = higher expression)
- **Spatial position** reflects the original tissue coordinates

**How to interpret each gene's pattern:**

- If a gene shows **bright/warm colors concentrated in a specific region**, that gene is spatially restricted to that area
- If a gene shows **uniform or scattered expression** across the tissue, that gene is broadly expressed or marks dispersed cell populations
- If two different marker genes show **identical spatial patterns**, they may mark the same cell type or co-localized biological process
- If a gene shows **little to no expression** (dark throughout), it is not active in this tissue or sample

These plots help identify where specific cell types or biological processes are localized within the tissue.

<img src="figures/visium_bc_spatial_MT-CYB.png" width="33%"><img src="figures/visium_bc_spatial_MT-ATP6.png" width="33%"><img src="figures/visium_bc_spatial_MT-CO2.png" width="33%">
<img src="figures/visium_bc_spatial_MT-CO1.png" width="33%"><img src="figures/visium_bc_spatial_LDHA.png" width="33%"><img src="figures/visium_bc_spatial_SLC2A1.png" width="33%">
<img src="figures/visium_bc_spatial_VEGFA.png" width="33%"><img src="figures/visium_bc_spatial_HIF1A.png" width="33%"><img src="figures/visium_bc_spatial_FCGR3A.png" width="33%">
<img src="figures/visium_bc_spatial_CSF1R.png" width="33%"><img src="figures/visium_bc_spatial_CD68.png" width="33%"><img src="figures/visium_bc_spatial_LYZ.png" width="33%">
<img src="figures/visium_bc_spatial_CD79A.png" width="33%"><img src="figures/visium_bc_spatial_MS4A1.png" width="33%"><img src="figures/visium_bc_spatial_TRAC.png" width="33%">

## Spatially Variable Genes (SVG)

Spatially variable genes are genes whose expression levels are not randomly distributed across the tissue but instead show structured spatial patterns. These genes can highlight functional tissue regions, spatial domains, or microenvironment niches.

### Spatial Autocorrelation and Moran's I

### What is spatial autocorrelation?

Spatial autocorrelation measures whether nearby locations in space have similar **gene expression values**. For a single gene, each spatial spot has a value — the level of expression of that gene (e.g., transcript count). The analysis asks: *"If a gene is highly expressed in one spot, are neighboring spots also likely to have high expression values for that same gene?"*

There are three possible outcomes:

- **Positive spatial autocorrelation** – Nearby spots have similar gene expression values. The gene forms clusters, patches, or gradients across the tissue.
- **Negative spatial autocorrelation** – Nearby spots have dissimilar gene expression values (like a checkerboard pattern). This is rare in biological tissues.
- **Zero spatial autocorrelation** – Expression value at one spot tells you nothing about its neighbors. The pattern is random.

### What is Moran's I?

Moran's I is a statistical score that quantifies spatial autocorrelation for each gene individually. It typically ranges from -1 to +1:

- **I > 0** → Positive autocorrelation (clustered pattern of expression values)
- **I = 0** → Random pattern of expression values
- **I < 0** → Negative autocorrelation (dispersed or alternating pattern of expression values)

The higher the Moran's I, the stronger the spatial clustering of that gene's expression values across the tissue.

### What this means biologically

**High Moran's I (genes with spatially clustered expression values):**

These genes mark biologically meaningful spatial structures. For example:
- A gene whose expression values are high only in a specific anatomical region (e.g., a layer of the tissue, a tumor core, a necrotic area)
- A gene marking a cell type that physically clusters together (e.g., epithelial cells forming glands)
- A gene involved in a localized process (e.g., inflammation at a wound site, hypoxia in a tumor center)

**Low Moran's I (genes with random or uniform expression values):**

These genes are either:
- Ubiquitously expressed (housekeeping genes with similar values everywhere)
- Expressed by cell types that are evenly dispersed throughout the tissue
- Truly random noise with no spatial organization in their expression values

**Why this matters:**

Identifying genes with high spatial autocorrelation helps you discover:
- Tissue zonation and compartment boundaries
- Disease hotspots or niches
- Genes that drive or mark spatial organization
- Candidate region-specific markers without needing prior knowledge

In short: Moran's I takes the expression values of each gene and tells you whether those values are spatially organized, revealing biological territories across your tissue.


<img src="figures/visium_bc_CRISP3.png?v=2" width="33%"><img src="figures/visium_bc_CPB1.png?v=2" width="33%"><img src="figures/visium_bc_IGKC.png?v=2" width="33%">
<img src="figures/visium_bc_IGLC2.png?v=2" width="33%"><img src="figures/visium_bc_IGHG3.png?v=2" width="33%"><img src="figures/visium_bc_MALAT1.png?v=2" width="33%">

[View Spatially Variable Genes (Moran’s I)](https://docs.google.com/spreadsheets/d/18yUKNfJoli7Hmwz5dQvJJVYpO5k_upeTbhsZso9qO24/edit?usp=sharing) 

## References

- **Scanpy** – single-cell analysis in Python: [Wolf et al., Genome Biology 2018](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1382-0)
- **Squidpy** – spatial single-cell analysis in Python: [Palla et al., Nature Methods 2022](https://www.nature.com/articles/s41592-021-01322-8)
