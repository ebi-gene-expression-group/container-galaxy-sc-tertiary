# Available Tools

HSCiAP relies on published analysis tools for single downstream analysis, starting
from a quantified matrix (usually in 10x format). In the process of wrapping
those tools, we have made an effort to isolate its main analysis steps, towards
an scenario where different analysis components from each tool can be used interchangeably
provided by that their adherence to exchange formats, or intermediate converters, allow so.
Tables contain all the tools currently wrapped, with links to their sorce codes and our wrappers.
The following section discusses the merits of each tool and the steps that we extracted.

# Scanpy

From its GitHub repository, Scanpy is defined as:

> Scanpy is a scalable toolkit for analyzing single-cell gene expression data built jointly with anndata. It includes preprocessing, visualization, clustering, trajectory inference and differential expression testing. The Python-based implementation efficiently deals with datasets of more than one million cells.

Scanpy is used as a Python library, so Python code needs to be written to use it directly. We have written
a CLI component named scanpy-scripts which allows its usage from the command line, avoiding the need
of having to write Python code to use it for the wrapped functionality. The section Scanpy Steps available
below details all the steps wrapped.

## Scanpy Exchange formats

- 10x: initial reading and as output in some steps.
- AnnData: read/write through most modules.
- Loom: read/write through most modules.

## Scanpy modules available

Table 1 below details all modules available as both a CLI component and as a Galaxy wrapper. The name of the module links to an active Galaxy instance where that can be used, and TS provides a second link to the module in the Galaxy Toolshed.
Each module is linked to one of the cli-layers to one or more of the relevant analysis areas: Clustering (**C**), Differential expression/Marker detection (**DE-MD**), Trajectories (**T**), Cell type alignment (**CT**), quality control (**QC**) and Dimensionality reduction (**DR**).

| Module | Description | cli-layer | Analysis areas |
|--------|-------------|-----------|----------------|
| [scanpy_filter_cells](https://humancellatlas.usegalaxy.eu/tool_runner?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Febi-gxa%2Fscanpy_filter_cells%2Fscanpy_filter_cells)<sup>[TS](https://toolshed.g2.bx.psu.edu/view/ebi-gxa/scanpy_filter_cells) | based on counts and numbers of genes expressed | [scanpy-scripts]([scanpy-scripts](https://github.com/ebi-gene-expression-group/scanpy-scripts/tree/master)) | QC |
| [scanpy_filter_genes](https://humancellatlas.usegalaxy.eu/tool_runner?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Febi-gxa%2Fscanpy_filter_genes%2Fscanpy_filter_genes)<sup>[TS](https://toolshed.g2.bx.psu.edu/view/ebi-gxa/scanpy_filter_genes) | based on counts and numbers of cells expressed | [scanpy-scripts]([scanpy-scripts](https://github.com/ebi-gene-expression-group/scanpy-scripts/tree/master)) | QC |
| [scanpy_find_cluster](https://humancellatlas.usegalaxy.eu/tool_runner?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Febi-gxa%2Fscanpy_find_cluster%2Fscanpy_find_cluster)<sup>[TS](https://toolshed.g2.bx.psu.edu/view/ebi-gxa/scanpy_find_cluster) | based on community detection on KNN graph | [scanpy-scripts]([scanpy-scripts](https://github.com/ebi-gene-expression-group/scanpy-scripts/tree/master)) | C |
| [scanpy_find_markers](https://humancellatlas.usegalaxy.eu/tool_runner?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Febi-gxa%2Fscanpy_find_markers%2Fscanpy_find_markers)<sup>[TS](https://toolshed.g2.bx.psu.edu/view/ebi-gxa/scanpy_find_markers) | to find differentially expressed genes between groups | [scanpy-scripts]([scanpy-scripts](https://github.com/ebi-gene-expression-group/scanpy-scripts/tree/master)) | DE-MD |
| [scanpy_find_variable_genes](https://humancellatlas.usegalaxy.eu/tool_runner?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Febi-gxa%2Fscanpy_find_variable_genes%2Fscanpy_find_variable_genes)<sup>[TS](https://toolshed.g2.bx.psu.edu/view/ebi-gxa/scanpy_find_variable_genes) | based on normalised dispersion of expression | [scanpy-scripts]([scanpy-scripts](https://github.com/ebi-gene-expression-group/scanpy-scripts/tree/master)) | P |
| [scanpy_compute_graph](https://humancellatlas.usegalaxy.eu/tool_runner?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Febi-gxa%2Fscanpy_compute_graph%2Fscanpy_compute_graph)<sup>[TS](https://toolshed.g2.bx.psu.edu/view/ebi-gxa/scanpy_compute_graph) | to derive kNN graph | [scanpy-scripts]([scanpy-scripts](https://github.com/ebi-gene-expression-group/scanpy-scripts/tree/master)) | C T |
| [scanpy_normalise_data](https://humancellatlas.usegalaxy.eu/tool_runner?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Febi-gxa%2Fscanpy_normalise_data%2Fscanpy_normalise_data)<sup>[TS](https://toolshed.g2.bx.psu.edu/view/ebi-gxa/scanpy_normalise_data) | to make all cells having the same total expression | [scanpy-scripts]([scanpy-scripts](https://github.com/ebi-gene-expression-group/scanpy-scripts/tree/master)) | P |
| [scanpy_parameter_iterator](https://humancellatlas.usegalaxy.eu/tool_runner?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Febi-gxa%2Fscanpy_parameter_iterator%2Fscanpy_parameter_iterator)<sup>[TS](https://toolshed.g2.bx.psu.edu/view/ebi-gxa/scanpy_parameter_iterator) | produce an iteration over a defined parameter | [scanpy-scripts]([scanpy-scripts](https://github.com/ebi-gene-expression-group/scanpy-scripts/tree/master)) |  |
| [scanpy_plot_embed](https://humancellatlas.usegalaxy.eu/tool_runner?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Febi-gxa%2Fscanpy_plot_embed%2Fscanpy_plot_embed)<sup>[TS](https://toolshed.g2.bx.psu.edu/view/ebi-gxa/scanpy_plot_embed) | visualise cell embeddings | [scanpy-scripts]([scanpy-scripts](https://github.com/ebi-gene-expression-group/scanpy-scripts/tree/master)) | C |
| [scanpy_plot_trajectory](https://humancellatlas.usegalaxy.eu/tool_runner?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Febi-gxa%2Fscanpy_plot_trajectory%2Fscanpy_plot_trajectory)<sup>[TS](https://toolshed.g2.bx.psu.edu/view/ebi-gxa/scanpy_plot_trajectory) | visualise cell trajectories | [scanpy-scripts]([scanpy-scripts](https://github.com/ebi-gene-expression-group/scanpy-scripts/tree/master)) | T |
| [scanpy_read_10x](https://humancellatlas.usegalaxy.eu/tool_runner?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Febi-gxa%2Fscanpy_read_10x%2Fscanpy_read_10x)<sup>[TS](https://toolshed.g2.bx.psu.edu/view/ebi-gxa/scanpy_read_10x) | into hdf5 object handled by scanpy | [scanpy-scripts]([scanpy-scripts](https://github.com/ebi-gene-expression-group/scanpy-scripts/tree/master)) |  |
| [scanpy_regress_variable](https://humancellatlas.usegalaxy.eu/tool_runner?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Febi-gxa%2Fscanpy_regress_variable%2Fscanpy_regress_variable)<sup>[TS](https://toolshed.g2.bx.psu.edu/view/ebi-gxa/scanpy_regress_variable) | variables that might introduce batch effect | [scanpy-scripts]([scanpy-scripts](https://github.com/ebi-gene-expression-group/scanpy-scripts/tree/master)) | P |
| [scanpy_run_diffmap](https://humancellatlas.usegalaxy.eu/tool_runner?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Febi-gxa%2Fscanpy_run_diffmap%2Fscanpy_run_diffmap)<sup>[TS](https://toolshed.g2.bx.psu.edu/view/ebi-gxa/scanpy_run_diffmap) | calculate diffusion components | [scanpy-scripts]([scanpy-scripts](https://github.com/ebi-gene-expression-group/scanpy-scripts/tree/master)) | T |
| [scanpy_run_dpt](https://humancellatlas.usegalaxy.eu/tool_runner?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Febi-gxa%2Fscanpy_run_dpt%2Fscanpy_run_dpt)<sup>[TS](https://toolshed.g2.bx.psu.edu/view/ebi-gxa/scanpy_run_dpt) | diffusion pseudotime inference | [scanpy-scripts]([scanpy-scripts](https://github.com/ebi-gene-expression-group/scanpy-scripts/tree/master)) | T |
| [scanpy_run_fdg](https://humancellatlas.usegalaxy.eu/tool_runner?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Febi-gxa%2Fscanpy_run_fdg%2Fscanpy_run_fdg)<sup>[TS](https://toolshed.g2.bx.psu.edu/view/ebi-gxa/scanpy_run_fdg) | visualise cell clusters using force-directed graph | [scanpy-scripts]([scanpy-scripts](https://github.com/ebi-gene-expression-group/scanpy-scripts/tree/master)) |  |
| [scanpy_run_paga](https://humancellatlas.usegalaxy.eu/tool_runner?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Febi-gxa%2Fscanpy_run_paga%2Fscanpy_run_paga)<sup>[TS](https://toolshed.g2.bx.psu.edu/view/ebi-gxa/scanpy_run_paga) | trajectory inference | [scanpy-scripts]([scanpy-scripts](https://github.com/ebi-gene-expression-group/scanpy-scripts/tree/master)) | T |
| [scanpy_run_pca](https://humancellatlas.usegalaxy.eu/tool_runner?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Febi-gxa%2Fscanpy_run_pca%2Fscanpy_run_pca)<sup>[TS](https://toolshed.g2.bx.psu.edu/view/ebi-gxa/scanpy_run_pca) | for dimensionality reduction | [scanpy-scripts]([scanpy-scripts](https://github.com/ebi-gene-expression-group/scanpy-scripts/tree/master)) | DR |
| [scanpy_run_tsne](https://humancellatlas.usegalaxy.eu/tool_runner?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Febi-gxa%2Fscanpy_run_tsne%2Fscanpy_run_tsne)<sup>[TS](https://toolshed.g2.bx.psu.edu/view/ebi-gxa/scanpy_run_tsne) | visualise cell clusters using tSNE | [scanpy-scripts]([scanpy-scripts](https://github.com/ebi-gene-expression-group/scanpy-scripts/tree/master)) | DR |
| [scanpy_run_umap](https://humancellatlas.usegalaxy.eu/tool_runner?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Febi-gxa%2Fscanpy_run_umap%2Fscanpy_run_umap)<sup>[TS](https://toolshed.g2.bx.psu.edu/view/ebi-gxa/scanpy_run_umap) | visualise cell clusters using UMAP | [scanpy-scripts]([scanpy-scripts](https://github.com/ebi-gene-expression-group/scanpy-scripts/tree/master)) | DR |
| [scanpy_scale_data](https://humancellatlas.usegalaxy.eu/tool_runner?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Febi-gxa%2Fscanpy_scale_data%2Fscanpy_scale_data)<sup>[TS](https://toolshed.g2.bx.psu.edu/view/ebi-gxa/scanpy_scale_data) | to make expression variance the same for all genes | [scanpy-scripts]([scanpy-scripts](https://github.com/ebi-gene-expression-group/scanpy-scripts/tree/master)) | P |

## Examples of Workflows using Scanpy

# Seurat

According to its website, Seurat is:

> ...an R package designed for QC, analysis, and exploration of single-cell RNA-seq data. Seurat aims to enable users to identify and interpret sources of heterogeneity from single-cell transcriptomic measurements, and to integrate diverse types of single-cell data.

Seurat, as any R package, requires the user to write R code in order to use it to analyze data. We have written a CLI component, Seurat scripts, which exposes relevant functionality through a command line interface, enabling its use without having to write R code.

## Seurat exchange formats

Seurat version 3 wrapped supports the following exchange formats:

- 10x (input)
- Loom (input and output in most steps)
- AnnData (input for most steps)
- SingleCellExperiment (input and output for most steps) as an RDS R object.
- Seurat native format as an RDS R object.

## Seurat modules available

Table 2 below details all modules available as both a CLI component and as a Galaxy wrapper. The name of the module links to an active Galaxy instance where that can be used, and TS provides a second link to the module in the Galaxy Toolshed.
Each module is linked to one of the cli-layers to one or more of the relevant analysis areas: Clustering (**C**), Differential expression/Marker detection (**DE-MD**), Trajectories (**T**), Cell type alignment (**CT**), quality control (**QC**) and Dimensionality reduction (**DR**).


| Module | Description | cli-layer | Analysis areas |
|--------|-------------|-----------|----------------|
| [seurat_dim_plot](https://humancellatlas.usegalaxy.eu/tool_runner?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Febi-gxa%2Fseurat_dim_plot%2Fseurat_dim_plot)<sup>[TS](https://toolshed.g2.bx.psu.edu/view/ebi-gxa/seurat_dim_plot) | graphs the output of a dimensional reduction technique (PCA by default). Cells are colored by their identity class. | [seurat-scripts]([seurat-scripts](https://github.com/ebi-gene-expression-group/seurat-scripts/tree/master)) | DR |
| [seurat_export_cellbrowser](https://humancellatlas.usegalaxy.eu/tool_runner?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Febi-gxa%2Fseurat_export_cellbrowser%2Fseurat_export_cellbrowser)<sup>[TS](https://toolshed.g2.bx.psu.edu/view/ebi-gxa/seurat_export_cellbrowser) | produces files for UCSC CellBrowser import. | [seurat-scripts]([seurat-scripts](https://github.com/ebi-gene-expression-group/seurat-scripts/tree/master)) |  |
| [seurat_filter_cells](https://humancellatlas.usegalaxy.eu/tool_runner?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Febi-gxa%2Fseurat_filter_cells%2Fseurat_filter_cells)<sup>[TS](https://toolshed.g2.bx.psu.edu/view/ebi-gxa/seurat_filter_cells) | filter cells in a Seurat object | [seurat-scripts]([seurat-scripts](https://github.com/ebi-gene-expression-group/seurat-scripts/tree/master)) |  |
| [seurat_find_clusters](https://humancellatlas.usegalaxy.eu/tool_runner?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Febi-gxa%2Fseurat_find_clusters%2Fseurat_find_clusters)<sup>[TS](https://toolshed.g2.bx.psu.edu/view/ebi-gxa/seurat_find_clusters) | find clusters of cells | [seurat-scripts]([seurat-scripts](https://github.com/ebi-gene-expression-group/seurat-scripts/tree/master)) | C |
| [seurat_find_markers](https://humancellatlas.usegalaxy.eu/tool_runner?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Febi-gxa%2Fseurat_find_markers%2Fseurat_find_markers)<sup>[TS](https://toolshed.g2.bx.psu.edu/view/ebi-gxa/seurat_find_markers) | find markers (differentially expressed genes) | [seurat-scripts]([seurat-scripts](https://github.com/ebi-gene-expression-group/seurat-scripts/tree/master)) | DE-MD |
| [seurat_find_variable_genes](https://humancellatlas.usegalaxy.eu/tool_runner?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Febi-gxa%2Fseurat_find_variable_genes%2Fseurat_find_variable_genes)<sup>[TS](https://toolshed.g2.bx.psu.edu/view/ebi-gxa/seurat_find_variable_genes) | identify variable genes | [seurat-scripts]([seurat-scripts](https://github.com/ebi-gene-expression-group/seurat-scripts/tree/master)) | C DE-MD |
| [seurat_normalise_data](https://humancellatlas.usegalaxy.eu/tool_runner?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Febi-gxa%2Fseurat_normalise_data%2Fseurat_normalise_data)<sup>[TS](https://toolshed.g2.bx.psu.edu/view/ebi-gxa/seurat_normalise_data) | normalise data | [seurat-scripts]([seurat-scripts](https://github.com/ebi-gene-expression-group/seurat-scripts/tree/master)) | C |
| [seurat_read10x](https://humancellatlas.usegalaxy.eu/tool_runner?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Febi-gxa%2Fseurat_read10x%2Fseurat_read10x)<sup>[TS](https://toolshed.g2.bx.psu.edu/view/ebi-gxa/seurat_read10x) | Loads 10x data into a serialized seurat R object | [seurat-scripts]([seurat-scripts](https://github.com/ebi-gene-expression-group/seurat-scripts/tree/master)) |  |
| [seurat_run_pca](https://humancellatlas.usegalaxy.eu/tool_runner?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Febi-gxa%2Fseurat_run_pca%2Fseurat_run_pca)<sup>[TS](https://toolshed.g2.bx.psu.edu/view/ebi-gxa/seurat_run_pca) | run a PCA dimensionality reduction | [seurat-scripts]([seurat-scripts](https://github.com/ebi-gene-expression-group/seurat-scripts/tree/master)) | DR |
| [seurat_run_tsne](https://humancellatlas.usegalaxy.eu/tool_runner?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Febi-gxa%2Fseurat_run_tsne%2Fseurat_run_tsne)<sup>[TS](https://toolshed.g2.bx.psu.edu/view/ebi-gxa/seurat_run_tsne) | run t-SNE dimensionality reduction | [seurat-scripts]([seurat-scripts](https://github.com/ebi-gene-expression-group/seurat-scripts/tree/master)) | C DR |
| [seurat_scale_data](https://humancellatlas.usegalaxy.eu/tool_runner?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Febi-gxa%2Fseurat_scale_data%2Fseurat_scale_data)<sup>[TS](https://toolshed.g2.bx.psu.edu/view/ebi-gxa/seurat_scale_data) | scale and center genes | [seurat-scripts]([seurat-scripts](https://github.com/ebi-gene-expression-group/seurat-scripts/tree/master)) | C |
| [seurat_find_neighbours](https://humancellatlas.usegalaxy.eu/tool_runner?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Febi-gxa%2Fseurat_find_neighbours%2Fseurat_find_neighbours)<sup>[TS](https://toolshed.g2.bx.psu.edu/view/ebi-gxa/seurat_find_neighbours) | constructs a Shared Nearest Neighbor (SNN) Graph | [seurat-scripts]([seurat-scripts](https://github.com/ebi-gene-expression-group/seurat-scripts/tree/master)) | C |

# SC3

According to its website, Single-Cell Consensus Clustering (SC3) is:

> ...a tool for unsupervised clustering of scRNA-seq data. SC3 achieves high accuracy and robustness by consistently integrating different clustering solutions through a consensus approach. An interactive graphical implementation makes SC3 accessible to a wide audience of users. In addition, SC3 also aids biological interpretation by identifying marker genes, differentially expressed genes and outlier cells. A manuscript describing SC3 in details is published in Nature Methods.

SC3, as any R package, requires the user to write R code in order to use it to analyse data. We have written a CLI component, SC3 scripts, which exposes relevant functionality through a command line interface, enabling its use without having to write R code.

## SC3 exchange formats

SC3 wrapped supports the following exchange formats:

- 10x (input)
- SingleCellExperiment (input and output for most steps) as an RDS R object.

## SC3 modules available

Table 3 below details all modules available as both a CLI component and as a Galaxy wrapper. The name of the module links to an active Galaxy instance where that can be used, and TS provides a second link to the module in the Galaxy Toolshed.
Each module is linked to one of the cli-layers to one or more of the relevant analysis areas: Clustering (**C**), Differential expression/Marker detection (**DE-MD**), Trajectories (**T**), Cell type alignment (**CT**), quality control (**QC**) and Dimensionality reduction (**DR**).

| Module | Description | cli-layer | Analysis areas |
|--------|-------------|-----------|----------------|
| [sc3_calc_biology](https://humancellatlas.usegalaxy.eu/tool_runner?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Febi-gxa%2Fsc3_calc_biology%2Fsc3_calc_biology)<sup>[TS](https://toolshed.g2.bx.psu.edu/view/ebi-gxa/sc3_calc_biology) | calculates DE genes, marker genes and cell outliers | [sc3-scripts]([sc3-scripts](https://github.com/ebi-gene-expression-group/sc3-scripts/tree/master)) | DE-MD |
| [sc3_calc_consens](https://humancellatlas.usegalaxy.eu/tool_runner?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Febi-gxa%2Fsc3_calc_consens%2Fsc3_calc_consens)<sup>[TS](https://toolshed.g2.bx.psu.edu/view/ebi-gxa/sc3_calc_consens) | from multiple runs of k-means clustering | [sc3-scripts]([sc3-scripts](https://github.com/ebi-gene-expression-group/sc3-scripts/tree/master)) | C |
| [sc3_calc_dists](https://humancellatlas.usegalaxy.eu/tool_runner?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Febi-gxa%2Fsc3_calc_dists%2Fsc3_calc_dists)<sup>[TS](https://toolshed.g2.bx.psu.edu/view/ebi-gxa/sc3_calc_dists) | between cells | [sc3-scripts]([sc3-scripts](https://github.com/ebi-gene-expression-group/sc3-scripts/tree/master)) |  |
| [sc3_calc_transfs](https://humancellatlas.usegalaxy.eu/tool_runner?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Febi-gxa%2Fsc3_calc_transfs%2Fsc3_calc_transfs)<sup>[TS](https://toolshed.g2.bx.psu.edu/view/ebi-gxa/sc3_calc_transfs) | of distances using PCA and graph Laplacian | [sc3-scripts]([sc3-scripts](https://github.com/ebi-gene-expression-group/sc3-scripts/tree/master)) | C |
| [sc3_estimate_k](https://humancellatlas.usegalaxy.eu/tool_runner?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Febi-gxa%2Fsc3_estimate_k%2Fsc3_estimate_k)<sup>[TS](https://toolshed.g2.bx.psu.edu/view/ebi-gxa/sc3_estimate_k) | the number of clusters for k-means clustering | [sc3-scripts]([sc3-scripts](https://github.com/ebi-gene-expression-group/sc3-scripts/tree/master)) | C |
| [sc3_kmeans](https://humancellatlas.usegalaxy.eu/tool_runner?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Febi-gxa%2Fsc3_kmeans%2Fsc3_kmeans)<sup>[TS](https://toolshed.g2.bx.psu.edu/view/ebi-gxa/sc3_kmeans) | perform k-means clustering | [sc3-scripts]([sc3-scripts](https://github.com/ebi-gene-expression-group/sc3-scripts/tree/master)) | C |
| [sc3_prepare](https://humancellatlas.usegalaxy.eu/tool_runner?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Febi-gxa%2Fsc3_prepare%2Fsc3_prepare)<sup>[TS](https://toolshed.g2.bx.psu.edu/view/ebi-gxa/sc3_prepare) | a sc3 SingleCellExperiment object | [sc3-scripts]([sc3-scripts](https://github.com/ebi-gene-expression-group/sc3-scripts/tree/master)) |  |

# Scater

According to its website, Scater:

> ...contains tools to help with the analysis of single-cell transcriptomic data, focusing on low-level steps such as quality control, normalization and visualization. It is based on the SingleCellExperiment class (from the SingleCellExperiment package), and thus is interoperable with many other Bioconductor packages such as scran, batchelor and iSEE.

Scater, as any R package, requires the user to write R code in order to use it to analyse data. We have written a CLI component, scater scripts, which exposes relevant functionality through a command line interface, enabling its use without having to write R code.

## Scater exchange formats

Scater wrapped supports the following exchange formats:

- 10x (input)
- SingleCellExperiment (input and output for most steps) as an RDS R object.

## Scater modules available

Table 4 below details all modules available as both a CLI component and as a Galaxy wrapper. The name of the module links to an active Galaxy instance where that can be used, and TS provides a second link to the module in the Galaxy Toolshed.
Each module is linked to one of the cli-layers to one or more of the relevant analysis areas: Clustering (**C**), Differential expression/Marker detection (**DE-MD**), Trajectories (**T**), Cell type alignment (**CT**), quality control (**QC**) and Dimensionality reduction (**DR**).

| Module | Description | cli-layer | Analysis areas |
|--------|-------------|-----------|----------------|
| [scater_calculate_cpm](https://humancellatlas.usegalaxy.eu/tool_runner?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Febi-gxa%2Fscater_calculate_cpm%2Fscater_calculate_cpm)<sup>[TS](https://toolshed.g2.bx.psu.edu/view/ebi-gxa/scater_calculate_cpm) | from raw counts | [scater-scripts]([scater-scripts](https://github.com/ebi-gene-expression-group/scater-scripts/tree/master)) |  |
| [scater_calculate_qc_metrics](https://humancellatlas.usegalaxy.eu/tool_runner?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Febi-gxa%2Fscater_calculate_qc_metrics%2Fscater_calculate_qc_metrics)<sup>[TS](https://toolshed.g2.bx.psu.edu/view/ebi-gxa/scater_calculate_qc_metrics) | based on expression values and experiment information | [scater-scripts]([scater-scripts](https://github.com/ebi-gene-expression-group/scater-scripts/tree/master)) |  |
| [scater_filter](https://humancellatlas.usegalaxy.eu/tool_runner?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Febi-gxa%2Fscater_filter%2Fscater_filter)<sup>[TS](https://toolshed.g2.bx.psu.edu/view/ebi-gxa/scater_filter) | cells and genes based on pre-calculated stats and QC metrics | [scater-scripts]([scater-scripts](https://github.com/ebi-gene-expression-group/scater-scripts/tree/master)) |  |
| [scater_is_outlier](https://humancellatlas.usegalaxy.eu/tool_runner?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Febi-gxa%2Fscater_is_outlier%2Fscater_is_outlier)<sup>[TS](https://toolshed.g2.bx.psu.edu/view/ebi-gxa/scater_is_outlier) | cells based on expression metrics | [scater-scripts]([scater-scripts](https://github.com/ebi-gene-expression-group/scater-scripts/tree/master)) |  |
| [scater_normalize](https://humancellatlas.usegalaxy.eu/tool_runner?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Febi-gxa%2Fscater_normalize%2Fscater_normalize)<sup>[TS](https://toolshed.g2.bx.psu.edu/view/ebi-gxa/scater_normalize) | expression values by library size in log2 scale | [scater-scripts]([scater-scripts](https://github.com/ebi-gene-expression-group/scater-scripts/tree/master)) | C |
| [scater_read_10x_results](https://humancellatlas.usegalaxy.eu/tool_runner?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Febi-gxa%2Fscater_read_10x_results%2Fscater_read_10x_results)<sup>[TS](https://toolshed.g2.bx.psu.edu/view/ebi-gxa/scater_read_10x_results) | Loads 10x data into a serialized scater R object | [scater-scripts]([scater-scripts](https://github.com/ebi-gene-expression-group/scater-scripts/tree/master)) |  |


# Monocle3

According to its website, Monocle3:

> ...can help you perform three main types of analysis:
- Clustering, classifying, and counting cells. Single-cell RNA-Seq experiments allow you to discover new (and possibly rare) subtypes of cells. Monocle 3 helps you identify them.
- Constructing single-cell trajectories. In development, disease, and throughout life, cells transition from one state to another. Monocle 3 helps you discover these transitions.
- Differential expression analysis. Characterizing new cell types and states begins with comparisons to other, better understood cells. Monocle 3 includes a sophisticated, but easy-to-use system for differential expression.

Monocle3, as any R package, requires the user to write R code in order to use it to analyse data. We have written a CLI component, monocle3 scripts, which exposes relevant functionality through a command line interface, enabling its use without having to write R code.

## Monocle3 exchange formats

Monocle3 wrapped supports the following exchange formats:

- TSV, CSV or RDS Gene expression matrix, plus metadata

## Monocle3 modules available

Table 5 below details all modules available as both a CLI component and as a Galaxy wrapper. The name of the module links to an active Galaxy instance where that can be used, and TS provides a second link to the module in the Galaxy Toolshed.
Each module is linked to one of the cli-layers to one or more of the relevant analysis areas: Clustering (**C**), Differential expression/Marker detection (**DE-MD**), Trajectories (**T**), Cell type alignment (**CT**), quality control (**QC**) and Dimensionality reduction (**DR**).

| Module | Description | cli-layer | Analysis areas |
|--------|-------------|-----------|----------------|
| [monocle3_create](https://humancellatlas.usegalaxy.eu/tool_runner?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Febi-gxa%2Fmonocle3_create%2Fmonocle3_create)<sup>[TS](https://toolshed.g2.bx.psu.edu/view/ebi-gxa/monocle3_create) | a Monocle3 object from input data | [monocle-scripts]([monocle-scripts](https://github.com/ebi-gene-expression-group/monocle-scripts/tree/master)) |  |
| [monocle3_diffExp](https://humancellatlas.usegalaxy.eu/tool_runner?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Febi-gxa%2Fmonocle3_diffExp%2Fmonocle3_diffExp)<sup>[TS](https://toolshed.g2.bx.psu.edu/view/ebi-gxa/monocle3_diffExp) |  of genes along a trajectory | [monocle-scripts]([monocle-scripts](https://github.com/ebi-gene-expression-group/monocle-scripts/tree/master)) | DE-MD |
| [monocle3_learnGraph](https://humancellatlas.usegalaxy.eu/tool_runner?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Febi-gxa%2Fmonocle3_learnGraph%2Fmonocle3_learnGraph)<sup>[TS](https://toolshed.g2.bx.psu.edu/view/ebi-gxa/monocle3_learnGraph) | between cells in dimensionality reduced space | [monocle-scripts]([monocle-scripts](https://github.com/ebi-gene-expression-group/monocle-scripts/tree/master)) | T |
| [monocle3_orderCells](https://humancellatlas.usegalaxy.eu/tool_runner?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Febi-gxa%2Fmonocle3_orderCells%2Fmonocle3_orderCells)<sup>[TS](https://toolshed.g2.bx.psu.edu/view/ebi-gxa/monocle3_orderCells) | along trajectories | [monocle-scripts]([monocle-scripts](https://github.com/ebi-gene-expression-group/monocle-scripts/tree/master)) | T |
| [monocle3_partition](https://humancellatlas.usegalaxy.eu/tool_runner?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Febi-gxa%2Fmonocle3_partition%2Fmonocle3_partition)<sup>[TS](https://toolshed.g2.bx.psu.edu/view/ebi-gxa/monocle3_partition) | of cells into groups | [monocle-scripts]([monocle-scripts](https://github.com/ebi-gene-expression-group/monocle-scripts/tree/master)) | C T |
| [monocle3_plotCells](https://humancellatlas.usegalaxy.eu/tool_runner?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Febi-gxa%2Fmonocle3_plotCells%2Fmonocle3_plotCells)<sup>[TS](https://toolshed.g2.bx.psu.edu/view/ebi-gxa/monocle3_plotCells) |  in the reduced dimensionality space | [monocle-scripts]([monocle-scripts](https://github.com/ebi-gene-expression-group/monocle-scripts/tree/master)) | T |
| [monocle3_preprocess](https://humancellatlas.usegalaxy.eu/tool_runner?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Febi-gxa%2Fmonocle3_preprocess%2Fmonocle3_preprocess)<sup>[TS](https://toolshed.g2.bx.psu.edu/view/ebi-gxa/monocle3_preprocess) | a Monocle3 object to an initially dimensionally reduced space | [monocle-scripts]([monocle-scripts](https://github.com/ebi-gene-expression-group/monocle-scripts/tree/master)) | P |
| [monocle3_reduceDim](https://humancellatlas.usegalaxy.eu/tool_runner?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Febi-gxa%2Fmonocle3_reduceDim%2Fmonocle3_reduceDim)<sup>[TS](https://toolshed.g2.bx.psu.edu/view/ebi-gxa/monocle3_reduceDim) | for downstream analysis | [monocle-scripts]([monocle-scripts](https://github.com/ebi-gene-expression-group/monocle-scripts/tree/master)) | DR |


# SCMap

According to its website, SCMap is:

> ...a method for projecting cells from a scRNA-seq experiment onto the cell-types or individual cells identified in other experiments (the application can be run for free, without restrictions, from http://www.hemberg-lab.cloud/scmap).

SCMap, as any R package, requires the user to write R code in order to use it to analyse data. We have written a CLI component, scmap scripts, which exposes relevant functionality through a command line interface, enabling its use without having to write R code.

## SCMap exchange formats

SCMap wrapped supports the following exchange formats:

- SingleCellExperiment

## SCMap modules available

Table 6 below details all modules available as both a CLI component and as a Galaxy wrapper. The name of the module links to an active Galaxy instance where that can be used, and TS provides a second link to the module in the Galaxy Toolshed.
Each module is linked to one of the cli-layers to one or more of the relevant analysis areas: Clustering (**C**), Differential expression/Marker detection (**DE-MD**), Trajectories (**T**), Cell type alignment (**CT**), quality control (**QC**) and Dimensionality reduction (**DR**).

| Module | Description | cli-layer | Analysis areas |
|--------|-------------|-----------|----------------|
| [scmap_get_std_output](https://humancellatlas.usegalaxy.eu/tool_runner?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Febi-gxa%2Fscmap_get_std_output%2Fscmap_get_std_output)<sup>[TS](https://toolshed.g2.bx.psu.edu/view/ebi-gxa/scmap_get_std_output) | Create final output in standard format to allow for downstream analysis of predicted labels by tools of the EBI gene expression group's cell-types-analysis package | [scmap-cli]([scmap-cli](https://github.com/ebi-gene-expression-group/scmap-cli/tree/master)) |  |
| [scmap_index_cell](https://humancellatlas.usegalaxy.eu/tool_runner?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Febi-gxa%2Fscmap_index_cell%2Fscmap_index_cell)<sup>[TS](https://toolshed.g2.bx.psu.edu/view/ebi-gxa/scmap_index_cell) | creates a cell index for a dataset to enable fast approximate nearest neighbour search | [scmap-cli]([scmap-cli](https://github.com/ebi-gene-expression-group/scmap-cli/tree/master)) | CT |
| [scmap_index_cluster](https://humancellatlas.usegalaxy.eu/tool_runner?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Febi-gxa%2Fscmap_index_cluster%2Fscmap_index_cluster)<sup>[TS](https://toolshed.g2.bx.psu.edu/view/ebi-gxa/scmap_index_cluster) | calculates centroids of each cell type and merges them into a single table | [scmap-cli]([scmap-cli](https://github.com/ebi-gene-expression-group/scmap-cli/tree/master)) | CT |
| [scmap_scmap_cell](https://humancellatlas.usegalaxy.eu/tool_runner?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Febi-gxa%2Fscmap_scmap_cell%2Fscmap_scmap_cell)<sup>[TS](https://toolshed.g2.bx.psu.edu/view/ebi-gxa/scmap_scmap_cell) | searches each cell in a query dataset for the nearest neighbours by cosine distance within a collection of reference datasets. | [scmap-cli]([scmap-cli](https://github.com/ebi-gene-expression-group/scmap-cli/tree/master)) | CT |
| [scmap_scmap_cluster](https://humancellatlas.usegalaxy.eu/tool_runner?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Febi-gxa%2Fscmap_scmap_cluster%2Fscmap_scmap_cluster)<sup>[TS](https://toolshed.g2.bx.psu.edu/view/ebi-gxa/scmap_scmap_cluster) | projects one dataset to another | [scmap-cli]([scmap-cli](https://github.com/ebi-gene-expression-group/scmap-cli/tree/master)) | CT |
| [scmap_select_features](https://humancellatlas.usegalaxy.eu/tool_runner?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Febi-gxa%2Fscmap_select_features%2Fscmap_select_features)<sup>[TS](https://toolshed.g2.bx.psu.edu/view/ebi-gxa/scmap_select_features) | finds the most informative features (genes/transcripts) for projection. | [scmap-cli]([scmap-cli](https://github.com/ebi-gene-expression-group/scmap-cli/tree/master)) | CT |

# SCCAF

According to its website, Single Cell Clustering Assessment Framework (SCCAF) is:

>...a novel method for automated identification of putative cell types from single cell RNA-seq (scRNA-seq) data. By iteratively applying clustering and a machine learning approach to gene expression profiles of a given set of cells, SCCAF simultaneously identifies distinct cell groups and a weighted list of feature genes for each group. The feature genes, which are overexpressed in the particular cell group, jointly discriminate the given cell group from other cells. Each such group of cells corresponds to a putative cell type or state, characterised by the feature genes as markers.

SCCAF provides a CLI component, which exposes relevant functionality through a command line interface, enabling its use without having to write Python code.

## SCCAF exchange formats

SCCAF wrapped supports the following exchange formats:

- AnnData (input and output)
- Loom (input)

## SCCAF modules available

Table 7 below details all modules available as both a CLI component and as a Galaxy wrapper. The name of the module links to an active Galaxy instance where that can be used, and TS provides a second link to the module in the Galaxy Toolshed.
Each module is linked to one of the cli-layers to one or more of the relevant analysis areas: Clustering (**C**), Differential expression/Marker detection (**DE-MD**), Trajectories (**T**), Cell type alignment (**CT**), quality control (**QC**) and Dimensionality reduction (**DR**).

| Module | Description | cli-layer | Analysis areas |
|--------|-------------|-----------|----------------|
| [run_sccaf](https://humancellatlas.usegalaxy.eu/tool_runner?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Febi-gxa%2Frun_sccaf%2Frun_sccaf)<sup>[TS](https://toolshed.g2.bx.psu.edu/view/ebi-gxa/run_sccaf) | to assess and optimise clustering | [None]([None](https://github.com/ebi-gene-expression-group/None/tree/master)) | C |
| [sccaf_asses](https://humancellatlas.usegalaxy.eu/tool_runner?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Febi-gxa%2Fsccaf_asses%2Fsccaf_asses)<sup>[TS](https://toolshed.g2.bx.psu.edu/view/ebi-gxa/sccaf_asses) | runs an assesment of an SCCAF optimisation result or an existing clustering. | [None]([None](https://github.com/ebi-gene-expression-group/None/tree/master)) |  |
| [sccaf_asses_merger](https://humancellatlas.usegalaxy.eu/tool_runner?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Febi-gxa%2Fsccaf_asses_merger%2Fsccaf_asses_merger)<sup>[TS](https://toolshed.g2.bx.psu.edu/view/ebi-gxa/sccaf_asses_merger) | brings together distributed assesments. | [None]([None](https://github.com/ebi-gene-expression-group/None/tree/master)) | C |
| [sccaf_regress_out](https://humancellatlas.usegalaxy.eu/tool_runner?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Febi-gxa%2Fsccaf_regress_out%2Fsccaf_regress_out)<sup>[TS](https://toolshed.g2.bx.psu.edu/view/ebi-gxa/sccaf_regress_out) | with multiple categorical keys on an AnnData object. | [None]([None](https://github.com/ebi-gene-expression-group/None/tree/master)) |  |
