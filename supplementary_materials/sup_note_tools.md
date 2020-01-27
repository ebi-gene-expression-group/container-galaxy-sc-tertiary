# Available Tools

HSCiAP relies on published analysis tools for single downstream analysis, starting
from a quantified matrix (usually in 10x format). In the process of wrapping
those tools, we have made an effort to isolate its main analysis steps, towards
an scenario where different analysis components from each tool can be used interchangeably
provided by that their adherence to exchange formats, or intermediate converters, allow so.
Table S1 contains all the tools currently wrapped, with links to their sorce codes and our wrappers.
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

## Examples of Workflows using Seurat

# SC3
