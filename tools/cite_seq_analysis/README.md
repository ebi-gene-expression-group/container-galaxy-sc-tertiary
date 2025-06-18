# CITE-Seq Analysis Tool

This Galaxy tool implements a comprehensive CITE-Seq analysis pipeline using the scverse ecosystem, including muon for multimodal data management and scvi-tools for advanced integration.

## Overview

CITE-seq (Cellular Indexing of Transcriptomes and Epitopes by Sequencing) is a method that allows for the simultaneous measurement of both RNA transcripts and surface proteins at the single-cell level. This tool provides a complete pipeline for analyzing CITE-seq data, from data loading to final visualization.

## Features

- Supports both separate RNA/ADT inputs and combined MuData objects
- Assumes RNA data is pre-QC'ed, performs quality control on ADT modality and aligns cells between modalities
- Multiple ADT normalization methods (CLR, DSB)
- RNA normalization, HVG selection, and scaling
- Integration methods:
  - Feature concatenation (RNA PCs + normalized ADT expression)
  - TotalVI (advanced joint modeling of RNA and protein data) with automatic GPU detection and acceleration
- Clustering using Leiden algorithm
- UMAP visualization
- Multiple output formats (MuData, AnnData, TSV)

## totalVI Integration Features

The totalVI integration method provides powerful joint modeling of RNA and protein data:

- **Batch Correction**: Specify a batch key (column in observation metadata) to perform batch correction
- **Covariate Handling**: Include continuous or categorical covariates to account for known sources of variation
- **Flexible Model Configuration**: Configure the latent space dimensions and training parameters
- **Parameter Validation**: Automatic checking of metadata columns to ensure they exist in both modalities

## GPU Acceleration for totalVI Integration

The totalVI integration method can be computationally intensive, especially for larger datasets. This tool includes automatic GPU detection and utilization:

- **Automatic Detection**: The tool automatically checks for CUDA-compatible GPUs
- **Resource Optimization**: Training batch size is adjusted based on available GPU memory
- **Fallback Mechanism**: If GPU is not available or an error occurs, processing continues on CPU
- **Transparent Reporting**: GPU usage information is logged and stored in the output MuData object
- **Device Selection**: When multiple GPUs are available, specify which one to use

GPU acceleration can significantly reduce processing time, especially for datasets with many cells or when training for many epochs.

## Dependencies

- scanpy (>=1.2.0)
- muon (>=0.7.0)
- anndata (>=0.10.1)
- numpy (>=1.4.4)
- pandas (>=1.5.3)
- scvi-tools (>=0.11.0) (optional, for totalVI integration)
- pytorch (>=1.13.0) (optional, for GPU-accelerated totalVI)

## Usage

### Inputs

- RNA count matrix (AnnData/h5ad)
- ADT count matrix (AnnData/h5ad)
- Or a combined MuData object containing both modalities

The tool can accept input data in several ways:
1. File paths to RNA and ADT data
2. File path to a pre-made MuData object
3. In-memory AnnData objects for RNA and ADT
4. In-memory MuData object
5. Mixed approach: in-memory RNA AnnData with ADT data loaded from file

### Using as a Module in Python

```python
import scanpy as sc
import muon as mu
from cite_seq_analysis import run_cite_seq_pipeline

# Option 1: Using file paths
mdata_result = run_cite_seq_pipeline(
    rna_input="path/to/rna_data.h5ad",
    adt_input="path/to/adt_data.h5ad",
    output_file="path/to/output.h5mu",
    output_format="mudata"
)

# Option 2: Using in-memory AnnData objects
rna_adata = sc.read_h5ad("path/to/rna_data.h5ad")
adt_adata = sc.read_h5ad("path/to/adt_data.h5ad")

mdata_result = run_cite_seq_pipeline(
    rna_adata=rna_adata,
    adt_adata=adt_adata,
    output_file="path/to/output.h5mu",
    output_format="mudata"
)

# Option 3: Using in-memory MuData object
mdata = mu.read("path/to/mudata.h5mu")
mdata_result = run_cite_seq_pipeline(
    mdata=mdata,
    output_file="path/to/output.h5mu",
    output_format="mudata"
)

# Option 4: Mixed approach - in-memory RNA AnnData with ADT from file
mdata_result = run_cite_seq_pipeline(
    rna_adata=rna_adata,
    adt_input="path/to/adt_data.h5ad",
    output_file="path/to/output.h5mu",
    output_format="mudata"
)
```

### Outputs

- MuData object with all analysis results
- AnnData object with protein expression in observation metadata
- TSV file with cell metadata, including clustering and protein expression

### Parameters

This tool offers extensive parameters to customize each step of the analysis. See Galaxy interface for details.

## References

- muon: Bredikhin et al., "Muon: multimodal omics analysis framework." Genome Biology, 2021
- scanpy: Wolf et al., "SCANPY: large-scale single-cell gene expression data analysis." Genome Biology, 2018
- scvi-tools: Gayoso et al., "A Python library for probabilistic analysis of single-cell omics data." Nature Biotechnology, 2022
- CITE-seq: Stoeckius et al., "Simultaneous epitope and transcriptome measurement in single cells." Nature Methods, 2017

## License

MIT License
