#!/usr/bin/env python
"""
CITE-Seq Analysis Pipeline.

This script implements a comprehensive CITE-Seq analysis pipeline using scverse packages,
including muon for multimodal data management and optionally scvi-tools for advanced integration.
"""

import argparse
import logging
import os
import sys
from typing import Optional, Dict, Tuple, Union, List

import anndata as ad
import muon as mu
import numpy as np
import pandas as pd
import scanpy as sc

# Set up logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger("cite-seq-analysis")


def setup_arg_parser():
    """Set up and return argument parser for command line usage."""
    parser = argparse.ArgumentParser(
        description="CITE-Seq Analysis Pipeline using scverse ecosystem."
    )

    # Input options
    input_group = parser.add_argument_group("Input Options")
    input_group.add_argument(
        "--rna-input", 
        type=str, 
        help="Path to RNA count matrix (AnnData, h5ad) or 10x directory"
    )
    input_group.add_argument(
        "--adt-input", 
        type=str, 
        help="Path to ADT count matrix (AnnData, h5ad) or 10x directory"
    )
    input_group.add_argument(
        "--mudata-input", 
        type=str, 
        help="Path to combined MuData object (instead of separate RNA and ADT inputs)"
    )

    # QC parameters
    qc_group = parser.add_argument_group("Quality Control Parameters")
    qc_group.add_argument(
        "--rna-min-genes", 
        type=int, 
        default=200, 
        help="Minimum number of genes per cell for RNA"
    )
    qc_group.add_argument(
        "--rna-min-cells", 
        type=int, 
        default=3, 
        help="Minimum number of cells per gene for RNA"
    )
    qc_group.add_argument(
        "--max-mito-pct", 
        type=float, 
        default=20.0, 
        help="Maximum percentage of mitochondrial genes"
    )
    qc_group.add_argument(
        "--mito-prefix", 
        type=str, 
        default="MT-", 
        help="Prefix for mitochondrial genes"
    )
    qc_group.add_argument(
        "--adt-min-counts", 
        type=int, 
        default=1, 
        help="Minimum counts per cell for ADT"
    )

    # Normalization parameters
    norm_group = parser.add_argument_group("Normalization Parameters")
    norm_group.add_argument(
        "--adt-norm-method", 
        type=str, 
        choices=["dsb", "clr"], 
        default="clr",
        help="Normalization method for ADT data"
    )
    norm_group.add_argument(
        "--use-isotype-control", 
        action="store_true", 
        help="Use isotype controls for DSB normalization"
    )
    norm_group.add_argument(
        "--isotype-controls", 
        nargs="+", 
        default=[],
        help="List of isotype control names (for DSB normalization)"
    )
    norm_group.add_argument(
        "--rna-target-sum", 
        type=float, 
        default=1e4, 
        help="Target sum for RNA normalization"
    )
    norm_group.add_argument(
        "--normalize-rna",
        action="store_true",
        help="Force normalization of RNA data even if it appears to be pre-processed"
    )
    norm_group.add_argument(
        "--raw-counts-location",
        type=str,
        choices=["X", "raw.X", "layer"],
        default="X",
        help="Location of raw counts in RNA AnnData: 'X' (default), 'raw.X' (backed up raw counts), or 'layer' (in a specific layer)"
    )
    norm_group.add_argument(
        "--raw-counts-layer",
        type=str,
        default="counts",
        help="Name of the layer containing raw counts (only used if raw-counts-location is 'layer')"
    )

    # HVG and scaling parameters
    hvg_group = parser.add_argument_group("Highly Variable Genes Parameters")
    hvg_group.add_argument(
        "--hvg-min-mean", 
        type=float, 
        default=0.0125, 
        help="Minimum mean expression for HVG"
    )
    hvg_group.add_argument(
        "--hvg-max-mean", 
        type=float, 
        default=3, 
        help="Maximum mean expression for HVG"
    )
    hvg_group.add_argument(
        "--hvg-min-disp", 
        type=float, 
        default=0.5, 
        help="Minimum dispersion for HVG"
    )
    hvg_group.add_argument(
        "--scale-max-value", 
        type=float, 
        default=10, 
        help="Maximum value after scaling"
    )

    # Integration parameters
    integration_group = parser.add_argument_group("Integration Parameters")
    integration_group.add_argument(
        "--integration-method", 
        type=str, 
        choices=["totalVI", "feature-concat"], 
        default="feature-concat",
        help="Method for integrating RNA and ADT data"
    )
    integration_group.add_argument(
        "--n-latent", 
        type=int, 
        default=20, 
        help="Number of latent dimensions for totalVI"
    )
    integration_group.add_argument(
        "--max-epochs", 
        type=int, 
        default=400, 
        help="Maximum number of training epochs for totalVI"
    )
    integration_group.add_argument(
        "--n-pcs", 
        type=int, 
        default=50, 
        help="Number of principal components for RNA"
    )
    integration_group.add_argument(
        "--gpu-device",
        type=int,
        default=0,
        help="GPU device ID to use for totalVI integration (if multiple GPUs are available)"
    )
    integration_group.add_argument(
        "--batch-key",
        type=str,
        default=None,
        help="Column name in obs metadata containing batch information for totalVI"
    )
    integration_group.add_argument(
        "--continuous-covariates",
        nargs="+",
        default=None,
        help="List of column names in obs metadata with continuous covariates for totalVI"
    )
    integration_group.add_argument(
        "--categorical-covariates",
        nargs="+",
        default=None,
        help="List of column names in obs metadata with categorical covariates for totalVI"
    )
    integration_group.add_argument(
        "--integrated-dim-reduction-key",
        type=str,
        default="X_integrated",
        help="Key name to use for the integrated representation in obsm"
    )

    # Clustering parameters
    cluster_group = parser.add_argument_group("Clustering Parameters")
    cluster_group.add_argument(
        "--n-neighbors", 
        type=int, 
        default=15, 
        help="Number of neighbors for neighborhood graph"
    )
    cluster_group.add_argument(
        "--resolution", 
        type=float, 
        default=0.5, 
        help="Resolution for Leiden clustering"
    )
    cluster_group.add_argument(
        "--umap-min-dist", 
        type=float, 
        default=0.3, 
        help="Minimum distance parameter for UMAP"
    )

    # Output parameters
    output_group = parser.add_argument_group("Output Parameters")
    output_group.add_argument(
        "--output-format", 
        type=str, 
        choices=["mudata", "anndata", "tsv"], 
        default="mudata",
        help="Output format type"
    )
    output_group.add_argument(
        "--output-file", 
        type=str, 
        required=True, 
        help="Path for output file"
    )

    return parser


def load_rna_data(rna_input: str) -> ad.AnnData:
    """
    Load RNA data from a file or directory.
    
    Args:
        rna_input: Path to RNA data (AnnData or 10x directory)
        
    Returns:
        AnnData object with RNA data
    """
    logger.info(f"Loading RNA data from {rna_input}")
    if os.path.isdir(rna_input):
        try:
            # Try loading as 10x directory
            rna_adata = sc.read_10x_mtx(rna_input)
            rna_adata.var_names_make_unique()
        except Exception as e:
            logger.error(f"Error loading RNA data from 10x directory: {e}")
            sys.exit(1)
    else:
        try:
            # Try loading as h5ad file
            rna_adata = sc.read_h5ad(rna_input)
        except Exception as e:
            logger.error(f"Error loading RNA data from file: {e}")
            sys.exit(1)
    
    # Check if RNA data appears to be QC'ed
    if 'n_genes' not in rna_adata.obs and 'n_counts' not in rna_adata.obs:
        logger.warning("RNA data does not appear to have QC metrics. This pipeline assumes pre-QC'ed RNA data.")
    
    return rna_adata


def load_adt_data(adt_input: str) -> ad.AnnData:
    """
    Load ADT data from a file or directory.
    
    Args:
        adt_input: Path to ADT data (AnnData or 10x directory)
        
    Returns:
        AnnData object with ADT data
    """
    logger.info(f"Loading ADT data from {adt_input}")
    if os.path.isdir(adt_input):
        try:
            # Try loading as 10x directory
            adt_adata = sc.read_10x_mtx(adt_input)
            adt_adata.var_names_make_unique()
        except Exception as e:
            logger.error(f"Error loading ADT data from 10x directory: {e}")
            sys.exit(1)
    else:
        try:
            # Try loading as h5ad file
            adt_adata = sc.read_h5ad(adt_input)
        except Exception as e:
            logger.error(f"Error loading ADT data from file: {e}")
            sys.exit(1)
    
    return adt_adata


def load_mudata(mudata_input: str) -> mu.MuData:
    """
    Load a MuData object from a file.
    
    Args:
        mudata_input: Path to a MuData object
        
    Returns:
        MuData object with RNA and ADT modalities
    """
    logger.info(f"Loading MuData from {mudata_input}")
    mdata = mu.read(mudata_input)
    
    # Check if it has the expected modalities
    if "rna" not in mdata.mod or "prot" not in mdata.mod:
        logger.error("MuData object must have 'rna' and 'prot' modalities")
        sys.exit(1)
    
    # Check if RNA modality appears to be QC'ed
    if 'n_genes' not in mdata.mod["rna"].obs and 'n_counts' not in mdata.mod["rna"].obs:
        logger.warning("RNA modality does not appear to have QC metrics. This pipeline assumes pre-QC'ed RNA data.")
    
    return mdata


def create_mudata(rna_adata: ad.AnnData, adt_adata: ad.AnnData) -> mu.MuData:
    """
    Create a MuData object from RNA and ADT AnnData objects.
    
    Args:
        rna_adata: AnnData object with RNA data
        adt_adata: AnnData object with ADT data
        
    Returns:
        MuData object with RNA and ADT modalities
    """
    logger.info("Creating MuData object")
    return mu.MuData({"rna": rna_adata, "prot": adt_adata})


def load_data(
    rna_input: Optional[str] = None,
    adt_input: Optional[str] = None,
    mudata_input: Optional[str] = None,
    rna_adata: Optional[ad.AnnData] = None,
    adt_adata: Optional[ad.AnnData] = None,
    mdata: Optional[mu.MuData] = None,
) -> Tuple[mu.MuData, Optional[mu.MuData]]:
    """
    Load data from files or use provided objects to create a MuData object.
    
    This function can work with:
    1. A path to a MuData file
    2. Paths to RNA and ADT files
    3. In-memory AnnData objects for RNA and ADT
    4. An in-memory MuData object
    
    Args:
        rna_input: Path to RNA data (AnnData or 10x directory)
        adt_input: Path to ADT data (AnnData or 10x directory)
        mudata_input: Path to a pre-made MuData object
        rna_adata: AnnData object with RNA data (in memory)
        adt_adata: AnnData object with ADT data (in memory)
        mdata: MuData object with RNA and ADT modalities (in memory)
        
    Returns:
        A tuple containing:
            - MuData object with modalities
            - Raw MuData object for DSB normalization (or None if not needed)
    """
    result_mdata = None
    
    # Case 1: MuData object is provided in memory
    if mdata is not None:
        logger.info("Using provided MuData object")
        result_mdata = mdata
        
        # Check if it has the expected modalities
        if "rna" not in result_mdata.mod or "prot" not in result_mdata.mod:
            logger.error("MuData object must have 'rna' and 'prot' modalities")
            sys.exit(1)
        
        # Check if RNA modality appears to be QC'ed
        if 'n_genes' not in result_mdata.mod["rna"].obs and 'n_counts' not in result_mdata.mod["rna"].obs:
            logger.warning("RNA modality does not appear to have QC metrics. This pipeline assumes pre-QC'ed RNA data.")
    
    # Case 2: MuData file path is provided
    elif mudata_input is not None:
        result_mdata = load_mudata(mudata_input)
    
    # Case 3: Both RNA and ADT AnnData objects are provided in memory
    elif rna_adata is not None and adt_adata is not None:
        logger.info("Using provided RNA and ADT AnnData objects")
        
        # Check if RNA data appears to be QC'ed
        if 'n_genes' not in rna_adata.obs and 'n_counts' not in rna_adata.obs:
            logger.warning("RNA data does not appear to have QC metrics. This pipeline assumes pre-QC'ed RNA data.")
        
        result_mdata = create_mudata(rna_adata, adt_adata)
    
    # Case 4: File paths for both RNA and ADT are provided
    elif rna_input is not None and adt_input is not None:
        rna_adata = load_rna_data(rna_input)
        adt_adata = load_adt_data(adt_input)
        result_mdata = create_mudata(rna_adata, adt_adata)
    
    # Case 5: RNA AnnData in memory and ADT file path
    elif rna_adata is not None and adt_input is not None:
        logger.info("Using provided RNA AnnData object and loading ADT data from file")
        
        # Check if RNA data appears to be QC'ed
        if 'n_genes' not in rna_adata.obs and 'n_counts' not in rna_adata.obs:
            logger.warning("RNA data does not appear to have QC metrics. This pipeline assumes pre-QC'ed RNA data.")
        
        adt_adata = load_adt_data(adt_input)
        result_mdata = create_mudata(rna_adata, adt_adata)
    
    # Case 6: RNA file path and ADT AnnData in memory
    elif rna_input is not None and adt_adata is not None:
        logger.info("Loading RNA data from file and using provided ADT AnnData object")
        rna_adata = load_rna_data(rna_input)
        result_mdata = create_mudata(rna_adata, adt_adata)
    
    else:
        logger.error("Insufficient data provided. Need either: MuData object/file, both RNA and ADT data (as files or AnnData objects), or a combination.")
        sys.exit(1)
    
    # Create a copy with raw counts for DSB normalization
    mdata_raw = None
    if result_mdata is not None:
        mdata_raw = result_mdata.copy()
        
    return result_mdata, mdata_raw


def perform_qc(
    mdata: mu.MuData,
    rna_min_genes: int = 200,
    rna_min_cells: int = 3,
    max_mito_pct: float = 20.0,
    mito_prefix: str = "MT-",
    adt_min_counts: int = 1,
) -> mu.MuData:
    """
    Perform quality control primarily on ADT modality, assuming RNA is pre-QC'ed.
    Only performs joint QC for cell alignment between modalities.
    
    Args:
        mdata: MuData object with RNA and ADT modalities
        rna_min_genes: Minimum number of genes per cell (used only if joint QC needed)
        rna_min_cells: Minimum number of cells per gene (used only if joint QC needed)
        max_mito_pct: Maximum percentage of mitochondrial genes (used only if joint QC needed)
        mito_prefix: Prefix for mitochondrial genes (used only if joint QC needed)
        adt_min_counts: Minimum counts per cell for ADT
        
    Returns:
        MuData object after QC filtering
    """
    # Log RNA assumption
    logger.info("Assuming RNA data has already been QC'ed")
    
    # QC for ADT
    logger.info("Performing QC on ADT data")
    sc.pp.filter_cells(mdata.mod["prot"], min_counts=adt_min_counts)
    
    # Ensure cell alignment between modalities
    common_cells = np.intersect1d(mdata.mod["rna"].obs_names, mdata.mod["prot"].obs_names)
    logger.info(f"Aligning cells between modalities: keeping {len(common_cells)} cells common to both")
    
    if len(common_cells) == 0:
        logger.error("No common cells between RNA and ADT modalities after QC")
        sys.exit(1)
    
    mdata.mod["rna"] = mdata.mod["rna"][common_cells]
    mdata.mod["prot"] = mdata.mod["prot"][common_cells]
    
    return mdata


def normalize_adt(
    mdata: mu.MuData,
    mdata_raw: Optional[mu.MuData] = None,
    method: str = "clr",
    use_isotype_control: bool = False,
    isotype_controls: List[str] = None,
) -> mu.MuData:
    """
    Normalize ADT data using specified method.
    
    Args:
        mdata: MuData object
        mdata_raw: MuData object with raw counts (for DSB)
        method: Normalization method ('dsb' or 'clr')
        use_isotype_control: Whether to use isotype controls for DSB
        isotype_controls: List of isotype control names
        
    Returns:
        MuData object with normalized ADT data
    """
    if method == "dsb":
        if mdata_raw is None:
            logger.warning("Raw data not provided for DSB normalization, falling back to CLR")
            method = "clr"
        else:
            try:
                import muon.prot.pp as mpp
                logger.info("Performing DSB normalization on ADT data")
                
                if use_isotype_control and isotype_controls:
                    isotype_mask = mdata.mod["prot"].var_names.isin(isotype_controls)
                    if sum(isotype_mask) > 0:
                        mpp.dsb(
                            mdata, 
                            raw_mu_data=mdata_raw, 
                            use_isotype_control=True,
                            isotype_controls=isotype_controls
                        )
                    else:
                        logger.warning("Specified isotype controls not found, using DSB without controls")
                        mpp.dsb(mdata, raw_mu_data=mdata_raw, use_isotype_control=False)
                else:
                    mpp.dsb(mdata, raw_mu_data=mdata_raw, use_isotype_control=False)
            except ImportError:
                logger.warning("muon.prot.pp not available, falling back to CLR normalization")
                method = "clr"
            except Exception as e:
                logger.warning(f"DSB normalization failed: {e}, falling back to CLR")
                method = "clr"
    
    if method == "clr":
        try:
            import muon.prot.pp as mpp
            logger.info("Performing CLR normalization on ADT data")
            mpp.clr(mdata.mod["prot"])
        except ImportError:
            logger.error("muon.prot.pp not available for CLR normalization")
            sys.exit(1)
        except Exception as e:
            logger.error(f"CLR normalization failed: {e}")
            sys.exit(1)
            
    return mdata


def process_rna(
    mdata: mu.MuData,
    target_sum: float = 1e4,
    hvg_min_mean: float = 0.0125,
    hvg_max_mean: float = 3,
    hvg_min_disp: float = 0.5,
    scale_max_value: float = 10,
    n_pcs: int = 50,
    normalize_rna: bool = False,
    raw_counts_location: str = "X",
    raw_counts_layer: str = "counts",
) -> mu.MuData:
    """
    Normalize and process RNA data, respecting existing preprocessing.
    
    Args:
        mdata: MuData object
        target_sum: Target sum for normalization
        hvg_min_mean: Minimum mean expression for HVG
        hvg_max_mean: Maximum mean expression for HVG
        hvg_min_disp: Minimum dispersion for HVG
        scale_max_value: Maximum value after scaling
        n_pcs: Number of principal components
        normalize_rna: Whether to force normalization of RNA data
        raw_counts_location: Location of raw counts in RNA AnnData ('X', 'raw.X', or 'layer')
        raw_counts_layer: Name of the layer containing raw counts (only used if raw_counts_location is 'layer')
        
    Returns:
        MuData object with processed RNA data
    """
    # Get RNA AnnData for easier reference
    rna_adata = mdata.mod["rna"]
    
    # Check if data appears to be already preprocessed
    has_hvg = 'highly_variable' in rna_adata.var
    needs_processing = normalize_rna or not has_hvg
    
    # Ensure raw counts are available for totalVI
    if raw_counts_location == "X":
        # If raw counts are in X and we're normalizing, we need to back them up
        if needs_processing:
            logger.info("Raw counts are in X. Backing up raw counts before normalization.")
            # Create a copy in .raw to preserve raw counts
            if rna_adata.raw is None:
                rna_adata.raw = rna_adata.copy()
                logger.info("Raw counts backed up in AnnData.raw")
    elif raw_counts_location == "raw.X":
        # Check if .raw exists
        if rna_adata.raw is None:
            logger.warning("Raw counts specified to be in .raw.X, but .raw not found. Using X as raw counts.")
    elif raw_counts_location == "layer":
        # Check if specified layer exists
        if raw_counts_layer not in rna_adata.layers:
            logger.warning(f"Raw counts layer '{raw_counts_layer}' not found. Using X as raw counts.")
    
    # Process RNA (normalize, log, find HVGs) if needed or explicitly requested
    if needs_processing:
        logger.info("Processing RNA data")
        
        if normalize_rna:
            logger.info(f"Normalizing RNA data to target sum of {target_sum}")
            sc.pp.normalize_total(rna_adata, target_sum=target_sum)
            sc.pp.log1p(rna_adata)
        
        if not has_hvg:
            logger.info("Finding highly variable genes")
            sc.pp.highly_variable_genes(
                rna_adata,
                min_mean=hvg_min_mean,
                max_mean=hvg_max_mean,
                min_disp=hvg_min_disp
            )
        
        logger.info("Scaling RNA data")
        sc.pp.scale(rna_adata, max_value=scale_max_value)
    else:
        logger.info("RNA data appears to be pre-processed, skipping normalization and scaling")
    
    # Count HVGs for reference
    n_hvg = sum(rna_adata.var.highly_variable)
    logger.info(f"Using {n_hvg} highly variable genes")
    
    # Run PCA if not already computed
    if "X_pca" not in rna_adata.obsm:
        logger.info(f"Running PCA with {n_pcs} components")
        sc.tl.pca(rna_adata, n_comps=min(n_pcs, n_hvg))
    else:
        logger.info("PCA already computed for RNA data")
    
    return mdata


def integrate_modalities(
    mdata: mu.MuData,
    method: str = "feature-concat",
    n_latent: int = 20,
    max_epochs: int = 400,
    gpu_device: int = 0,
    batch_key: Optional[str] = None,
    continuous_covariate_keys: Optional[List[str]] = None,
    categorical_covariate_keys: Optional[List[str]] = None,
    integrated_dim_reduction_key: str = "X_integrated",
    raw_counts_location: str = "X",
    raw_counts_layer: str = "counts",
) -> mu.MuData:
    """
    Integrate RNA and ADT modalities.
    
    Args:
        mdata: MuData object with processed RNA and ADT data
        method: Integration method ('totalVI' or 'feature-concat')
        n_latent: Number of latent dimensions for totalVI
        max_epochs: Maximum number of training epochs for totalVI
        gpu_device: GPU device ID to use (if available)
        batch_key: Column name in obs metadata containing batch information
        continuous_covariate_keys: List of columns in obs metadata with continuous covariates
        categorical_covariate_keys: List of columns in obs metadata with categorical covariates
        integrated_dim_reduction_key: Key name to use for the integrated representation in obsm
        raw_counts_location: Location of raw counts in RNA AnnData ('X', 'raw.X', or 'layer')
        raw_counts_layer: Name of the layer containing raw counts (only used if raw_counts_location is 'layer')
        
    Returns:
        MuData object with integrated data
    """
    if method == "totalVI":
        try:
            import scvi
            logger.info("Using totalVI for integration")
            
            # Check for GPU availability
            try:
                import torch
                gpu_available = torch.cuda.is_available()
                gpu_count = torch.cuda.device_count() if gpu_available else 0
                
                if gpu_available:
                    # Check if specified GPU device is valid
                    if gpu_device >= 0 and gpu_device < gpu_count:
                        # Set the device for PyTorch
                        torch.cuda.set_device(gpu_device)
                        device_name = torch.cuda.get_device_name(gpu_device)
                        logger.info(f"Using GPU device {gpu_device}: {device_name}")
                    else:
                        if gpu_device >= gpu_count:
                            logger.warning(f"Specified GPU device {gpu_device} not available. "
                                          f"Using default device (0). Available devices: {gpu_count}")
                            gpu_device = 0
                            
                    gpu_names = [torch.cuda.get_device_name(i) for i in range(gpu_count)]
                    logger.info(f"GPU available: {gpu_available}, Count: {gpu_count}, Devices: {', '.join(gpu_names)}")
                    logger.info("Training totalVI with GPU acceleration")
                else:
                    logger.warning("No GPU detected. totalVI will run on CPU, which may be significantly slower.")
            except ImportError:
                logger.warning("Could not check GPU availability (torch not installed). Proceeding with totalVI.")
                gpu_available = False
            
            # Prepare raw counts for totalVI
            # totalVI expects raw counts for RNA in X, and normalized data for protein in X
            logger.info("Preparing raw counts for totalVI...")
            
            # Get direct reference to the RNA AnnData
            rna_adata = mdata.mod["rna"]
            
            # Backup the current X to a temporary layer to restore later
            # Only create backup if we'll actually modify X
            need_to_swap = False
            orig_X = None
            
            if raw_counts_location == "raw.X" and rna_adata.raw is not None:
                need_to_swap = True
                logger.info("Temporarily swapping raw counts from .raw.X to X for totalVI")
                # Make a deep copy of the original data matrix
                if isinstance(rna_adata.X, np.ndarray):
                    orig_X = rna_adata.X.copy()
                else:
                    # For sparse matrices
                    orig_X = rna_adata.X.copy()
            elif raw_counts_location == "layer" and raw_counts_layer in rna_adata.layers:
                need_to_swap = True
                logger.info(f"Temporarily swapping raw counts from layer '{raw_counts_layer}' to X for totalVI")
                # Make a deep copy of the original data matrix
                if isinstance(rna_adata.X, np.ndarray):
                    orig_X = rna_adata.X.copy()
                else:
                    # For sparse matrices
                    orig_X = rna_adata.X.copy()
            elif raw_counts_location == "X" and rna_adata.raw is not None:
                need_to_swap = True
                logger.info("Restoring raw counts from .raw.X to X for totalVI")
                # Make a deep copy of the original data matrix
                if isinstance(rna_adata.X, np.ndarray):
                    orig_X = rna_adata.X.copy()
                else:
                    # For sparse matrices
                    orig_X = rna_adata.X.copy()
                
            # Only swap data if needed
            if need_to_swap:
                if raw_counts_location == "raw.X" and rna_adata.raw is not None:
                    # Swap raw counts from raw.X to X for totalVI
                    rna_adata.X = rna_adata.raw.X
                elif raw_counts_location == "layer" and raw_counts_layer in rna_adata.layers:
                    # Swap raw counts from specified layer to X
                    rna_adata.X = rna_adata.layers[raw_counts_layer]
                elif raw_counts_location == "X" and rna_adata.raw is not None:
                    # Use raw counts from raw.X
                    rna_adata.X = rna_adata.raw.X
            else:
                # Using current X as is
                if raw_counts_location == "raw.X":
                    logger.warning("raw.X specified but not found. Using current X for totalVI (may not be raw counts)")
                elif raw_counts_location == "layer":
                    logger.warning(f"Layer '{raw_counts_layer}' not found. Using current X for totalVI (may not be raw counts)")
                else:
                    logger.info("Using current X as raw counts for totalVI")
            
            # Store information about raw counts source for metadata
            totalvi_raw_counts_info = {
                "source": raw_counts_location,
                "layer_name": raw_counts_layer if raw_counts_location == "layer" else None
            }
            
            # Validate and prepare batch and covariate parameters
            validated_batch_key = None
            validated_continuous_covs = []
            validated_categorical_covs = []
            
            # Check if batch key exists in both modalities
            if batch_key is not None:
                if batch_key in mdata.obs.columns:
                    logger.info(f"Using '{batch_key}' as batch key for totalVI")
                    validated_batch_key = batch_key
                else:
                    logger.warning(f"Specified batch key '{batch_key}' not found in data. Proceeding without batch correction.")
            
            # Check for continuous covariates
            if continuous_covariate_keys:
                for key in continuous_covariate_keys:
                    if key in mdata.obs.columns:
                        # Check if the column contains numeric data
                        if pd.api.types.is_numeric_dtype(mdata.obs[key]):
                            validated_continuous_covs.append(key)
                            logger.info(f"Using '{key}' as continuous covariate for totalVI")
                        else:
                            logger.warning(f"Specified continuous covariate '{key}' is not numeric. Skipping.")
                    else:
                        logger.warning(f"Specified continuous covariate '{key}' not found in data. Skipping.")
            
            # Check for categorical covariates
            if categorical_covariate_keys:
                for key in categorical_covariate_keys:
                    if key in mdata.obs.columns:
                        validated_categorical_covs.append(key)
                        logger.info(f"Using '{key}' as categorical covariate for totalVI")
                    else:
                        logger.warning(f"Specified categorical covariate '{key}' not found in data. Skipping.")
            
            # Setup anndata with validated parameters
            logger.info("Setting up AnnData for totalVI")
            scvi.model.TOTALVI.setup_anndata(
                mdata,
                batch_key=validated_batch_key,
                continuous_covariate_keys=validated_continuous_covs if validated_continuous_covs else None,
                categorical_covariate_keys=validated_categorical_covs if validated_categorical_covs else None
            )
            
            # Train totalVI model
            logger.info(f"Training totalVI model with {n_latent} latent dimensions")
            
            # If using GPU, get recommended GPU memory settings
            if gpu_available:
                try:
                    # Use the specified GPU device
                    gpu_memory_mb = torch.cuda.get_device_properties(gpu_device).total_memory / (1024 * 1024)
                    logger.info(f"GPU {gpu_device} memory: {gpu_memory_mb:.1f} MB")
                    
                    # Adjust batch size based on available GPU memory
                    if gpu_memory_mb > 10000:  # More than 10GB
                        batch_size = 256
                    elif gpu_memory_mb > 6000:  # More than 6GB
                        batch_size = 128
                    else:  # Less memory
                        batch_size = 64
                        
                    logger.info(f"Using batch size of {batch_size} based on available GPU memory")
                    
                    # Set CUDA_VISIBLE_DEVICES to focus on just this GPU
                    import os
                    os.environ["CUDA_VISIBLE_DEVICES"] = str(gpu_device)
                    logger.info(f"Set CUDA_VISIBLE_DEVICES={gpu_device}")
                    
                except Exception as e:
                    batch_size = 128
                    logger.info(f"Could not determine optimal batch size, using default: {batch_size}. Error: {e}")
            else:
                batch_size = 128
            
            logger.info("Training totalVI model with raw counts")
            model = scvi.model.TOTALVI(mdata, n_latent=n_latent)
            model.train(max_epochs=max_epochs, batch_size=batch_size)
            
            # Get latent representation and add to MuData
            logger.info(f"Extracting latent representation and storing as '{integrated_dim_reduction_key}'")
            latent_representation = model.get_latent_representation()
            
            # Restore original X if we swapped it
            if need_to_swap and orig_X is not None:
                logger.info("Restoring original normalized data to X matrix")
                rna_adata.X = orig_X
                # Clear reference to free memory
                orig_X = None
            
            # Transfer latent representation to the MuData object
            mdata.obsm[integrated_dim_reduction_key] = latent_representation
            
            # Store model for later use
            mdata.uns["totalVI_model"] = model
            
            # Add integration info to metadata
            mdata.uns["totalVI_info"] = {
                # GPU information
                "gpu_available": gpu_available,
                "gpu_count": gpu_count if gpu_available else 0,
                "gpu_names": gpu_names if gpu_available else [],
                "gpu_device_used": gpu_device if gpu_available else None,
                "batch_size": batch_size,
                "n_epochs": max_epochs,
                
                # Batch and covariate information
                "batch_key": validated_batch_key,
                "continuous_covariates": validated_continuous_covs,
                "categorical_covariates": validated_categorical_covs,
                
                # Raw counts information
                "raw_counts_location": raw_counts_location,
                "raw_counts_layer": raw_counts_layer if raw_counts_location == "layer" else None,
                "temporary_swap_performed": need_to_swap
            }
            
        except ImportError:
            logger.warning("scvi-tools not installed, falling back to feature concatenation")
            method = "feature-concat"
        except Exception as e:
            logger.warning(f"totalVI integration failed: {e}, falling back to feature concatenation")
            method = "feature-concat"
    
    if method == "feature-concat":
        logger.info("Using feature concatenation for integration")
        
        # Ensure PCA has been run on RNA
        if "X_pca" not in mdata.mod["rna"].obsm:
            logger.warning("PCA not found in RNA modality, running PCA")
            sc.tl.pca(mdata.mod["rna"])
        
        # Get RNA PCs and normalized protein expression
        rna_pcs = mdata.mod["rna"].obsm["X_pca"]
        
        if not isinstance(mdata.mod["prot"].X, np.ndarray):
            prot_matrix = mdata.mod["prot"].X.toarray()
        else:
            prot_matrix = mdata.mod["prot"].X
        
        # Concatenate matrices
        logger.info(f"Concatenating RNA PCs and normalized protein expression and storing as '{integrated_dim_reduction_key}'")
        integrated_matrix = np.concatenate([rna_pcs, prot_matrix], axis=1)
        mdata.obsm[integrated_dim_reduction_key] = integrated_matrix
    
    return mdata


def cluster_and_visualize(
    mdata: mu.MuData,
    n_neighbors: int = 15,
    resolution: float = 0.5,
    umap_min_dist: float = 0.3,
    integrated_dim_reduction_key: str = "X_integrated",
) -> mu.MuData:
    """
    Perform clustering and UMAP visualization on integrated data.
    
    Args:
        mdata: MuData object with integrated data
        n_neighbors: Number of neighbors for neighborhood graph
        resolution: Resolution for Leiden clustering
        umap_min_dist: Minimum distance parameter for UMAP
        integrated_dim_reduction_key: Key name of the integrated representation in obsm
        
    Returns:
        MuData object with clustering and UMAP results
    """
    # Check if integrated representation exists
    if integrated_dim_reduction_key not in mdata.obsm:
        logger.error(f"Integrated representation '{integrated_dim_reduction_key}' not found. Run integration first.")
        sys.exit(1)
    
    logger.info(f"Computing neighbors with n_neighbors={n_neighbors} using '{integrated_dim_reduction_key}'")
    sc.pp.neighbors(mdata, n_neighbors=n_neighbors, use_rep=integrated_dim_reduction_key)
    
    logger.info(f"Running UMAP with min_dist={umap_min_dist}")
    sc.tl.umap(mdata, min_dist=umap_min_dist)
    
    logger.info(f"Performing Leiden clustering with resolution={resolution}")
    sc.tl.leiden(mdata, resolution=resolution)
    
    return mdata


def save_results(
    mdata: mu.MuData, 
    output_file: str, 
    output_format: str = "mudata"
) -> None:
    """
    Save results in the specified format.
    
    Args:
        mdata: MuData object with analysis results
        output_file: Path for output file
        output_format: Output format ('mudata', 'anndata', or 'tsv')
    """
    if output_format == "mudata":
        logger.info(f"Saving MuData to {output_file}")
        mdata.write(output_file)
    
    elif output_format == "anndata":
        logger.info(f"Saving as AnnData to {output_file}")
        # Create AnnData with RNA as base and add protein data to obs
        adata = mdata.mod["rna"].copy()
        
        # Add protein expression as obs
        if isinstance(mdata.mod["prot"].X, np.ndarray):
            prot_df = pd.DataFrame(
                mdata.mod["prot"].X,
                index=mdata.mod["prot"].obs_names,
                columns=mdata.mod["prot"].var_names
            )
        else:
            prot_df = pd.DataFrame(
                mdata.mod["prot"].X.toarray(),
                index=mdata.mod["prot"].obs_names,
                columns=mdata.mod["prot"].var_names
            )
        
        # Add protein columns with prefix
        for col in prot_df.columns:
            adata.obs[f"ADT_{col}"] = prot_df[col]
        
        # Add clustering and UMAP
        if "leiden" in mdata.obs:
            adata.obs["leiden"] = mdata.obs["leiden"]
        
        if "X_umap" in mdata.obsm:
            adata.obsm["X_umap"] = mdata.obsm["X_umap"]
        
        if "X_integrated" in mdata.obsm:
            adata.obsm["X_integrated"] = mdata.obsm["X_integrated"]
        
        # Save
        adata.write(output_file)
    
    elif output_format == "tsv":
        logger.info(f"Saving metadata as TSV to {output_file}")
        # Create metadata DataFrame
        metadata = mdata.obs.copy()
        
        # Add protein expression
        if isinstance(mdata.mod["prot"].X, np.ndarray):
            prot_df = pd.DataFrame(
                mdata.mod["prot"].X,
                index=mdata.mod["prot"].obs_names,
                columns=mdata.mod["prot"].var_names
            )
        else:
            prot_df = pd.DataFrame(
                mdata.mod["prot"].X.toarray(),
                index=mdata.mod["prot"].obs_names,
                columns=mdata.mod["prot"].var_names
            )
        
        # Add protein columns with prefix
        for col in prot_df.columns:
            metadata[f"ADT_{col}"] = prot_df[col]
        
        # Add UMAP coordinates
        if "X_umap" in mdata.obsm:
            metadata["UMAP_1"] = mdata.obsm["X_umap"][:, 0]
            metadata["UMAP_2"] = mdata.obsm["X_umap"][:, 1]
        
        # Save
        metadata.to_csv(output_file, sep="\t")


def run_cite_seq_pipeline(
    # Input parameters - file paths
    rna_input: Optional[str] = None,
    adt_input: Optional[str] = None,
    mudata_input: Optional[str] = None,
    
    # Input parameters - in-memory objects
    rna_adata: Optional[ad.AnnData] = None,
    adt_adata: Optional[ad.AnnData] = None, 
    mdata: Optional[mu.MuData] = None,
    
    # QC parameters
    rna_min_genes: int = 200,
    rna_min_cells: int = 3,
    max_mito_pct: float = 20.0,
    mito_prefix: str = "MT-",
    adt_min_counts: int = 1,
    
    # Normalization parameters
    adt_norm_method: str = "clr",
    use_isotype_control: bool = False,
    isotype_controls: List[str] = None,
    rna_target_sum: float = 1e4,
    normalize_rna: bool = False,
    raw_counts_location: str = "X",
    raw_counts_layer: str = "counts",
    
    # HVG and scaling parameters
    hvg_min_mean: float = 0.0125,
    hvg_max_mean: float = 3,
    hvg_min_disp: float = 0.5,
    scale_max_value: float = 10,
    
    # Integration parameters
    integration_method: str = "feature-concat",
    n_latent: int = 20,
    max_epochs: int = 400,
    n_pcs: int = 50,
    gpu_device: int = 0,
    batch_key: Optional[str] = None,
    continuous_covariate_keys: Optional[List[str]] = None,
    categorical_covariate_keys: Optional[List[str]] = None,
    integrated_dim_reduction_key: str = "X_integrated",
    
    # Clustering parameters
    n_neighbors: int = 15,
    resolution: float = 0.5,
    umap_min_dist: float = 0.3,
    
    # Output parameters
    output_format: str = "mudata",
    output_file: str = None,
) -> mu.MuData:
    """
    Run the complete CITE-Seq analysis pipeline.
    
    This pipeline assumes that RNA data has already been QC'ed. It will perform QC
    on the ADT modality and align cells between modalities.
    
    Args:
        rna_input: Path to RNA data file (pre-QC'ed)
        adt_input: Path to ADT data file
        mudata_input: Path to MuData object file
        rna_adata: AnnData object with RNA data (in memory)
        adt_adata: AnnData object with ADT data (in memory)
        mdata: MuData object with RNA and ADT modalities (in memory)
        Various parameters for each step (see individual functions)
        
    Returns:
        MuData object with analysis results
    """
    # 1. Load data
    mdata, mdata_raw = load_data(
        rna_input=rna_input, 
        adt_input=adt_input, 
        mudata_input=mudata_input,
        rna_adata=rna_adata,
        adt_adata=adt_adata,
        mdata=mdata
    )
    
    # 2. Perform QC
    mdata = perform_qc(
        mdata=mdata,
        rna_min_genes=rna_min_genes,
        rna_min_cells=rna_min_cells,
        max_mito_pct=max_mito_pct,
        mito_prefix=mito_prefix,
        adt_min_counts=adt_min_counts
    )
    
    # 3. Normalize ADT data
    mdata = normalize_adt(
        mdata=mdata,
        mdata_raw=mdata_raw,
        method=adt_norm_method,
        use_isotype_control=use_isotype_control,
        isotype_controls=isotype_controls
    )
    
    # 4. Process RNA data
    mdata = process_rna(
        mdata=mdata,
        target_sum=rna_target_sum,  # This is ok - function expects target_sum
        hvg_min_mean=hvg_min_mean,
        hvg_max_mean=hvg_max_mean,
        hvg_min_disp=hvg_min_disp,
        scale_max_value=scale_max_value,
        n_pcs=n_pcs,
        normalize_rna=normalize_rna,
        raw_counts_location=raw_counts_location,
        raw_counts_layer=raw_counts_layer
    )
    
    # 5. Integrate modalities
    mdata = integrate_modalities(
        mdata=mdata,
        method=integration_method,
        n_latent=n_latent,
        max_epochs=max_epochs,
        gpu_device=gpu_device,
        batch_key=batch_key,
        continuous_covariate_keys=continuous_covariate_keys,
        categorical_covariate_keys=categorical_covariate_keys,
        integrated_dim_reduction_key=integrated_dim_reduction_key,
        raw_counts_location=raw_counts_location,
        raw_counts_layer=raw_counts_layer
    )
    
    # 6. Cluster and visualize
    mdata = cluster_and_visualize(
        mdata=mdata,
        n_neighbors=n_neighbors,
        resolution=resolution,
        umap_min_dist=umap_min_dist,
        integrated_dim_reduction_key=integrated_dim_reduction_key
    )
    
    # 7. Save results if output file is specified
    if output_file:
        save_results(
            mdata=mdata,
            output_file=output_file,
            output_format=output_format
        )
    
    return mdata


def main():
    """Main function to run from command line."""
    parser = setup_arg_parser()
    args = parser.parse_args()
    
    # Convert args to dict and pass to the pipeline function
    args_dict = vars(args)
    
    # When called from command line, we're only using file paths, not in-memory objects
    # So we explicitly set the in-memory object parameters to None
    args_dict['rna_adata'] = None
    args_dict['adt_adata'] = None
    args_dict['mdata'] = None
    
    # Rename parameters to match function signature
    if 'continuous_covariates' in args_dict:
        args_dict['continuous_covariate_keys'] = args_dict.pop('continuous_covariates')
    if 'categorical_covariates' in args_dict:
        args_dict['categorical_covariate_keys'] = args_dict.pop('categorical_covariates')
    
    # Run the pipeline
    run_cite_seq_pipeline(**args_dict)
    
    logger.info("CITE-Seq analysis completed successfully")


if __name__ == "__main__":
    main()
