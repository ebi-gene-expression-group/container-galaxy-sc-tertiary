#!/usr/bin/env python
"""
Script to prepare test data for CITE-seq analysis.
This script loads a 10x HDF5 file, separates RNA and ADT data, performs QC on RNA data,
and saves them as separate AnnData files for testing.
"""

import os
import sys
import logging
import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

# Set up logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger("prepare-test-data")

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Prepare test data for CITE-seq analysis"
    )
    parser.add_argument(
        "--input", "-i", 
        type=str, 
        required=True,
        help="Path to 10x HDF5 file containing both RNA and ADT data"
    )
    parser.add_argument(
        "--output-dir", "-o", 
        type=str, 
        default="./test_data",
        help="Directory to save output files (default: ./test_data)"
    )
    parser.add_argument(
        "--rna-prefix", 
        type=str, 
        default="GEX",
        help="Prefix or key for RNA data in HDF5 file (default: GEX)"
    )
    parser.add_argument(
        "--adt-prefix", 
        type=str, 
        default="ADT",
        help="Prefix or key for ADT data in HDF5 file (default: ADT)"
    )
    parser.add_argument(
        "--min-genes", 
        type=int, 
        default=200,
        help="Minimum number of genes per cell for RNA QC (default: 200)"
    )
    parser.add_argument(
        "--min-cells", 
        type=int, 
        default=3,
        help="Minimum number of cells per gene for RNA QC (default: 3)"
    )
    parser.add_argument(
        "--max-mito-pct", 
        type=float, 
        default=25.0,
        help="Maximum percentage of mitochondrial genes (default: 25.0)"
    )
    parser.add_argument(
        "--mito-prefix", 
        type=str, 
        default="MT-",
        help="Prefix for mitochondrial genes (default: MT-)"
    )
    parser.add_argument(
        "--n-hvg", 
        type=int, 
        default=2000,
        help="Number of highly variable genes to select (default: 2000)"
    )
    return parser.parse_args()

def load_and_split_10x_h5(file_path, rna_prefix="GEX", adt_prefix="ADT"):
    """
    Load a 10x HDF5 file and split into RNA and ADT AnnData objects.
    
    Args:
        file_path: Path to the 10x HDF5 file
        rna_prefix: Prefix or key for RNA data in the HDF5 file 
        adt_prefix: Prefix or key for ADT data in the HDF5 file
    
    Returns:
        Tuple of (rna_adata, adt_adata)
    """
    logger.info(f"Loading 10x HDF5 file: {file_path}")
    
    try:
        # Try to load as a mudata with modalities
        import muon as mu
        try:
            mdata = mu.read_10x_h5(file_path)
            if len(mdata.mod) >= 2:
                logger.info(f"Successfully loaded as MuData with {len(mdata.mod)} modalities")
                # Find RNA and ADT modalities by name
                if rna_prefix in mdata.mod:
                    rna_adata = mdata.mod[rna_prefix].copy()
                    logger.info(f"Found RNA modality with {rna_adata.shape[0]} cells and {rna_adata.shape[1]} genes")
                else:
                    logger.error(f"RNA modality with prefix '{rna_prefix}' not found in file")
                    sys.exit(1)
                
                if adt_prefix in mdata.mod:
                    adt_adata = mdata.mod[adt_prefix].copy()
                    logger.info(f"Found ADT modality with {adt_adata.shape[0]} cells and {adt_adata.shape[1]} features")
                else:
                    logger.error(f"ADT modality with prefix '{adt_prefix}' not found in file")
                    sys.exit(1)
                
                return rna_adata, adt_adata
        except:
            logger.info("Could not load as MuData, trying alternative methods")
            pass
            
        # Try to load as a single AnnData and split based on gene names
        adata = sc.read_10x_h5(file_path)
        logger.info(f"Loaded data with {adata.shape[0]} cells and {adata.shape[1]} features")
        
        # Check if we have a column that distinguishes between RNA and ADT
        feature_types = None
        if 'feature_types' in adata.var:
            feature_types = adata.var['feature_types']
        elif 'feature_type' in adata.var:
            feature_types = adata.var['feature_type']
        
        if feature_types is not None:
            # Split based on feature types
            rna_mask = np.array([t.lower() in ["gene", "gene expression"] for t in feature_types])
            adt_mask = np.array([t.lower() in ["antibody", "antibody capture", "adt"] for t in feature_types])
            
            if np.sum(rna_mask) > 0 and np.sum(adt_mask) > 0:
                logger.info(f"Splitting data based on feature types: {np.sum(rna_mask)} RNA genes, {np.sum(adt_mask)} ADT features")
                rna_adata = adata[:, rna_mask].copy()
                adt_adata = adata[:, adt_mask].copy()
                return rna_adata, adt_adata
        
        # If we couldn't split by feature types, try to infer by gene naming patterns
        logger.info("Trying to split data based on gene naming patterns")
        var_names = adata.var_names.tolist()
        
        # Common prefixes for antibodies
        adt_prefixes = ["CD", "HLA-", "AB-", "ADT-"]
        
        adt_mask = np.zeros(len(var_names), dtype=bool)
        for prefix in adt_prefixes:
            adt_mask = adt_mask | np.array([name.startswith(prefix) for name in var_names])
        
        # Assume everything else is RNA
        rna_mask = ~adt_mask
        
        if np.sum(rna_mask) > 0 and np.sum(adt_mask) > 0:
            logger.info(f"Splitting data based on gene naming: {np.sum(rna_mask)} RNA genes, {np.sum(adt_mask)} ADT features")
            rna_adata = adata[:, rna_mask].copy()
            adt_adata = adata[:, adt_mask].copy()
            return rna_adata, adt_adata
        
        # If we still couldn't split, raise an error
        logger.error("Could not determine how to split RNA and ADT data")
        sys.exit(1)
        
    except Exception as e:
        logger.error(f"Error loading 10x HDF5 file: {e}")
        sys.exit(1)

def qc_rna_data(adata, min_genes=200, min_cells=3, max_mito_pct=25.0, mito_prefix="MT-"):
    """
    Perform quality control on RNA data.
    
    Args:
        adata: AnnData object containing RNA data
        min_genes: Minimum number of genes per cell
        min_cells: Minimum number of cells per gene
        max_mito_pct: Maximum percentage of mitochondrial genes
        mito_prefix: Prefix for mitochondrial genes
    
    Returns:
        AnnData object with QC metrics and filtered
    """
    logger.info(f"Performing RNA QC with min_genes={min_genes}, min_cells={min_cells}, max_mito_pct={max_mito_pct}")
    
    # Calculate QC metrics
    sc.pp.calculate_qc_metrics(
        adata, 
        qc_vars=[mito_prefix], 
        inplace=True, 
        percent_top=None
    )
    
    # Filter cells
    mito_key = f"pct_{mito_prefix}_genes" if f"pct_{mito_prefix}_genes" in adata.obs else f"pct_{mito_prefix.lower()}_genes"
    if mito_key not in adata.obs:
        logger.warning(f"Mitochondrial metric '{mito_key}' not found in obs. Available keys: {list(adata.obs.columns)}")
        # Try to calculate manually
        if any(g.startswith(mito_prefix) for g in adata.var_names):
            logger.info(f"Calculating mitochondrial percentage manually using prefix '{mito_prefix}'")
            mito_genes = [g.startswith(mito_prefix) for g in adata.var_names]
            adata.obs[mito_key] = np.sum(adata[:, mito_genes].X, axis=1) / np.sum(adata.X, axis=1) * 100
        else:
            logger.warning(f"No mitochondrial genes found with prefix '{mito_prefix}'. Skipping mitochondrial filtering.")
            mito_key = None
    
    # Filter cells based on QC metrics
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)
    
    # Filter by mitochondrial percentage if available
    if mito_key is not None:
        adata = adata[adata.obs[mito_key] < max_mito_pct, :]
    
    logger.info(f"After QC: {adata.shape[0]} cells and {adata.shape[1]} genes")
    return adata

def process_rna_data(adata, n_hvg=2000):
    """
    Process RNA data: normalize, log-transform, find HVGs, scale.
    
    Args:
        adata: AnnData object containing RNA data
        n_hvg: Number of highly variable genes to select
    
    Returns:
        Processed AnnData object
    """
    logger.info(f"Processing RNA data")
    
    # Normalize to 10,000 reads per cell
    sc.pp.normalize_total(adata, target_sum=1e4)
    
    # Log-transform
    sc.pp.log1p(adata)
    
    # Find highly variable genes
    sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg)
    
    # Set raw counts
    adata.raw = adata
    
    # Scale data (zero mean and unit variance)
    adata = adata[:, adata.var.highly_variable]
    sc.pp.scale(adata, max_value=10)
    
    logger.info(f"Selected {np.sum(adata.var.highly_variable)} highly variable genes")
    return adata

def main():
    """Main function to prepare test data for CITE-seq analysis."""
    args = parse_args()
    
    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Load and split data
    rna_adata, adt_adata = load_and_split_10x_h5(
        args.input,
        rna_prefix=args.rna_prefix,
        adt_prefix=args.adt_prefix
    )
    
    # QC and process RNA data
    rna_adata = qc_rna_data(
        rna_adata,
        min_genes=args.min_genes,
        min_cells=args.min_cells,
        max_mito_pct=args.max_mito_pct,
        mito_prefix=args.mito_prefix
    )
    
    rna_adata = process_rna_data(rna_adata, n_hvg=args.n_hvg)
    
    # Ensure cells are aligned between RNA and ADT
    common_cells = np.intersect1d(rna_adata.obs_names, adt_adata.obs_names)
    logger.info(f"Found {len(common_cells)} cells common to both RNA and ADT data")
    
    if len(common_cells) == 0:
        logger.error("No common cells found between RNA and ADT data")
        sys.exit(1)
    
    # Subset both datasets to common cells
    rna_adata = rna_adata[common_cells, :].copy()
    adt_adata = adt_adata[common_cells, :].copy()
    
    # Save the data
    rna_output = os.path.join(args.output_dir, "rna_preprocessed.h5ad")
    adt_output = os.path.join(args.output_dir, "adt_raw.h5ad")
    
    logger.info(f"Saving RNA data to {rna_output}")
    rna_adata.write_h5ad(rna_output)
    
    logger.info(f"Saving ADT data to {adt_output}")
    adt_adata.write_h5ad(adt_output)
    
    logger.info("Done!")
    
    # Print instructions for using the data with cite_seq_analysis.py
    print("\nTest data preparation complete!")
    print("You can now use these files with cite_seq_analysis.py as follows:")
    print(f"\npython cite_seq_analysis.py --rna-input {rna_output} --adt-input {adt_output} [other options]")
    
if __name__ == "__main__":
    main()
