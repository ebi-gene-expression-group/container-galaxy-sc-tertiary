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
import scipy.sparse
import anndata as ad
import traceback

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
        # Check if the file exists
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"File not found: {file_path}")
            
        # Strategy 1: Try to load as a MuData with modalities using muon
        try:
            import muon as mu
            try:
                logger.info("Attempting to load with muon.read_10x_h5")
                mdata = mu.read_10x_h5(file_path)
                if len(mdata.mod) >= 2:
                    logger.info(f"Successfully loaded as MuData with {len(mdata.mod)} modalities: {list(mdata.mod.keys())}")
                    
                    # Find RNA modality (try common naming conventions)
                    rna_keys = [rna_prefix, 'rna', 'RNA', 'GEX', 'Gene Expression']
                    rna_key = next((k for k in rna_keys if k in mdata.mod), None)
                    
                    # Find ADT modality (try common naming conventions)
                    adt_keys = [adt_prefix, 'adt', 'ADT', 'prot', 'protein', 'Antibody', 'Antibody Capture']
                    adt_key = next((k for k in adt_keys if k in mdata.mod), None)
                    
                    if rna_key and adt_key:
                        rna_adata = mdata.mod[rna_key].copy()
                        adt_adata = mdata.mod[adt_key].copy()
                        logger.info(f"Found RNA modality '{rna_key}' with {rna_adata.shape[0]} cells and {rna_adata.shape[1]} genes")
                        logger.info(f"Found ADT modality '{adt_key}' with {adt_adata.shape[0]} cells and {adt_adata.shape[1]} features")
                        return rna_adata, adt_adata
                    
                    if rna_key is None:
                        logger.warning(f"RNA modality not found in MuData. Available modalities: {list(mdata.mod.keys())}")
                    if adt_key is None:
                        logger.warning(f"ADT modality not found in MuData. Available modalities: {list(mdata.mod.keys())}")
            except Exception as e:
                logger.info(f"Could not load as MuData with muon: {e}")
        except ImportError:
            logger.info("muon package not available, skipping MuData loading attempt")
            
        # Strategy 2: Try to load as a single AnnData and split based on metadata
        try:
            logger.info("Attempting to load with scanpy.read_10x_h5")
            adata = sc.read_10x_h5(file_path)
            logger.info(f"Loaded data with {adata.shape[0]} cells and {adata.shape[1]} features")
            
            # Make var_names unique if they're not already
            if not adata.var_names.is_unique:
                logger.warning("Variable names are not unique, making them unique")
                adata.var_names_make_unique()
            
            # Method 1: Check for feature type annotations in var
            feature_type_columns = ['feature_types', 'feature_type', 'type', 'Type', 'modality', 'Modality']
            
            for col in feature_type_columns:
                if col in adata.var:
                    logger.info(f"Found feature type column: {col}")
                    feature_types = adata.var[col]
                    
                    # Different naming conventions for RNA
                    rna_types = ["gene", "gene expression", "rna", "count", "Gene Expression", "RNA"]
                    rna_mask = np.array([str(t).lower() in [rt.lower() for rt in rna_types] for t in feature_types])
                    
                    # Different naming conventions for ADT
                    adt_types = ["antibody", "antibody capture", "adt", "protein", "Antibody", "Antibody Capture", "CITE", "prot"]
                    adt_mask = np.array([str(t).lower() in [at.lower() for at in adt_types] for t in feature_types])
                    
                    if np.sum(rna_mask) > 0 and np.sum(adt_mask) > 0:
                        logger.info(f"Splitting data based on feature types: {np.sum(rna_mask)} RNA genes, {np.sum(adt_mask)} ADT features")
                        rna_adata = adata[:, rna_mask].copy()
                        adt_adata = adata[:, adt_mask].copy()
                        return rna_adata, adt_adata
            
            # Method 2: Try to infer by gene naming patterns
            logger.info("Trying to split data based on gene name patterns")
            var_names = [str(name) for name in adata.var_names.tolist()]
            
            # Common patterns for antibodies and surface proteins
            adt_patterns = {
                "prefixes": ["CD", "HLA-", "AB-", "ADT-", "TotalSeq", "PD-", "PD1", "PDL1", "TIGIT", "Ki67"],
                "contains": ["antibody", "isotype", "control", "IgG"]
            }
            
            # Create a mask for ADT features
            adt_mask = np.zeros(len(var_names), dtype=bool)
            
            # Check prefixes
            for prefix in adt_patterns["prefixes"]:
                prefix_mask = np.array([name.startswith(prefix) for name in var_names])
                adt_mask = adt_mask | prefix_mask
                if np.sum(prefix_mask) > 0:
                    logger.info(f"Found {np.sum(prefix_mask)} features with ADT prefix '{prefix}'")
            
            # Check substrings
            for substring in adt_patterns["contains"]:
                substring_mask = np.array([substring.lower() in name.lower() for name in var_names])
                adt_mask = adt_mask | substring_mask
                if np.sum(substring_mask) > 0:
                    logger.info(f"Found {np.sum(substring_mask)} features containing '{substring}'")
            
            # If ADT features were found, assume the rest are RNA
            if np.sum(adt_mask) > 0:
                rna_mask = ~adt_mask
                if np.sum(rna_mask) > 0:
                    logger.info(f"Splitting data based on gene naming patterns: {np.sum(rna_mask)} RNA genes, {np.sum(adt_mask)} ADT features")
                    rna_adata = adata[:, rna_mask].copy()
                    adt_adata = adata[:, adt_mask].copy()
                    return rna_adata, adt_adata
                else:
                    logger.warning("All features classified as ADT, cannot identify RNA features")
            
            # Method 3: If we have a very large number of genes and small number of ADTs, use size heuristic
            n_features = adata.shape[1]
            if n_features > 1000:  # Typical for combined data
                # Assume 95-99% of features are RNA genes, rest are ADTs
                num_adt = min(500, int(n_features * 0.05))  # Guess: at most 5% or 500 features are ADTs
                
                # Look for features that might be ADTs based on number of UMIs per feature
                # ADTs typically have much higher counts
                if scipy.sparse.issparse(adata.X):
                    feature_means = np.array(adata.X.mean(axis=0)).flatten()
                else:
                    feature_means = np.mean(adata.X, axis=0)
                
                # Sort features by their average counts
                sorted_indices = np.argsort(-feature_means)  # Descending order
                
                # Check if there's a clear separation in the distribution of feature means
                # that might indicate ADTs vs genes
                potential_adt_indices = sorted_indices[:min(500, int(n_features * 0.1))]
                potential_adt_means = feature_means[potential_adt_indices]
                
                if len(potential_adt_means) > 10:
                    # Try to find a natural break in the distribution
                    deltas = np.diff(potential_adt_means)
                    if np.max(deltas) > np.median(deltas) * 5:  # Significant drop
                        cutoff_idx = np.argmax(deltas) + 1
                        adt_indices = potential_adt_indices[:cutoff_idx]
                        adt_mask = np.zeros(n_features, dtype=bool)
                        adt_mask[adt_indices] = True
                        rna_mask = ~adt_mask
                        
                        logger.info(f"Splitting data based on expression level heuristic: {np.sum(rna_mask)} RNA genes, {np.sum(adt_mask)} ADT features")
                        rna_adata = adata[:, rna_mask].copy()
                        adt_adata = adata[:, adt_mask].copy()
                        return rna_adata, adt_adata
            
        except Exception as e:
            logger.warning(f"Error when trying to load with scanpy: {e}")
        
        # Strategy 3: Try to load two separate files if they exist
        rna_filepath = file_path.replace('.h5', '_rna.h5').replace('.h5ad', '_rna.h5ad')
        adt_filepath = file_path.replace('.h5', '_adt.h5').replace('.h5ad', '_adt.h5ad')
        
        if os.path.exists(rna_filepath) and os.path.exists(adt_filepath):
            logger.info(f"Found separate RNA and ADT files, loading them directly")
            rna_adata = sc.read(rna_filepath)
            adt_adata = sc.read(adt_filepath)
            logger.info(f"Loaded RNA data with {rna_adata.shape[0]} cells and {rna_adata.shape[1]} genes")
            logger.info(f"Loaded ADT data with {adt_adata.shape[0]} cells and {adt_adata.shape[1]} features")
            return rna_adata, adt_adata
            
        # If we've tried everything and still couldn't split the data
        logger.error("Could not determine how to split RNA and ADT data from the input file")
        logger.error("Please ensure your file contains both RNA and ADT data with appropriate feature annotations")
        sys.exit(1)
        
    except Exception as e:
        logger.error(f"Error loading 10x HDF5 file: {e}")
        sys.exit(1)

def detect_mito_pattern(adata):
    """
    Automatically detect mitochondrial gene pattern in the dataset.
    
    Args:
        adata: AnnData object
    
    Returns:
        Tuple of (pattern_type, pattern_value) where:
        - pattern_type: 'prefix', 'substring', or None if no pattern found
        - pattern_value: The detected prefix/substring or None
    """
    common_mito_patterns = [
        # Prefixes
        ("prefix", "MT-"),   # Human
        ("prefix", "mt-"),   # Mouse
        ("prefix", "Mt-"),   # Alternative
        ("prefix", "MT."),   # Ensembl style
        ("prefix", "mt."),   # Ensembl style
        ("prefix", "mt"),    # Short form
        ("prefix", "MT"),    # Short form
        # Substrings
        ("substring", "MT"),    # Generic
        ("substring", "MITO"),  # Generic
        ("substring", "mito"),  # Generic
    ]
    
    # Check if we have gene symbols
    if adata.var_names.nlevels == 1:  # Basic index
        # Try common patterns
        for pattern_type, pattern in common_mito_patterns:
            if pattern_type == "prefix":
                matches = [g for g in adata.var_names if str(g).startswith(pattern)]
                if len(matches) > 5:  # Arbitrary threshold to avoid false positives
                    return pattern_type, pattern
            elif pattern_type == "substring":
                matches = [g for g in adata.var_names if pattern in str(g)]
                if len(matches) > 5:  # Arbitrary threshold to avoid false positives
                    return pattern_type, pattern
    
    # Check if we have a gene symbol column
    for col in ['gene_symbols', 'symbol', 'gene_symbol', 'gene_name', 'genes']:
        if col in adata.var.columns:
            for pattern_type, pattern in common_mito_patterns:
                if pattern_type == "prefix":
                    matches = [g for g in adata.var[col] if str(g).startswith(pattern)]
                    if len(matches) > 5:  # Arbitrary threshold to avoid false positives
                        return pattern_type, pattern, col
                elif pattern_type == "substring":
                    matches = [g for g in adata.var[col] if pattern in str(g)]
                    if len(matches) > 5:  # Arbitrary threshold to avoid false positives
                        return pattern_type, pattern, col
    
    # If ensembl IDs are used, check for mitochondrial chromosome
    for col in ['chromosome', 'chrom', 'chr']:
        if col in adata.var.columns:
            matches = [g for g in adata.var[col] if str(g).lower() in ['m', 'mt', 'chrm', 'chrmt']]
            if len(matches) > 5:  # Arbitrary threshold to avoid false positives
                return "chromosome", col
    
    # If gene IDs contain the ensembl mitochondrial pattern (e.g., ENSMUSG00000064336)
    if any('MT' in str(g) for g in adata.var_names):
        return "custom", 'MT'
    
    return None, None

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
    
    # Filter cells and genes first (basic QC)
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)
    
    logger.info(f"After basic filtering: {adata.shape[0]} cells and {adata.shape[1]} genes")
    
    # Try to identify mitochondrial genes
    mito_key = None
    use_mito_filtering = True
    
    try:
        # Auto-detect mitochondrial pattern
        logger.info("Attempting to auto-detect mitochondrial gene pattern")
        pattern_type, pattern_value, *extra_args = detect_mito_pattern(adata) + [None]
        
        if pattern_type == "prefix":
            logger.info(f"Detected mitochondrial genes with prefix: '{pattern_value}'")
            mito_prefix = pattern_value
            qc_vars = [pattern_value]
            mito_key = f"pct_{pattern_value}_genes".replace("-", "_").replace(".", "_").lower()
            
        elif pattern_type == "substring":
            logger.info(f"Detected mitochondrial genes with substring: '{pattern_value}'")
            # For substrings, we need to create a mask
            mito_mask = np.array([pattern_value in str(g) for g in adata.var_names])
            adata.var['mito'] = mito_mask
            qc_vars = ['mito']
            mito_key = f"pct_mito_genes"
            
        elif pattern_type == "chromosome" and pattern_value:
            logger.info(f"Detected mitochondrial genes on chromosome column: '{pattern_value}'")
            # For chromosome annotations
            col = pattern_value
            adata.var['mito'] = [str(g).lower() in ['m', 'mt', 'chrm', 'chrmt'] for g in adata.var[col]]
            qc_vars = ['mito']
            mito_key = f"pct_mito_genes"
            
        elif pattern_type == "custom":
            logger.info(f"Using custom mitochondrial detection with pattern: '{pattern_value}'")
            # Custom detection - create mask
            mito_mask = np.array([pattern_value in str(g) for g in adata.var_names])
            adata.var['mito'] = mito_mask
            qc_vars = ['mito']
            mito_key = f"pct_mito_genes"
            
        else:
            # If no pattern detected, try user-provided prefix
            logger.info(f"No mitochondrial pattern auto-detected, trying user-provided prefix: '{mito_prefix}'")
            if any(str(g).startswith(mito_prefix) for g in adata.var_names):
                logger.info(f"Found mitochondrial genes with prefix: '{mito_prefix}'")
                qc_vars = [mito_prefix]
                mito_key = f"pct_{mito_prefix}_genes".replace("-", "_").replace(".", "_").lower()
            else:
                logger.warning(f"No mitochondrial genes found with prefix '{mito_prefix}'")
                qc_vars = None
                use_mito_filtering = False
    
        # Calculate QC metrics
        if qc_vars:
            try:
                sc.pp.calculate_qc_metrics(
                    adata, 
                    qc_vars=qc_vars, 
                    inplace=True, 
                    percent_top=None
                )
                logger.info(f"Calculated QC metrics with qc_vars={qc_vars}")
            except Exception as e:
                logger.error(f"Error calculating QC metrics: {e}")
                use_mito_filtering = False
        else:
            logger.info("Calculating QC metrics without mitochondrial genes")
            sc.pp.calculate_qc_metrics(
                adata,
                inplace=True,
                percent_top=None
            )
            use_mito_filtering = False
    
    except Exception as e:
        logger.error(f"Error during mitochondrial detection: {e}")
        logger.info("Calculating basic QC metrics without mitochondrial genes")
        sc.pp.calculate_qc_metrics(
            adata,
            inplace=True,
            percent_top=None
        )
        use_mito_filtering = False
    
    # Filter by mitochondrial percentage if applicable
    if use_mito_filtering and mito_key and mito_key in adata.obs:
        logger.info(f"Filtering cells with > {max_mito_pct}% mitochondrial content using key: '{mito_key}'")
        adata = adata[adata.obs[mito_key] < max_mito_pct, :].copy()
        logger.info(f"After mitochondrial filtering: {adata.shape[0]} cells and {adata.shape[1]} genes")
    else:
        logger.warning("Skipping mitochondrial filtering due to no detected mitochondrial genes or metrics")
    
    logger.info(f"Final RNA dataset after QC: {adata.shape[0]} cells and {adata.shape[1]} genes")
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
    
    try:
        # Store original raw counts before any processing in a layer
        logger.info("Storing original raw counts in 'counts' layer")
        if scipy.sparse.issparse(adata.X):
            adata.layers['counts'] = adata.X.copy()
        else:
            adata.layers['counts'] = np.array(adata.X)
        
        # Normalize to 10,000 reads per cell
        logger.info("Normalizing data to 10,000 reads per cell")
        sc.pp.normalize_total(adata, target_sum=1e4)
        
        # Log-transform
        logger.info("Performing log1p transformation")
        sc.pp.log1p(adata)
        
        # Find highly variable genes - handle errors if they occur
        try:
            logger.info(f"Identifying {n_hvg} highly variable genes")
            sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg)
            
            # Check if HVGs were identified correctly
            if 'highly_variable' in adata.var:
                n_hvgs = np.sum(adata.var.highly_variable)
                if n_hvgs > 0:
                    logger.info(f"Selected {n_hvgs} highly variable genes")
                    adata_hvg = adata[:, adata.var.highly_variable].copy()
                    
                    # Scale data (zero mean and unit variance)
                    logger.info("Scaling data to zero mean and unit variance")
                    try:
                        sc.pp.scale(adata_hvg, max_value=10)
                        return adata_hvg
                    except Exception as e:
                        logger.warning(f"Error during scaling: {e}")
                        logger.info("Proceeding with unscaled data")
                        return adata_hvg
                else:
                    logger.warning(f"No highly variable genes identified, using all genes")
                    return adata
            else:
                logger.warning("Highly variable genes not computed properly, using all genes")
                return adata
        except Exception as e:
            logger.warning(f"Error computing highly variable genes: {e}")
            logger.info("Proceeding with all genes")
            return adata
            
    except Exception as e:
        logger.error(f"Error during RNA processing: {e}")
        logger.info("Returning minimally processed data")
        # Return the original data if processing fails
        return adata

def ensure_anndata_integrity(adata, dataset_name="dataset"):
    """
    Ensure AnnData object has the necessary properties and structure.
    
    Args:
        adata: AnnData object
        dataset_name: Name of dataset for logging
    
    Returns:
        Cleaned AnnData object
    """
    logger.info(f"Checking integrity of {dataset_name} AnnData")
    
    # Ensure var_names and obs_names are unique
    if not adata.var_names.is_unique:
        logger.warning(f"{dataset_name} var_names are not unique, making them unique")
        adata.var_names_make_unique()
    
    if not adata.obs_names.is_unique:
        logger.warning(f"{dataset_name} obs_names are not unique, making them unique")
        adata.obs_names_make_unique()
    
    # Clean up problematic columns that might cause issues later
    for attr in [adata.obs, adata.var]:
        cols_to_remove = []
        for col in attr.columns:
            # Check for columns with None values or mixed types that could cause issues
            if attr[col].isna().all():
                cols_to_remove.append(col)
            elif attr[col].dtype == 'object' and not isinstance(attr[col].iloc[0], str):
                try:
                    # Try to convert to string
                    attr[col] = attr[col].astype(str)
                except:
                    cols_to_remove.append(col)
        
        if cols_to_remove:
            logger.warning(f"Removing problematic columns from {dataset_name}: {cols_to_remove}")
            for col in cols_to_remove:
                del attr[col]
    
    return adata

def generate_synthetic_test_data(n_cells=1000, n_genes=2000, n_adt=50):
    """
    Generate synthetic test data if loading fails.
    
    Args:
        n_cells: Number of cells to generate
        n_genes: Number of RNA genes
        n_adt: Number of ADT features
    
    Returns:
        Tuple of (rna_adata, adt_adata)
    """
    logger.info(f"Generating synthetic test data with {n_cells} cells, {n_genes} genes, and {n_adt} ADT features")
    
    # Create RNA data
    rna_counts = np.random.negative_binomial(5, 0.3, size=(n_cells, n_genes))
    gene_names = [f"Gene_{i}" for i in range(n_genes)]
    
    # Add some mitochondrial genes
    mt_indices = np.random.choice(n_genes, 20, replace=False)
    for i in mt_indices:
        gene_names[i] = f"MT-{gene_names[i]}"
    
    cell_names = [f"Cell_{i}" for i in range(n_cells)]
    
    rna_adata = ad.AnnData(
        X=rna_counts,
        obs=pd.DataFrame(index=cell_names),
        var=pd.DataFrame(index=gene_names)
    )
    
    # Create ADT data
    adt_counts = np.random.negative_binomial(10, 0.2, size=(n_cells, n_adt))
    adt_names = [f"CD{i}" if i % 3 else f"HLA-DR{i}" for i in range(n_adt)]
    
    adt_adata = ad.AnnData(
        X=adt_counts,
        obs=pd.DataFrame(index=cell_names),
        var=pd.DataFrame(index=adt_names)
    )
    
    logger.info("Successfully generated synthetic test data")
    return rna_adata, adt_adata

def main():
    """Main function to prepare test data for CITE-seq analysis."""
    args = parse_args()
    
    try:
        # Create output directory if it doesn't exist
        os.makedirs(args.output_dir, exist_ok=True)
        
        # Load and split data
        try:
            rna_adata, adt_adata = load_and_split_10x_h5(
                args.input,
                rna_prefix=args.rna_prefix,
                adt_prefix=args.adt_prefix
            )
        except Exception as e:
            logger.error(f"Failed to load data: {e}")
            logger.info("Creating synthetic test data as fallback")
            rna_adata, adt_adata = generate_synthetic_test_data()
            
        # Ensure AnnData integrity before processing
        rna_adata = ensure_anndata_integrity(rna_adata, "RNA")
        adt_adata = ensure_anndata_integrity(adt_adata, "ADT")
        
        # QC and process RNA data
        try:
            rna_adata = qc_rna_data(
                rna_adata,
                min_genes=args.min_genes,
                min_cells=args.min_cells,
                max_mito_pct=args.max_mito_pct,
                mito_prefix=args.mito_prefix
            )
        except Exception as e:
            logger.error(f"Error during RNA QC: {e}")
            logger.info("Proceeding with minimal RNA QC")
            # Basic filtering as fallback
            sc.pp.filter_cells(rna_adata, min_genes=args.min_genes)
            sc.pp.filter_genes(rna_adata, min_cells=args.min_cells)
        
        try:
            rna_adata = process_rna_data(rna_adata, n_hvg=args.n_hvg)
        except Exception as e:
            logger.error(f"Error during RNA processing: {e}")
            logger.info("Using minimally processed RNA data")
        
        # Ensure cells are aligned between RNA and ADT
        common_cells = np.intersect1d(rna_adata.obs_names, adt_adata.obs_names)
        logger.info(f"Found {len(common_cells)} cells common to both RNA and ADT data")
        
        if len(common_cells) == 0:
            logger.error("No common cells found between RNA and ADT data")
            logger.info("Creating new synthetic dataset with aligned cells")
            rna_adata, adt_adata = generate_synthetic_test_data()
            common_cells = rna_adata.obs_names.tolist()
        
        # Subset both datasets to common cells
        rna_adata = rna_adata[common_cells, :].copy()
        adt_adata = adt_adata[common_cells, :].copy()
        
        # One final check to ensure everything is in order
        rna_adata = ensure_anndata_integrity(rna_adata, "Final RNA")
        adt_adata = ensure_anndata_integrity(adt_adata, "Final ADT")
        
        # Add some basic cell metadata to RNA and ADT data if not present
        if 'n_genes' not in rna_adata.obs:
            rna_adata.obs['n_genes'] = (rna_adata.X > 0).sum(axis=1)
        if 'n_counts' not in rna_adata.obs:
            rna_adata.obs['n_counts'] = rna_adata.X.sum(axis=1)
        if 'n_counts' not in adt_adata.obs:
            adt_adata.obs['n_counts'] = adt_adata.X.sum(axis=1)
        
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
        print(f"RNA data: {rna_adata.shape[0]} cells, {rna_adata.shape[1]} genes")
        print(f"ADT data: {adt_adata.shape[0]} cells, {adt_adata.shape[1]} features")
        print("\nNote: Raw RNA counts are stored in the 'counts' layer")
        print("\nYou can now use these files with cite_seq_analysis.py as follows:")
        print(f"\npython cite_seq_analysis.py \\\n    --rna-input {rna_output} \\\n    --adt-input {adt_output} \\\n    --raw-counts-location layer \\\n    --raw-counts-layer counts \\\n    [other options]")
    
    except Exception as e:
        logger.error(f"Unhandled error in main function: {e}")
        import traceback
        logger.error(traceback.format_exc())
        sys.exit(1)
    
if __name__ == "__main__":
    main()
