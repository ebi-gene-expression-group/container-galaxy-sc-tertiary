#!/usr/bin/env python
"""
Example script showing how to use the CITE-Seq analysis pipeline as a module.

This script demonstrates different ways to provide data to the pipeline,
especially using in-memory AnnData and MuData objects.
"""

import os
import sys
import logging
import scanpy as sc
import muon as mu
import anndata as ad
import numpy as np
import pandas as pd

# Add the parent directory to the path so we can import the module
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from cite_seq_analysis.cite_seq_analysis import run_cite_seq_pipeline, load_rna_data, load_adt_data

# Set up logging
logging.basicConfig(
    level=logging.INFO, 
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger("cite-seq-example")


def example_with_file_paths(rna_file, adt_file, output_file):
    """Example using file paths for RNA and ADT data."""
    logger.info("Running example with file paths")
    
    # Run the pipeline with file paths
    mdata_result = run_cite_seq_pipeline(
        rna_input=rna_file,
        adt_input=adt_file,
        output_file=output_file,
        output_format="mudata"
    )
    
    return mdata_result


def example_with_anndata_objects(rna_file, adt_file, output_file):
    """Example using in-memory AnnData objects."""
    logger.info("Running example with AnnData objects")
    
    # Load the AnnData objects
    rna_adata = sc.read_h5ad(rna_file)
    adt_adata = sc.read_h5ad(adt_file)
    
    # Run the pipeline with AnnData objects
    mdata_result = run_cite_seq_pipeline(
        rna_adata=rna_adata,
        adt_adata=adt_adata,
        output_file=output_file,
        output_format="mudata"
    )
    
    return mdata_result


def example_mixed_approach(rna_file, adt_file, output_file):
    """Example using in-memory RNA AnnData object and ADT file path."""
    logger.info("Running example with mixed approach")
    
    # Load only the RNA AnnData object
    rna_adata = sc.read_h5ad(rna_file)
    
    # Run the pipeline with RNA AnnData object and ADT file path
    mdata_result = run_cite_seq_pipeline(
        rna_adata=rna_adata,
        adt_input=adt_file,
        output_file=output_file,
        output_format="mudata"
    )
    
    return mdata_result


def example_create_objects(output_dir):
    """Example creating objects from scratch."""
    logger.info("Running example with objects created from scratch")
    
    # Create a simple RNA count matrix
    rna_counts = np.random.negative_binomial(5, 0.3, size=(100, 1000))
    rna_adata = ad.AnnData(rna_counts)
    rna_adata.var_names = [f"gene_{i}" for i in range(rna_adata.shape[1])]
    rna_adata.obs_names = [f"cell_{i}" for i in range(rna_adata.shape[0])]
    
    # Add QC metrics to show that RNA data has been pre-QC'ed
    rna_adata.obs['n_genes'] = np.random.randint(500, 2000, size=rna_adata.shape[0])
    rna_adata.obs['n_counts'] = np.random.randint(1000, 10000, size=rna_adata.shape[0])
    
    # Create a simple ADT count matrix
    adt_counts = np.random.negative_binomial(2, 0.5, size=(100, 30))
    adt_adata = ad.AnnData(adt_counts)
    adt_adata.var_names = [f"protein_{i}" for i in range(adt_adata.shape[1])]
    adt_adata.obs_names = [f"cell_{i}" for i in range(rna_adata.shape[0])]
    
    # Add isotype controls
    isotype_controls = ["isotype_ctrl_1", "isotype_ctrl_2"]
    adt_adata.var_names = pd.Index(
        [isotype_controls[i % 2] if i < 2 else f"protein_{i}" 
         for i in range(adt_adata.shape[1])]
    )
    
    # Run the pipeline with the created objects
    output_file = os.path.join(output_dir, "output_created_objects.h5mu")
    mdata_result = run_cite_seq_pipeline(
        rna_adata=rna_adata,
        adt_adata=adt_adata,
        output_file=output_file,
        output_format="mudata"
    )
    
    return mdata_result


def example_advanced_preprocessing(rna_file, adt_file, output_file):
    """Example with advanced preprocessing before running the pipeline."""
    logger.info("Running example with advanced preprocessing")
    
    # Load RNA data and perform additional preprocessing
    rna_adata = sc.read_h5ad(rna_file)
    
    # Add cell cycle scores
    sc.pp.highly_variable_genes(rna_adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    sc.tl.score_genes_cell_cycle(rna_adata, s_genes=['MCM5', 'PCNA', 'TYMS'], 
                                g2m_genes=['UBE2C', 'BIRC5', 'TPX2'])
    
    # Load ADT data and perform additional preprocessing
    adt_adata = sc.read_h5ad(adt_file)
    
    # Run the pipeline with preprocessed AnnData objects
    mdata_result = run_cite_seq_pipeline(
        rna_adata=rna_adata,
        adt_adata=adt_adata,
        output_file=output_file,
        output_format="mudata"
    )
    
    return mdata_result


def example_with_gpu_acceleration(rna_file, adt_file, output_file, gpu_device=0):
    """Example using totalVI with GPU acceleration."""
    logger.info(f"Running example with totalVI integration using GPU device {gpu_device}")
    
    # Load the AnnData objects
    rna_adata = sc.read_h5ad(rna_file)
    adt_adata = sc.read_h5ad(adt_file)
    
    # Run the pipeline with totalVI and GPU acceleration
    mdata_result = run_cite_seq_pipeline(
        rna_adata=rna_adata,
        adt_adata=adt_adata,
        integration_method="totalVI",
        gpu_device=gpu_device,
        output_file=output_file,
        output_format="mudata"
    )
    
    return mdata_result


if __name__ == "__main__":
    # Define file paths
    script_dir = os.path.dirname(os.path.abspath(__file__))
    test_data_dir = os.path.join(script_dir, "test-data")
    output_dir = os.path.join(script_dir, "example_output")
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Check if test data exists, otherwise use paths
    rna_file = os.path.join(test_data_dir, "test_rna.h5ad")
    adt_file = os.path.join(test_data_dir, "test_adt.h5ad")
    
    # If test files don't exist, try using the test_cite_seq_analysis script to create them
    if not os.path.exists(rna_file) or not os.path.exists(adt_file):
        logger.info("Test data not found. Trying to create test data...")
        try:
            from test_cite_seq_analysis import create_test_data
            rna_file, adt_file, _ = create_test_data(test_data_dir)
        except ImportError:
            logger.error("Could not import create_test_data function.")
            logger.error("Please provide paths to existing RNA and ADT data files.")
            sys.exit(1)
    
    # Example 1: Using file paths
    mdata1 = example_with_file_paths(
        rna_file, 
        adt_file, 
        os.path.join(output_dir, "example1_output.h5mu")
    )
    
    # Example 2: Using in-memory AnnData objects
    mdata2 = example_with_anndata_objects(
        rna_file, 
        adt_file, 
        os.path.join(output_dir, "example2_output.h5mu")
    )
    
    # Example 3: Using mixed approach
    mdata3 = example_mixed_approach(
        rna_file, 
        adt_file, 
        os.path.join(output_dir, "example3_output.h5mu")
    )
    
    # Example 4: Creating objects from scratch
    mdata4 = example_create_objects(output_dir)
    
    # Example 5: Using totalVI with GPU acceleration (if available)
    try:
        import torch
        if torch.cuda.is_available():
            mdata5 = example_with_gpu_acceleration(
                rna_file,
                adt_file,
                os.path.join(output_dir, "example5_output_gpu.h5mu"),
                gpu_device=0  # Use first GPU
            )
            logger.info("GPU acceleration example completed")
        else:
            logger.info("Skipping GPU acceleration example - no GPU available")
    except ImportError:
        logger.info("Skipping GPU acceleration example - PyTorch not available")
    
    logger.info("All examples completed successfully!")
