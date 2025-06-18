#!/usr/bin/env python
"""
Test script for CITE-Seq analysis pipeline.

This script demonstrates how to use the CITE-Seq analysis pipeline
as an imported module, rather than a command-line tool.
"""

import os
import sys
import logging
import anndata as ad
import muon as mu
import numpy as np
import pandas as pd
import scanpy as sc
import argparse

# Add the parent directory to the path so we can import the module
sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
from cite_seq_analysis import run_cite_seq_pipeline

# Set up logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger("cite-seq-test")


def create_test_data(output_dir):
    """Create test data for CITE-Seq analysis."""
    # Create a simple RNA count matrix that's already QC'ed
    rna_counts = np.random.negative_binomial(5, 0.3, size=(100, 1000))
    rna_adata = ad.AnnData(rna_counts)
    rna_adata.var_names = [f"gene_{i}" for i in range(rna_adata.shape[1])]
    rna_adata.obs_names = [f"cell_{i}" for i in range(rna_adata.shape[0])]
    
    # Add some genes with MT prefix
    rna_adata.var_names = pd.Index(
        [f"MT-gene_{i}" if i < 50 else f"gene_{i}" for i in range(rna_adata.shape[1])]
    )
    
    # Add QC metrics to show that RNA data has been pre-QC'ed
    rna_adata.obs['n_genes'] = np.random.randint(500, 2000, size=rna_adata.shape[0])
    rna_adata.obs['n_counts'] = np.random.randint(1000, 10000, size=rna_adata.shape[0])
    rna_adata.obs['pct_counts_mito'] = np.random.uniform(0, 10, size=rna_adata.shape[0])
    rna_adata.uns['qc_done'] = True

    # Create a simple ADT count matrix
    adt_counts = np.random.negative_binomial(2, 0.5, size=(100, 30))
    adt_adata = ad.AnnData(adt_counts)
    adt_adata.var_names = [f"protein_{i}" for i in range(adt_adata.shape[1])]
    adt_adata.obs_names = [f"cell_{i}" for i in range(adt_adata.shape[0])]
    
    # Add isotype controls
    isotype_controls = ["isotype_ctrl_1", "isotype_ctrl_2"]
    adt_adata.var_names = pd.Index(
        [isotype_controls[i % 2] if i < 2 else f"protein_{i}" 
         for i in range(adt_adata.shape[1])]
    )

    # Save the test data
    rna_file = os.path.join(output_dir, "test_rna.h5ad")
    adt_file = os.path.join(output_dir, "test_adt.h5ad")
    
    logger.info(f"Saving test RNA data to {rna_file} (pre-QC'ed)")
    rna_adata.write(rna_file)
    
    logger.info(f"Saving test ADT data to {adt_file}")
    adt_adata.write(adt_file)
    
    # Create MuData object
    mdata = mu.MuData({"rna": rna_adata, "prot": adt_adata})
    mudata_file = os.path.join(output_dir, "test_mudata.h5mu")
    
    logger.info(f"Saving test MuData to {mudata_file}")
    mdata.write(mudata_file)
    
    return rna_file, adt_file, mudata_file


def run_test_pipeline(rna_file, adt_file, mudata_file, output_dir):
    """Run the CITE-Seq pipeline on test data."""
    # Test 1: Using separate RNA and ADT file inputs
    logger.info("Test 1: Using separate RNA and ADT file inputs")
    output_file_separate = os.path.join(output_dir, "output_separate.h5mu")
    mdata_separate = run_cite_seq_pipeline(
        rna_input=rna_file,
        adt_input=adt_file,
        output_file=output_file_separate,
        output_format="mudata"
    )

    # Test 2: Using MuData file input
    logger.info("Test 2: Using MuData file input")
    output_file_mudata = os.path.join(output_dir, "output_mudata.h5mu")
    mdata_combined = run_cite_seq_pipeline(
        mudata_input=mudata_file,
        output_file=output_file_mudata,
        output_format="mudata"
    )

    # Test 3: Using in-memory AnnData objects
    logger.info("Test 3: Using in-memory AnnData objects")
    rna_adata = sc.read_h5ad(rna_file)
    adt_adata = sc.read_h5ad(adt_file)
    output_file_inmem = os.path.join(output_dir, "output_inmem.h5mu")
    mdata_inmem = run_cite_seq_pipeline(
        rna_adata=rna_adata,
        adt_adata=adt_adata,
        output_file=output_file_inmem,
        output_format="mudata"
    )
    
    # Test 4: Using in-memory MuData object
    logger.info("Test 4: Using in-memory MuData object")
    mdata_obj = mu.read(mudata_file)
    output_file_mdata = os.path.join(output_dir, "output_mdata_obj.h5mu")
    mdata_result = run_cite_seq_pipeline(
        mdata=mdata_obj,
        output_file=output_file_mdata,
        output_format="mudata"
    )
    
    # Test 5: Using in-memory RNA AnnData with ADT file
    logger.info("Test 5: Using in-memory RNA AnnData with ADT file")
    output_file_mixed = os.path.join(output_dir, "output_mixed.h5mu")
    mdata_mixed = run_cite_seq_pipeline(
        rna_adata=rna_adata,
        adt_input=adt_file,
        output_file=output_file_mixed,
        output_format="mudata"
    )

    # Test 6: Feature-concat integration and TSV output
    logger.info("Test 6: Feature-concat integration and TSV output")
    output_file_tsv = os.path.join(output_dir, "output_metadata.tsv")
    run_cite_seq_pipeline(
        rna_input=rna_file,
        adt_input=adt_file,
        integration_method="feature-concat",
        output_file=output_file_tsv,
        output_format="tsv"
    )

    # Test 7: AnnData output
    logger.info("Test 7: AnnData output")
    output_file_anndata = os.path.join(output_dir, "output_anndata.h5ad")
    run_cite_seq_pipeline(
        rna_input=rna_file,
        adt_input=adt_file,
        output_file=output_file_anndata,
        output_format="anndata"
    )

    logger.info("All tests completed successfully!")


def main():
    """Main function to run the test script."""
    parser = argparse.ArgumentParser(description="Test CITE-Seq analysis pipeline")
    parser.add_argument(
        "--output-dir", 
        type=str, 
        default="./test_output", 
        help="Output directory for test files"
    )
    parser.add_argument(
        "--create-data-only", 
        action="store_true", 
        help="Only create test data, don't run pipeline"
    )
    
    args = parser.parse_args()
    
    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Create test data
    rna_file, adt_file, mudata_file = create_test_data(args.output_dir)
    
    if not args.create_data_only:
        # Run pipeline on test data
        run_test_pipeline(rna_file, adt_file, mudata_file, args.output_dir)


if __name__ == "__main__":
    main()
