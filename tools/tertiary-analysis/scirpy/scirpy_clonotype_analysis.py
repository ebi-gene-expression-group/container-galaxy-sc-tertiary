#!/usr/bin/env python

"""
TCR Clonotype Analysis script for single-cell data

This script performs TCR/BCR repertoire analysis on single-cell data using scirpy.
It supports both command-line usage and can be imported as a module.
Data can be loaded from combined AnnData/MuData files or from separate
gene expression and VDJ data files in various formats.

Author: GitHub Copilot
"""

import argparse
import os
import logging
import warnings
import sys
from typing import Tuple, Optional, List, Union, Dict, Any

import muon as mu
import numpy as np
import pandas as pd
import scanpy as sc
import scirpy as ir
import matplotlib.pyplot as plt
from matplotlib import cm as mpl_cm
from cycler import cycler
from anndata import AnnData
from mudata import MuData

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('scirpy_clonotype_analysis')

# Suppress warnings
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)


def read_input_data(
    input_file: str = None,
    gex_file: str = None,
    vdj_file: str = None,
    vdj_format: str = None,
    sample_name: str = "sample",
    **vdj_kwargs
) -> MuData:
    """
    Read input data from file(s) and convert to MuData format if necessary.
    Supports reading combined data from a MuData file, separate GEX and VDJ data,
    or VDJ data from various supported formats.
    
    Parameters
    ----------
    input_file : str, optional
        Path to the input AnnData or MuData file that contains combined data.
    gex_file : str, optional
        Path to the gene expression AnnData file.
    vdj_file : str or list of str, optional
        Path(s) to VDJ data file(s) in various formats.
    vdj_format : str, optional
        Format of the VDJ data. Options: '10x', 'tracer', 'bracer', 'airr', 'bd_rhapsody'.
    sample_name : str, default="sample"
        Sample name to use when combining separate GEX and VDJ data.
    **vdj_kwargs : dict
        Additional keyword arguments to pass to the VDJ reader function.
    
    Returns
    -------
    MuData
        The loaded data as a MuData object.
    """
    if input_file is not None:
        # Load combined data from a single file
        logger.info(f"Loading combined data from {input_file}...")
        if input_file.endswith('.h5mu'):
            mdata = mu.read(input_file)
        elif input_file.endswith('.h5ad'):
            # Create MuData object from AnnData object
            adata = sc.read(input_file)
            # Check if it contains TCR data
            if "airr" not in adata.obsm:
                raise ValueError("Input AnnData object does not contain TCR data in adata.obsm['airr']")
            # Create MuData object with separate modalities for gene expression and TCR data
            mdata = MuData({"gex": adata})
        else:
            raise ValueError(f"Unsupported file format: {input_file}")
        
        return mdata
    
    elif gex_file is None and vdj_file is None:
        raise ValueError("Either input_file, or at least one of gex_file or vdj_file must be provided")
    
    # Handle GEX data
    if gex_file is not None:
        logger.info(f"Loading gene expression data from {gex_file}...")
        if gex_file.endswith('.h5ad'):
            adata_gex = sc.read(gex_file)
        elif gex_file.endswith('.h5'):
            # Assume 10X h5 format
            adata_gex = sc.read_10x_h5(gex_file)
            adata_gex.var_names_make_unique()
        elif gex_file.endswith('.mtx') or gex_file.endswith('.mtx.gz'):
            # Assume matrix market format
            path = os.path.dirname(gex_file)
            adata_gex = sc.read_mtx(gex_file)
            try:
                genes = pd.read_csv(os.path.join(path, 'genes.tsv'), sep='\t', header=None)[0]
                barcodes = pd.read_csv(os.path.join(path, 'barcodes.tsv'), sep='\t', header=None)[0]
                adata_gex.var_names = genes
                adata_gex.obs_names = barcodes
            except:
                logger.warning("Could not find genes.tsv or barcodes.tsv. Using auto-generated names.")
        elif gex_file.endswith('.csv') or gex_file.endswith('.csv.gz'):
            # Assume CSV format with genes in columns
            adata_gex = sc.read_csv(gex_file)
        elif gex_file.endswith('.tsv') or gex_file.endswith('.tsv.gz'):
            # Assume TSV format with genes in columns
            adata_gex = sc.read_text(gex_file, delimiter='\t')
        else:
            raise ValueError(f"Unsupported GEX file format: {gex_file}")
    else:
        # Create empty AnnData if no GEX data provided
        adata_gex = AnnData(np.zeros((0, 0)))
    
    # Handle VDJ data
    if vdj_file is not None:
        logger.info(f"Loading VDJ data from {vdj_file}...")
        
        if vdj_format == '10x':
            adata_vdj = ir.io.read_10x_vdj(vdj_file, **vdj_kwargs)
        elif vdj_format == 'tracer':
            adata_vdj = ir.io.read_tracer(vdj_file, **vdj_kwargs)
        elif vdj_format == 'bracer':
            adata_vdj = ir.io.read_bracer(vdj_file, **vdj_kwargs)
        elif vdj_format == 'airr':
            adata_vdj = ir.io.read_airr(vdj_file, **vdj_kwargs)
        elif vdj_format == 'bd_rhapsody':
            adata_vdj = ir.io.read_bd_rhapsody(vdj_file, **vdj_kwargs)
        elif vdj_format == 'dandelion':
            # You may need to provide a loaded dandelion object to vdj_file in this case
            adata_vdj = ir.io.from_dandelion(vdj_file, **vdj_kwargs)
        elif vdj_file.endswith('.h5ad'):
            adata_vdj = sc.read(vdj_file)
        else:
            raise ValueError(f"Unsupported VDJ format: {vdj_format}. Please specify a valid format.")
    else:
        # No VDJ data provided
        adata_vdj = None
    
    # Create MuData object
    if adata_vdj is not None:
        # If we have both GEX and VDJ data, merge them
        if adata_gex.shape[0] > 0:
            # Process cell barcodes to ensure they can be matched
            if len(adata_gex.obs_names.intersection(adata_vdj.obs_names)) == 0:
                # Try to harmonize barcodes if there's no direct overlap
                logger.warning("No overlapping cell barcodes found between GEX and VDJ data.")
                logger.warning("Will try to create a combined MuData object anyway.")
            
            # Create MuData
            mdata = MuData({"gex": adata_gex, "airr": adata_vdj})
            
            # Add sample name if not present
            if "sample" not in mdata.obs:
                mdata.obs["sample"] = sample_name
        else:
            # Only VDJ data available
            mdata = MuData({"airr": adata_vdj})
            if "sample" not in mdata.obs:
                mdata.obs["sample"] = sample_name
    else:
        # Only GEX data available - this is not valid for TCR analysis
        raise ValueError("No VDJ/TCR data provided. Cannot perform TCR analysis.")
    
    return mdata


def perform_tcr_analysis(
    mdata: MuData,
    output_dir: str,
    input_file_name: str = "input_data",
    min_cells: int = 2,
    filter_multichain: bool = True,
    filter_orphan_chains: bool = True,
    tcrdist_cutoff: int = 15,
    plot_format: str = "png",
    plot_dpi: int = 120,
    metadata_obskey: Optional[str] = None,
    clustering_obskey: Optional[str] = None,
    perform_epitope_analysis: bool = False,
    verbose: bool = False
) -> MuData:
    """
    Perform TCR repertoire analysis on single-cell data.

    Parameters
    ----------
    mdata : MuData
        Input MuData object with TCR data in the AIRR modality.
    output_dir : str
        Directory where output files will be saved.
    input_file_name : str, default="input_data"
        Name of the input file (used for reporting).
    min_cells : int, default=2
        Minimum number of cells for clonotype visualization.
    filter_multichain : bool, default=True
        Whether to filter out multichain cells (potential doublets).
    filter_orphan_chains : bool, default=True
        Whether to filter out cells with orphan chains.
    tcrdist_cutoff : int, default=15
        Distance cutoff for TCRdist clustering.
    plot_format : str, default="png"
        Format for output plots ("png", "pdf", "svg").
    plot_dpi : int, default=120
        DPI for output plots.
    metadata_obskey : str, optional
        Column in .obs containing metadata for visualization (e.g., "sample", "patient").
    clustering_obskey : str, optional
        Column in .obs containing cell type or cluster information.
    perform_epitope_analysis : bool, default=False
        Whether to perform epitope analysis using VDJDB.
    verbose : bool, default=False
        Whether to print verbose information.

    Returns
    -------
    MuData
        The processed data object.
    """
    # Configure verbosity
    sc.settings.verbosity = 2 if verbose else 0
    
    # Set up figure parameters
    sc.set_figure_params(figsize=(4, 4))
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Log function call
    logger.info(f"Starting TCR analysis")
    logger.info(f"Results will be saved to {output_dir}")

    # Create chain indices for receptor data
    logger.info("Creating chain indices...")
    ir.pp.index_chains(mdata)

    # Perform chain QC
    logger.info("Performing chain QC...")
    ir.tl.chain_qc(mdata)

    # Plot chain pairing information
    logger.info("Plotting chain pairing information...")
    fig = plt.figure(figsize=(6, 4))
    ir.pl.group_abundance(mdata, groupby="airr:chain_pairing", target_col=metadata_obskey or "gex:source")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"chain_pairing.{plot_format}"), dpi=plot_dpi)
    plt.close()

    # Plot receptor subtype information
    logger.info("Plotting receptor subtype information...")
    fig = plt.figure(figsize=(6, 4))
    ir.pl.group_abundance(mdata, groupby="airr:receptor_subtype", target_col=metadata_obskey or "gex:source")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"receptor_subtype.{plot_format}"), dpi=plot_dpi)
    plt.close()

    # Filter out multichain cells if requested
    if filter_multichain:
        logger.info("Filtering out multichain cells...")
        before_count = mdata.n_obs
        mu.pp.filter_obs(mdata, "airr:chain_pairing", lambda x: x != "multichain")
        after_count = mdata.n_obs
        logger.info(f"Filtered out {before_count - after_count} multichain cells")

    # Filter out cells with orphan chains if requested
    if filter_orphan_chains:
        logger.info("Filtering out cells with orphan chains...")
        before_count = mdata.n_obs
        mu.pp.filter_obs(mdata, "airr:chain_pairing", lambda x: ~np.isin(x, ["orphan VDJ", "orphan VJ"]))
        after_count = mdata.n_obs
        logger.info(f"Filtered out {before_count - after_count} cells with orphan chains")

    logger.info(f"Remaining cells after filtering: {mdata.n_obs}")

    # Compute CDR3 neighborhood graph and define clonotypes
    logger.info("Computing CDR3 neighborhood graph and defining clonotypes...")
    ir.pp.ir_dist(mdata)
    ir.tl.define_clonotypes(mdata, receptor_arms="all", dual_ir="primary_only")

    # Compute clonotype network
    logger.info("Computing clonotype network...")
    ir.tl.clonotype_network(mdata, min_cells=min_cells)

    # Plot clonotype network
    logger.info("Plotting clonotype network...")
    fig = plt.figure(figsize=(8, 8))
    ir.pl.clonotype_network(
        mdata, 
        color=metadata_obskey or "gex:source",
        base_size=20,
        label_fontsize=9,
        panel_size=(7, 7)
    )
    plt.savefig(os.path.join(output_dir, f"clonotype_network.{plot_format}"), dpi=plot_dpi)
    plt.close()

    # Re-compute CDR3 neighborhood graph using amino acid similarity and define clonotype clusters
    logger.info("Computing CDR3 neighborhood graph using amino acid similarity...")
    ir.pp.ir_dist(mdata, metric="tcrdist", sequence="aa", cutoff=tcrdist_cutoff)
    ir.tl.define_clonotype_clusters(mdata, sequence="aa", metric="tcrdist", receptor_arms="all", dual_ir="any")
    ir.tl.clonotype_network(mdata, min_cells=min_cells, sequence="aa", metric="tcrdist")

    # Plot clonotype network with amino acid similarity
    logger.info("Plotting clonotype network with amino acid similarity...")
    fig = plt.figure(figsize=(8, 8))
    ir.pl.clonotype_network(
        mdata, 
        color=metadata_obskey or "gex:patient",
        base_size=20,
        label_fontsize=9,
        panel_size=(7, 7)
    )
    plt.savefig(os.path.join(output_dir, f"clonotype_network_aa_similarity.{plot_format}"), dpi=plot_dpi)
    plt.close()

    # Compute V-gene based clonotype clusters
    logger.info("Computing V-gene based clonotype clusters...")
    ir.tl.define_clonotype_clusters(
        mdata,
        sequence="aa",
        metric="tcrdist",
        receptor_arms="all",
        dual_ir="any",
        same_v_gene=True,
        key_added="cc_aa_tcrdist_same_v"
    )

    # Identify clonotypes with different V genes
    try:
        ct_different_v = mdata.obs.groupby("airr:cc_aa_tcrdist").apply(lambda x: x["airr:cc_aa_tcrdist_same_v"].nunique() > 1)
        ct_different_v = ct_different_v[ct_different_v].index.values.tolist()
        logger.info(f"Found {len(ct_different_v)} clonotype clusters with different V genes")
    except Exception as e:
        logger.warning(f"Could not identify clonotypes with different V genes: {e}")
        ct_different_v = []

    # Analyze clonal expansion
    logger.info("Analyzing clonal expansion...")
    ir.tl.clonal_expansion(mdata)

    # Plot clonal expansion
    if clustering_obskey:
        logger.info("Plotting clonal expansion by cluster...")
        fig = plt.figure(figsize=(8, 6))
        ir.pl.clonal_expansion(mdata, target_col="clone_id", groupby=clustering_obskey, breakpoints=(1, 2, 5))
        plt.savefig(os.path.join(output_dir, f"clonal_expansion_normalized.{plot_format}"), dpi=plot_dpi)
        plt.close()
        
        fig = plt.figure(figsize=(8, 6))
        ir.pl.clonal_expansion(mdata, target_col="clone_id", groupby=clustering_obskey, breakpoints=(1, 2, 5), normalize=False)
        plt.savefig(os.path.join(output_dir, f"clonal_expansion_absolute.{plot_format}"), dpi=plot_dpi)
        plt.close()
    
    # Calculate alpha diversity
    if clustering_obskey:
        logger.info("Calculating alpha diversity...")
        fig = plt.figure(figsize=(8, 6))
        ir.pl.alpha_diversity(mdata, metric="normalized_shannon_entropy", groupby=clustering_obskey)
        plt.savefig(os.path.join(output_dir, f"alpha_diversity.{plot_format}"), dpi=plot_dpi)
        plt.close()
    
    # Plot clonotype abundance
    if clustering_obskey:
        logger.info("Plotting clonotype abundance...")
        fig = plt.figure(figsize=(8, 6))
        ir.pl.group_abundance(mdata, groupby="airr:clone_id", target_col=clustering_obskey, max_cols=10)
        plt.savefig(os.path.join(output_dir, f"clonotype_abundance.{plot_format}"), dpi=plot_dpi)
        plt.close()
    
    # Generate spectratype plots
    logger.info("Generating spectratype plots...")
    fig = plt.figure(figsize=(8, 6))
    ir.pl.spectratype(mdata, color=clustering_obskey or "gex:cluster", viztype="bar", fig_kws={"dpi": plot_dpi})
    plt.savefig(os.path.join(output_dir, f"spectratype_bar.{plot_format}"), dpi=plot_dpi)
    plt.close()
    
    fig = plt.figure(figsize=(8, 6))
    ir.pl.spectratype(
        mdata,
        color=clustering_obskey or "gex:cluster",
        viztype="curve",
        curve_layout="shifted",
        fig_kws={"dpi": plot_dpi},
        kde_kws={"kde_norm": False}
    )
    plt.savefig(os.path.join(output_dir, f"spectratype_curve.{plot_format}"), dpi=plot_dpi)
    plt.close()
    
    # Calculate repertoire overlap
    if metadata_obskey:
        logger.info("Calculating repertoire overlap...")
        df, dst, lk = ir.tl.repertoire_overlap(mdata, metadata_obskey or "gex:sample", inplace=False)
        
        # Plot repertoire overlap
        fig = plt.figure(figsize=(8, 8))
        ir.pl.repertoire_overlap(
            mdata,
            metadata_obskey or "gex:sample",
            yticklabels=True,
            xticklabels=True
        )
        plt.savefig(os.path.join(output_dir, f"repertoire_overlap.{plot_format}"), dpi=plot_dpi)
        plt.close()
    
    # Detect convergent evolution
    logger.info("Detecting convergent evolution...")
    ir.tl.clonotype_convergence(mdata, key_coarse="cc_aa_tcrdist", key_fine="clone_id")
    
    # Perform epitope analysis if requested
    if perform_epitope_analysis:
        logger.info("Performing epitope analysis using VDJDB...")
        try:
            vdjdb = ir.datasets.vdjdb()
            ir.pp.ir_dist(mdata, vdjdb, metric="identity", sequence="aa")
            ir.tl.ir_query(
                mdata,
                vdjdb,
                metric="identity",
                sequence="aa",
                receptor_arms="any",
                dual_ir="any"
            )
            ir.tl.ir_query_annotate(
                mdata,
                vdjdb,
                metric="identity",
                sequence="aa",
                include_ref_cols=["antigen.species", "antigen.gene"],
                strategy="most-frequent"
            )
            
            # Save query results
            epitope_df = ir.tl.ir_query_annotate_df(
                mdata,
                vdjdb,
                metric="identity",
                sequence="aa",
                include_ref_cols=["antigen.species", "antigen.gene"]
            )
            if not epitope_df.empty:
                epitope_df.to_csv(os.path.join(output_dir, "epitope_matches.csv"), index=True)
                logger.info(f"Saved epitope matches to {os.path.join(output_dir, 'epitope_matches.csv')}")
                
        except Exception as e:
            logger.error(f"Error during epitope analysis: {e}")
    
    # Save the processed data
    logger.info("Saving processed data...")
    mdata.write(os.path.join(output_dir, "processed_tcr_data.h5mu"))
    
    # Generate summary report
    logger.info("Generating summary report...")
    with open(os.path.join(output_dir, "tcr_analysis_summary.txt"), "w") as f:
        f.write("TCR Analysis Summary\n")
        f.write("===================\n\n")
        f.write(f"Input data: {input_file_name}\n")
        f.write(f"Output directory: {output_dir}\n\n")
        
        f.write("Data statistics:\n")
        f.write(f"- Total cells after filtering: {mdata.n_obs}\n")
        f.write(f"- Unique clonotypes: {mdata.obs['airr:clone_id'].nunique()}\n")
        f.write(f"- Unique clonotype clusters: {mdata.obs['airr:cc_aa_tcrdist'].nunique()}\n\n")
        
        receptor_subtypes = mdata.obs["airr:receptor_subtype"].value_counts().to_dict()
        f.write("Receptor subtypes:\n")
        for subtype, count in receptor_subtypes.items():
            f.write(f"- {subtype}: {count} cells\n")
        f.write("\n")
        
        chain_pairings = mdata.obs["airr:chain_pairing"].value_counts().to_dict()
        f.write("Chain pairings:\n")
        for pairing, count in chain_pairings.items():
            f.write(f"- {pairing}: {count} cells\n")
        f.write("\n")
        
        if ct_different_v:
            f.write(f"Found {len(ct_different_v)} clonotype clusters with different V genes.\n")
            f.write(f"Clonotype clusters with different V genes: {', '.join(ct_different_v)}\n\n")
    
    logger.info("TCR analysis completed successfully!")
    return mdata


def main():
    """
    Parse command line arguments and run the TCR analysis.
    """
    parser = argparse.ArgumentParser(description='Perform TCR repertoire analysis on single-cell data')
    
    # Input options - mutually exclusive groups
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument('--input', '-i', type=str,
                          help='Path to input AnnData or MuData file containing both GEX and TCR data')
    input_group.add_argument('--gex', type=str,
                          help='Path to gene expression AnnData file (.h5ad, .h5, .mtx, .csv, .tsv)')

    # VDJ data options
    vdj_group = parser.add_argument_group('VDJ data options')
    vdj_group.add_argument('--vdj', type=str,
                         help='Path to VDJ/TCR data file')
    vdj_group.add_argument('--vdj-format', type=str, choices=['10x', 'tracer', 'bracer', 'airr', 'bd_rhapsody', 'dandelion'],
                         help='Format of the VDJ data')
    vdj_group.add_argument('--sample-name', type=str, default='sample',
                         help='Sample name to use when combining GEX and VDJ data (default: sample)')
    
    # Output parameters
    parser.add_argument('--output-dir', '-o', type=str, required=True,
                        help='Directory where output files will be saved')
    
    # Analysis parameters
    parser.add_argument('--min-cells', type=int, default=2,
                        help='Minimum number of cells for clonotype visualization (default: 2)')
    parser.add_argument('--filter-multichain', action='store_true',
                        help='Filter out multichain cells (potential doublets)')
    parser.add_argument('--filter-orphan-chains', action='store_true',
                        help='Filter out cells with orphan chains')
    parser.add_argument('--tcrdist-cutoff', type=int, default=15,
                        help='Distance cutoff for TCRdist clustering (default: 15)')
    
    # Metadata parameters
    parser.add_argument('--metadata-obskey', type=str, default=None,
                        help='Column in .obs containing metadata for visualization (e.g., "sample", "patient")')
    parser.add_argument('--clustering-obskey', type=str, default=None,
                        help='Column in .obs containing cell type or cluster information')
    
    # Output parameters
    parser.add_argument('--plot-format', type=str, default='png', choices=['png', 'pdf', 'svg'],
                        help='Format for output plots (default: png)')
    parser.add_argument('--plot-dpi', type=int, default=120,
                        help='DPI for output plots (default: 120)')
    
    # Advanced options
    parser.add_argument('--perform-epitope-analysis', action='store_true',
                        help='Perform epitope analysis using VDJDB')
    parser.add_argument('--verbose', '-v', action='store_true',
                        help='Print verbose information')
    
    args = parser.parse_args()
    
    # Validate arguments
    if args.vdj is not None and args.vdj_format is None:
        parser.error("--vdj-format is required when --vdj is specified")
    
    # First read the input data
    mdata = read_input_data(
        input_file=args.input,
        gex_file=args.gex,
        vdj_file=args.vdj,
        vdj_format=args.vdj_format,
        sample_name=args.sample_name
    )
    
    # Determine input_file_name for reporting
    if args.input is not None:
        input_file_name = args.input
    elif args.gex is not None and args.vdj is not None:
        input_file_name = f"GEX: {args.gex}, VDJ: {args.vdj}"
    elif args.vdj is not None:
        input_file_name = f"VDJ: {args.vdj}"
    else:
        input_file_name = "Unknown input"
    
    # Then run the TCR analysis
    perform_tcr_analysis(
        mdata=mdata,
        output_dir=args.output_dir,
        input_file_name=input_file_name,
        min_cells=args.min_cells,
        filter_multichain=args.filter_multichain,
        filter_orphan_chains=args.filter_orphan_chains,
        tcrdist_cutoff=args.tcrdist_cutoff,
        plot_format=args.plot_format,
        plot_dpi=args.plot_dpi,
        metadata_obskey=args.metadata_obskey,
        clustering_obskey=args.clustering_obskey,
        perform_epitope_analysis=args.perform_epitope_analysis,
        verbose=args.verbose
    )


if __name__ == "__main__":
    main()
