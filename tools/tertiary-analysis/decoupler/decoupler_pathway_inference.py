# import the necessary packages
import argparse

import anndata as ad
import decoupler as dc
import pandas as pd
import hdf5plugin

# define arguments for the script
parser = argparse.ArgumentParser()

# add AnnData input file option
parser.add_argument("-i", "--input_anndata", help="AnnData input file", required=True)

# add network input file option
parser.add_argument("--species", help="Species name to download network. Default: Human", default="human")

# output file prefix
parser.add_argument(
    "--output",
    help="output files prefix",
    default=None,
)

# optional obsm key to store mlm estimate
parser.add_argument(
    "-m", "--mlm_key", help="key to store mlm estimate", default="new_mlm_estimate"
)

# optional obsm key to store mlm p-values
parser.add_argument(
    "-p", "--mlm_p_key", help="key to store mlm p-values", default="new_mlm_p_values"
)

# figure size option
parser.add_argument("-f", "--figure_size", help="figure size", default="10,10")

# path to save Activities AnnData file
parser.add_argument(
    "-a", "--activities_path", help="Path to save Activities AnnData file", default=None
)

# Column name in net with source nodes
parser.add_argument(
    "-s", "--source", help="Column name in net with source nodes.", default="source"
)

# Column name in net with target nodes
parser.add_argument(
    "-t", "--target", help="Column name in net with target nodes.", default="target"
)

# Column name in net with weights.
parser.add_argument(
    "-w", "--weight", help="Column name in net with weights.", default="weight"
)


args = parser.parse_args()

# check that either -o or --output is specified
if args.output is None:
    raise ValueError("Please specify either -o or --output")

# read in the AnnData input file
adata = ad.read_h5ad(args.input_anndata)

# Downloads the network file from prpgeny
network =  dc.get_progeny(organism=args.species, top=500)

if (
    args.source not in network.columns
    or args.target not in network.columns
    or args.weight not in network.columns
):
    raise ValueError(
        "Source, target, and weight columns are not present in the network"
    )



dc.run_mlm(
    mat=adata,
    net=network,
    source=args.source,
    target=args.target,
    weight=args.weight,
    verbose=True
)

mlm_key = "mlm_estimate"
if args.mlm_key is not None and args.mlm_p_key is not mlm_key:
    adata.obsm[args.mlm_key] = adata.obsm[mlm_key].copy()
    # delete adata.obsm[mlm_key]
    del adata.obsm[mlm_key]
    mlm_key = args.mlm_key

mlm_pvals_key = "mlm_pvals"
if args.mlm_p_key is not None and args.mlm_p_key is not mlm_pvals_key:
    adata.obsm[args.mlm_p_key] = adata.obsm[mlm_pvals_key].copy()
    # delete adata.obsm[mlm_pvals_key]
    del adata.obsm[mlm_pvals_key]
    mlm_pvals_key = args.mlm_p_key

if args.output is not None:
    # write adata.obsm[ulm_key] and adata.obsm[ulm_pvals_key] to the output network files
    adata.obsm[mlm_key].to_csv(args.output + "_mlm.tsv", sep="\t")
    adata.obsm[mlm_pvals_key].to_csv(args.output + "_mlm_pvals.tsv", sep="\t")

# if args.activities_path is specified, generate the activities AnnData and save the AnnData object to the specified path
if args.activities_path is not None:
    acts = dc.get_acts(adata, obsm_key=mlm_key)
    adata.write_h5ad(args.activities_path)