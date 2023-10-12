# import the necessary packages
import argparse

import anndata as ad
import decoupler as dc
import pandas as pd

# define arguments for the script
parser = argparse.ArgumentParser()
# add AnnData input file option
parser.add_argument("-i", "--input_anndata", help="AnnData input file", required=True)
# add network input file option
parser.add_argument("-n", "--input_network", help="network input file", required=True)
# add source option
parser.add_argument(
    "-s", "--source", help="Source column name in the network", default="source"
)
# add target
parser.add_argument(
    "-t", "--target", help="Target column name in the network", default="target"
)
# add weights option
parser.add_argument(
    "-w", "--weight", help="Weight column name in the network", default="weight"
)
# optional AnnData formatted output file option
parser.add_argument(
    "-o", "--output_anndata", help="AnnData formatted output file", default=None
)
# optional network formatted output file option
parser.add_argument(
    "-on",
    "--output_network",
    help="network formatted output files prefix",
    default=None,
)
# optional obsm key to store mlm estimate
parser.add_argument(
    "-m", "--ulm_key", help="key to store ulm estimate", default="ulm_estimate"
)
# optional obsm key to store mlm p-values
parser.add_argument(
    "-p", "--ulm_p_key", help="key to store ulm p-values", default="ulm_p_values"
)
# figure size option
parser.add_argument("-f", "--figure_size", help="figure size", default="10,10")
# path to save Activities AnnData file
parser.add_argument(
    "-a", "--activities_path", help="Path to save Activities AnnData file", default=None
)
# argument for number of marker genes
parser.add_argument(
    "-g",
    "--num_marker_genes",
    help="Number of marker genes in plotting",
    type=int,
    default=1000,
)
# argument for grouping on marker genes
parser.add_argument(
    "-G",
    "--grouping_key",
    help="Column name to group on",
    default="louvain",
)


args = parser.parse_args()
# read the input file

# check that either -o or -on is specified
if args.output_anndata is None and args.output_network is None:
    raise ValueError("Please specify either -o or -on")

# read in the network and check that source, target, and weight columns are present
network = pd.read_csv(args.input_network, index_col=0)
if (
    args.source not in network.columns
    or args.target not in network.columns
    or args.weight not in network.columns
):
    raise ValueError(
        "Source, target, and weight columns are not present in the network"
    )

# read in the AnnData input file
adata = ad.read_h5ad(args.input_anndata)

dc.run_ulm(
    mat=adata,
    net=network,
    source=args.source,
    target=args.target,
    weight=args.weight,
    verbose=True,
)

# if -o is specified, and the --mlm_key or --mlm_p_key is specified, store the mlm estimates and p-values
ulm_key = "ulm_estimate"
if args.ulm_key is not None and args.ulm_p_key is not ulm_key:
    adata.obsm[args.mlm_key] = adata.obsm[ulm_key].copy()
    # delete adata.obsm[ulm_key]
    del adata.obsm[ulm_key]
    ulm_key = args.ulm_key

ulm_pvals_key = "ulm_pvals"
if args.ulm_p_key is not None and args.ulm_p_key is not ulm_pvals_key:
    adata.obsm[args.ulm_p_key] = adata.obsm[ulm_pvals_key].copy()
    del adata.obsm[ulm_pvals_key]
    ulm_pvals_key = args.ulm_p_key

if args.output_network is not None:
    # write adata.obsm[ulm_key] and adata.obsm[ulm_pvals_key] to the output network files
    adata.obsm[ulm_key].to_csv(args.output_network + "_ulm.tsv", sep="\t")
    adata.obsm[ulm_pvals_key].to_csv(args.output_network + "_ulm_pvals.tsv", sep="\t")

# if args.activities_path is specified, generate the activities AnnData and save the AnnData object to the specified path
if args.activities_path is not None:
    acts = dc.get_acts(adata, obsm_key=ulm_key)
    adata.write_h5ad(args.activities_path)
