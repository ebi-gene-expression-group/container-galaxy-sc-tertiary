# import the necessary packages
import argparse

import anndata as ad
import decoupler as dc
import pandas as pd

# define arguments for the script
parser = argparse.ArgumentParser()

# add AnnData input file option
parser.add_argument(
    "-i", "--input_anndata", help="AnnData input file", required=True
)

# add network input file option
parser.add_argument(
    "-n", "--input_network", help="Network input file", required=True
)

# output file prefix
parser.add_argument(
    "-o",
    "--output",
    help="output files prefix",
    default=None,
)

# path to save Activities AnnData file
parser.add_argument(
    "-a",
    "--activities_path",
    help="Path to save Activities AnnData file",
    default=None,
)

# Column name in net with source nodes
parser.add_argument(
    "-s",
    "--source",
    help="Column name in net with source nodes.",
    default="source",
)

# Column name in net with target nodes
parser.add_argument(
    "-t",
    "--target",
    help="Column name in net with target nodes.",
    default="target",
)

# Column name in net with weights.
parser.add_argument(
    "-w", "--weight", help="Column name in net with weights.", default="weight"
)

# add boolean argument for use_raw
parser.add_argument(
    "--use_raw",
    action="store_true",
    default=False,
    help="Whether to use the raw part of the AnnData object",
)

# add argument for min_cells
parser.add_argument(
    "--min_n",
    help="Minimum of targets per source. If less, sources are removed.",
    default=5,
    type=int,
)

# add activity inference method option
parser.add_argument(
    "-m",
    "--method",
    help="Activity inference method",
    default="mlm",
    required=True,
)
args = parser.parse_args()

# check that either -o or --output is specified
if args.output is None:
    raise ValueError("Please specify either -o or --output")

# read in the AnnData input file
adata = ad.read_h5ad(args.input_anndata)

# read in the input file network input file
network = pd.read_csv(args.input_network, sep="\t")

if (
    args.source not in network.columns
    or args.target not in network.columns
    or args.weight not in network.columns
):
    raise ValueError(
        "Source, target, and weight columns are not present in the network"
    )


print(type(args.min_n))

if args.method == "mlm":
    dc.run_mlm(
        mat=adata,
        net=network,
        source=args.source,
        target=args.target,
        weight=args.weight,
        verbose=True,
        min_n=args.min_n,
        use_raw=args.use_raw,
    )

    if args.output is not None:
        # write adata.obsm[mlm_key] and adata.obsm[mlm_pvals_key] to the
        # output network files
        combined_df = pd.concat(
            [adata.obsm["mlm_estimate"], adata.obsm["mlm_pvals"]], axis=1
        )

        # Save the combined dataframe to a file
        combined_df.to_csv(args.output + ".tsv", sep="\t")

    # if args.activities_path is specified, generate the activities AnnData
    # and save the AnnData object to the specified path
    if args.activities_path is not None:
        acts = dc.get_acts(adata, obsm_key="mlm_estimate")
        acts.write_h5ad(args.activities_path)

elif args.method == "ulm":
    dc.run_ulm(
        mat=adata,
        net=network,
        source=args.source,
        target=args.target,
        weight=args.weight,
        verbose=True,
        min_n=args.min_n,
        use_raw=args.use_raw,
    )

    if args.output is not None:
        # write adata.obsm[mlm_key] and adata.obsm[mlm_pvals_key] to the
        # output network files
        combined_df = pd.concat(
            [adata.obsm["ulm_estimate"], adata.obsm["ulm_pvals"]], axis=1
        )

        # Save the combined dataframe to a file
        combined_df.to_csv(args.output + ".tsv", sep="\t")

    # if args.activities_path is specified, generate the activities AnnData
    # and save the AnnData object to the specified path
    if args.activities_path is not None:
        acts = dc.get_acts(adata, obsm_key="ulm_estimate")
        acts.write_h5ad(args.activities_path)
