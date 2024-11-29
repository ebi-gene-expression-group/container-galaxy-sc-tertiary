# import the necessary packages
import argparse

import anndata as ad
import decoupler as dc
import pandas as pd

# define arguments for the script
parser = argparse.ArgumentParser()

# add AnnData input file option
parser.add_argument(
    "-i", "--input",
    help=(
        "AnnData or table input file. Table input is meant for a single "
        "comparison, not gene x cells"
    ),
    required=True
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

# Column for stat to use when providing a table
parser.add_argument(
    "--stat",
    help="Stat to use when providing a table. Default is 'log2FC'.",
    default="log2FC",
)
# Optional column for p-value or FDR in the table
parser.add_argument(
    "--p_value_column",
    required=False,
    help="Column name in the table with p-values or FDRs.",
    default=None,
)
# Optional column for FDR threshold when given a table
parser.add_argument(
    "--p_value_threshold",
    required=False,
    type=float,
    help="Column name in the table with FDRs.",
    default=0.05
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

# add activity inference method option
parser.add_argument(
    "-g",
    "--var_gene_symbols_field",
    help="Gene symbols field",
    default=None,
)
args = parser.parse_args()

# check that either -o or --output is specified
if args.output is None:
    raise ValueError("Please specify either -o or --output")

# detect based on input file extension if the input file is AnnData or matrix
if args.input.endswith(".h5ad"):
    input_type = "AnnData"
elif args.input.endswith(".tsv") or args.input.endswith(".csv"):
    input_type = "matrix"
else:
    raise ValueError(
        "Invalid input file. Please provide a valid AnnData or matrix file."
    )

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

if input_type == "AnnData":
    # read in the AnnData input file
    adata = ad.read_h5ad(args.input)

    if args.var_gene_symbols_field and args.var_gene_symbols_field in adata.var.columns:
        # Storing index in a column called 'index_bak'
        adata.var['index_bak'] = adata.var.index
        adata.var.set_index(args.var_gene_symbols_field, inplace=True)
else:
    # read in the matrix input file, genes in rows and columns for stats
    adata = pd.read_csv(args.input, sep="\t", index_col=0)
    if args.stat not in adata.columns:
        raise ValueError(f"Stat column {args.stat} not found in input table header.")
    if args.p_value_column and args.p_value_column not in adata.columns:
        raise ValueError(
            f"P-value column {args.p_value_column} not found in input table header."
        )
    if args.p_value_column and args.p_value_threshold is not None:
        adata = adata[adata[args.p_value_column] <= args.p_value_threshold]

    adata = adata[[args.stat]].T

if args.method == "mlm":
    res = dc.run_mlm(
        mat=adata,
        net=network,
        source=args.source,
        target=args.target,
        weight=args.weight,
        verbose=True,
        min_n=args.min_n,
        use_raw=args.use_raw,
    )

elif args.method == "ulm":
    res = dc.run_ulm(
        mat=adata,
        net=network,
        source=args.source,
        target=args.target,
        weight=args.weight,
        verbose=True,
        min_n=args.min_n,
        use_raw=args.use_raw,
    )
elif args.method == "consensus":
    res = dc.run_consensus(
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
    if input_type == "AnnData":
        combined_df = pd.concat(
            [adata.obsm[f"{args.method}_estimate"],
             adata.obsm[f"{args.method}_pvals"]], axis=1
        )
    else:
        tf_est, tf_pvals = res
        combined_df = pd.DataFrame(
            {
                # index is written, so no need for the set names
                f"{args.method}_estimate": tf_est.iloc[0],
                f"{args.method}_pvals": tf_pvals.iloc[0],
            }
        )
    # sort ascending on the p-values
    combined_df.sort_values(by=f"{args.method}_pvals", inplace=True)

    # Save the combined dataframe to a file
    combined_df.to_csv(args.output + ".tsv", sep="\t")

# if args.activities_path is specified, generate the activities AnnData
# and save the AnnData object to the specified path
if args.activities_path is not None and input_type == "AnnData":
    acts = dc.get_acts(adata, obsm_key=f"{args.method}_estimate")
    acts.write_h5ad(args.activities_path)
