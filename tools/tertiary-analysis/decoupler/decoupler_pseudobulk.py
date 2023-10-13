import argparse

import anndata
import decoupler
import pandas as pd


def get_pseudobulk(
    adata,
    sample_col,
    groups_col,
    layer=None,
    mode="sum",
    min_cells=10,
    min_counts=1000,
    use_raw=False,
):
    """
    >>> import scanpy as sc
    >>> adata = sc.datasets.pbmc68k_reduced()
    >>> adata.X = abs(adata.X).astype(int)
    >>> pseudobulk = get_pseudobulk(adata, "bulk_labels", "louvain")
    """

    return decoupler.get_pseudobulk(
        adata,
        sample_col=sample_col,
        groups_col=groups_col,
        layer=layer,
        mode=mode,
        use_raw=use_raw,
        min_cells=min_cells,
        min_counts=min_counts,
    )


def prepend_c_to_index(index_value):
    if index_value and index_value[0].isdigit():
        return "C" + index_value
    return index_value


# write results for loading into DESeq2
def write_DESeq2_inputs(pdata, layer=None, output_dir="", factor_fields=None):
    """
    >>> import scanpy as sc
    >>> adata = sc.datasets.pbmc68k_reduced()
    >>> adata.X = abs(adata.X).astype(int)
    >>> pseudobulk = get_pseudobulk(adata, "bulk_labels", "louvain")
    >>> write_DESeq2_inputs(pseudobulk)
    """
    # add / to output_dir if is not empty or if it doesn't end with /
    if output_dir != "" and not output_dir.endswith("/"):
        output_dir = output_dir + "/"
    obs_for_deseq = pdata.obs.copy()
    # replace any index starting with digits to start with C instead.
    obs_for_deseq.rename(index=prepend_c_to_index, inplace=True)
    # avoid dash that is read as point on R colnames.
    obs_for_deseq.index = obs_for_deseq.index.str.replace("-", "_")
    obs_for_deseq.index = obs_for_deseq.index.str.replace(" ", "_")
    col_metadata_file = f"{output_dir}col_metadata.tsv"
    # write obs to a col_metadata file
    if factor_fields:
        # only output the index plus the columns in factor_fields in that order
        obs_for_deseq[factor_fields].to_csv(col_metadata_file, sep="\t", index=True)
    else:
        obs_for_deseq.to_csv(col_metadata_file, sep="\t", index=True)
    # write var to a gene_metadata file
    pdata.var.to_csv(f"{output_dir}gene_metadata.tsv", sep="\t", index=True)
    # write the counts matrix of a specified layer to file
    if layer is None:
        # write the X numpy matrix transposed to file
        df = pd.DataFrame(pdata.X.T, index=pdata.var.index, columns=obs_for_deseq.index)
    else:
        df = pd.DataFrame(
            pdata.layers[layer].T, index=pdata.var.index, columns=obs_for_deseq.index
        )
    df.to_csv(f"{output_dir}counts_matrix.tsv", sep="\t", index_label="")


def plot_pseudobulk_samples(
    pseudobulk_data,
    groupby,
    figsize=(10, 10),
    save_path=None,
):
    """
    >>> import scanpy as sc
    >>> adata = sc.datasets.pbmc68k_reduced()
    >>> adata.X = abs(adata.X).astype(int)
    >>> pseudobulk = get_pseudobulk(adata, "bulk_labels", "louvain")
    >>> plot_pseudobulk_samples(pseudobulk, groupby=["bulk_labels", "louvain"], figsize=(10, 10))
    """
    fig = decoupler.plot_psbulk_samples(
        pseudobulk_data, groupby=groupby, figsize=figsize, return_fig=True
    )
    if save_path:
        fig.savefig(f"{save_path}/pseudobulk_samples.png")
    else:
        fig.show()


def plot_filter_by_expr(
    pseudobulk_data, group, min_count=None, min_total_count=None, save_path=None
):
    """
    >>> import scanpy as sc
    >>> adata = sc.datasets.pbmc68k_reduced()
    >>> adata.X = abs(adata.X).astype(int)
    >>> pseudobulk = get_pseudobulk(adata, "bulk_labels", "louvain")
    >>> plot_filter_by_expr(pseudobulk, group="bulk_labels", min_count=10, min_total_count=200)
    """
    fig = decoupler.plot_filter_by_expr(
        pseudobulk_data,
        group=group,
        min_count=min_count,
        min_total_count=min_total_count,
        return_fig=True,
    )
    if save_path:
        fig.savefig(f"{save_path}/filter_by_expr.png")
    else:
        fig.show()


def filter_by_expr(pdata, min_count=None, min_total_count=None):
    """
    >>> import scanpy as sc
    >>> adata = sc.datasets.pbmc68k_reduced()
    >>> adata.X = abs(adata.X).astype(int)
    >>> pseudobulk = get_pseudobulk(adata, "bulk_labels", "louvain")
    >>> pdata_filt = filter_by_expr(pseudobulk, min_count=10, min_total_count=200)
    """
    genes = decoupler.filter_by_expr(
        pdata, min_count=min_count, min_total_count=min_total_count
    )
    return pdata[:, genes].copy()


def check_fields(fields, adata, obs=True, context=None):
    """
    >>> import scanpy as sc
    >>> adata = sc.datasets.pbmc68k_reduced()
    >>> check_fields(["bulk_labels", "louvain"], adata, obs=True)
    """

    legend = ""
    if context:
        legend = f", passed in {context},"
    if obs:
        if not set(fields).issubset(set(adata.obs.columns)):
            raise ValueError(
                f"Some of the following fields {legend} are not present in adata.obs: {fields}. Possible fields are: {list(set(adata.obs.columns))}"
            )
    else:
        if not set(fields).issubset(set(adata.var.columns)):
            raise ValueError(
                f"Some of the following fields {legend} are not present in adata.var: {fields}. Possible fields are: {list(set(adata.var.columns))}"
            )


def main(args):
    # Load AnnData object from file
    adata = anndata.read_h5ad(args.adata_file)

    # Merge adata.obs fields specified in args.adata_obs_fields_to_merge
    if args.adata_obs_fields_to_merge:
        # first split potential groups by ";" and iterate over them
        for group in args.adata_obs_fields_to_merge.split(";"):
            fields = group.split(",")
            check_fields(fields, adata)
            adata = merge_adata_obs_fields(fields, adata)

    check_fields([args.groupby, args.sample_key], adata)

    factor_fields = None
    if args.factor_fields:
        factor_fields = args.factor_fields.split(",")
        check_fields(factor_fields, adata)

    print(f"Using mode: {args.mode}")
    # Perform pseudobulk analysis
    pseudobulk_data = get_pseudobulk(
        adata,
        sample_col=args.sample_key,
        groups_col=args.groupby,
        layer=args.layer,
        mode=args.mode,
        use_raw=args.use_raw,
        min_cells=args.min_cells,
        min_counts=args.min_counts,
    )

    # Plot pseudobulk samples
    plot_pseudobulk_samples(
        pseudobulk_data,
        args.groupby,
        save_path=args.save_path,
        figsize=args.plot_samples_figsize,
    )

    plot_filter_by_expr(
        pseudobulk_data,
        group=args.groupby,
        min_count=args.min_counts,
        min_total_count=args.min_total_counts,
        save_path=args.save_path,
    )

    # Filter by expression if enabled
    if args.filter_expr:
        filtered_adata = filter_by_expr(
            pseudobulk_data,
            min_count=args.min_counts,
            min_total_count=args.min_total_counts,
        )

        pseudobulk_data = filtered_adata

    # Save the pseudobulk data
    if args.anndata_output_path:
        pseudobulk_data.write_h5ad(args.anndata_output_path, compression="gzip")

    write_DESeq2_inputs(
        pseudobulk_data, output_dir=args.deseq2_output_path, factor_fields=factor_fields
    )


def merge_adata_obs_fields(obs_fields_to_merge, adata):
    """
    Merge adata.obs fields specified in args.adata_obs_fields_to_merge

    Parameters
    ----------
    obs_fields_to_merge : str
        Fields in adata.obs to merge, comma separated
    adata : anndata.AnnData
        The AnnData object

    Returns
    -------
    anndata.AnnData
        The merged AnnData object

    docstring tests:
    >>> import scanpy as sc
    >>> ad = sc.datasets.pbmc68k_reduced()
    >>> ad = merge_adata_obs_fields(["bulk_labels","louvain"], ad)
    >>> ad.obs.columns
    Index(['bulk_labels', 'n_genes', 'percent_mito', 'n_counts', 'S_score',
           'G2M_score', 'phase', 'louvain', 'bulk_labels_louvain'],
          dtype='object')
    """
    field_name = "_".join(obs_fields_to_merge)
    for field in obs_fields_to_merge:
        if field not in adata.obs.columns:
            raise ValueError(f"The '{field}' column is not present in adata.obs.")
        if field_name not in adata.obs.columns:
            adata.obs[field_name] = adata.obs[field].astype(str)
        else:
            adata.obs[field_name] = (
                adata.obs[field_name] + "_" + adata.obs[field].astype(str)
            )
    return adata


if __name__ == "__main__":
    # Create argument parser
    parser = argparse.ArgumentParser(
        description="Perform pseudobulk analysis on an AnnData object"
    )

    # Add arguments
    parser.add_argument("adata_file", type=str, help="Path to the AnnData file")
    parser.add_argument(
        "-m",
        "--adata_obs_fields_to_merge",
        type=str,
        help="Fields in adata.obs to merge, comma separated. You can have more than one set of fields, separated by semi-colon ;",
    )
    parser.add_argument(
        "--groupby",
        type=str,
        required=True,
        help="The column in adata.obs that defines the groups",
    )
    parser.add_argument(
        "--sample_key",
        required=True,
        type=str,
        help="The column in adata.obs that defines the samples",
    )
    # add argument for layer
    parser.add_argument(
        "--layer",
        type=str,
        default=None,
        help="The name of the layer of the AnnData object to use",
    )
    # add argument for mode
    parser.add_argument(
        "--mode",
        type=str,
        default="sum",
        help="The mode for Decoupler pseudobulk analysis",
        choices=["sum", "mean", "median"],
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
        "--min_cells",
        type=int,
        default=10,
        help="Minimum number of cells for pseudobulk analysis",
    )
    parser.add_argument(
        "--save_path", type=str, help="Path to save the plot (optional)"
    )
    parser.add_argument(
        "--min_counts",
        type=int,
        help="Minimum count threshold for filtering by expression",
    )
    parser.add_argument(
        "--min_total_counts",
        type=int,
        help="Minimum total count threshold for filtering by expression",
    )
    parser.add_argument(
        "--anndata_output_path",
        type=str,
        help="Path to save the filtered AnnData object or pseudobulk data",
    )
    parser.add_argument(
        "--filter_expr", action="store_true", help="Enable filtering by expression"
    )
    parser.add_argument(
        "--factor_fields",
        type=str,
        help="Comma separated list of fields for the factors",
    )
    parser.add_argument(
        "--deseq2_output_path",
        type=str,
        help="Path to save the DESeq2 inputs",
        required=True,
    )
    parser.add_argument(
        "--plot_samples_figsize",
        type=int,
        default=[10, 10],
        nargs=2,
        help="Size of the samples plot as a tuple (two arguments)",
    )
    parser.add_argument("--plot_filtering_figsize", type=int, default=[10, 10], nargs=2)

    # Parse the command line arguments
    args = parser.parse_args()

    # Call the main function
    main(args)
