import argparse

import anndata
import decoupler
import numpy as np
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


def create_pseudo_replicates(adata, sample_key, num_replicates, seed=None):
    """
    Create pseudo replicates for each sample in the sample_key groups.

    Parameters
    ----------
    adata : anndata.AnnData
        The AnnData object.
    sample_key : str
        The column in adata.obs that defines the samples.
    num_replicates : int
        Number of pseudo replicates to create per sample.

    Returns
    -------
    anndata.AnnData
        The AnnData object with pseudo replicates.

    Examples
    --------
    >>> import anndata
    >>> import pandas as pd
    >>> import numpy as np
    >>> data = {
    ...     'obs': pd.DataFrame({'sample': ['A', 'A', 'B', 'B']}),
    ...     'X': np.array([[1, 0], [0, 1], [1, 1], [0, 0]])
    ... }
    >>> adata = anndata.AnnData(X=data['X'], obs=data['obs'])
    >>> adata = create_pseudo_replicates(adata, 'sample', 2)
    >>> adata.obs['sample_pseudo'].tolist()
    ['A_rep1', 'A_rep2', 'B_rep1', 'B_rep2']
    """
    if seed is not None:
        np.random.seed(seed)

    new_sample_key = f"{sample_key}_pseudo"
    adata.obs[new_sample_key] = adata.obs[sample_key].astype(str)

    for sample in adata.obs[sample_key].unique():
        sample_indices = adata.obs[
            adata.obs[sample_key] == sample].index.to_numpy()
        np.random.shuffle(sample_indices)  # Shuffle the indices to randomize
        replicate_size = int(len(sample_indices) / num_replicates)
        for i in range(num_replicates):
            start_idx = i * replicate_size
            end_idx = start_idx + replicate_size
            replicate_indices = sample_indices[start_idx:end_idx]
            adata.obs.loc[replicate_indices, new_sample_key] = (
                adata.obs.loc[replicate_indices, new_sample_key] + f"_rep{i+1}"
            )

    return adata


def prepend_c_to_index(index_value):
    if index_value and index_value[0].isdigit():
        return "C" + index_value
    return index_value


def genes_to_ignore_per_contrast_field(
    count_matrix_df,
    samples_metadata,
    sample_metadata_col_contrasts,
    min_counts_per_sample=5,
    use_cpms=False,
):
    """
    # This function calculates the genes to ignore per contrast field
    # (e.g., bulk_labels, louvain).
    # It does this by first getting the count matrix for each group,
    # then identifying genes with a count below a specified threshold.
    # The genes to ignore are those that are present in more than a specified
    # number of groups.

    >>> import pandas as pd
    >>> samples_metadata = pd.DataFrame({'sample':
    ...                                    ['S1', 'S2', 'S3',
    ...                                     'S4', 'S5', 'S6'],
    ...                                  'contrast_field':
    ...                                    ['A', 'A', 'A', 'B', 'B', 'B']})
    >>> count_matrix_df = pd.DataFrame(
    ...                       {'S1':
    ...                          [30, 1, 40, 50, 30],
    ...                        'S2':
    ...                          [40, 2, 60, 50, 80],
    ...                        'S3':
    ...                          [80, 1, 60, 50, 50],
    ...                        'S4': [1, 50, 50, 50, 2],
    ...                        'S5': [3, 40, 40, 40, 2],
    ...                        'S6': [0, 50, 50, 50, 1]})
    >>> count_matrix_df.index = ['Gene1', 'Gene2', 'Gene3', 'Gene4', 'Gene5']
    >>> df = genes_to_ignore_per_contrast_field(count_matrix_df,
    ...             samples_metadata, min_counts_per_sample=5,
    ...             sample_metadata_col_contrasts='contrast_field')
    >>> df[df['contrast_field'] == 'A'].genes_to_ignore.tolist()[0]
    'Gene2'
    >>> df[df['contrast_field'] == 'B'].genes_to_ignore.tolist()[0]
    'Gene1'
    >>> df[df['contrast_field'] == 'B'].genes_to_ignore.tolist()[1]
    'Gene5'
    """

    # Initialize a dictionary to store the genes to ignore per contrast field
    contrast_fields = []
    genes_to_ignore = []

    # Iterate over the contrast fields
    for contrast_field in samples_metadata[
        sample_metadata_col_contrasts
    ].unique():
        # Get the count matrix for the current contrast field
        count_matrix_field = count_matrix_df.loc[
            :,
            (
                samples_metadata[sample_metadata_col_contrasts]
                == contrast_field
            ).tolist(),
        ]

        # We derive min_counts from the number of samples with that
        # contrast_field value
        min_counts = count_matrix_field.shape[1] * min_counts_per_sample

        if use_cpms:
            # Convert counts to counts per million (CPM)
            count_matrix_field = (
                count_matrix_field.div(count_matrix_field.sum(axis=1), axis=0)
                * 1e6
            )
            min_counts = 1  # use 1 CPM

        # Calculate the total number of cells in the current contrast field
        # (this produces a vector of counts per gene)
        total_counts_per_gene = count_matrix_field.sum(axis=1)

        # Identify genes with a count below the specified threshold
        genes = total_counts_per_gene[
            total_counts_per_gene < min_counts
        ].index.tolist()
        if len(genes) > 0:
            # genes_to_ignore[contrast_field] = " ".join(genes)
            for gene in genes:
                genes_to_ignore.append(gene)
                contrast_fields.append(contrast_field)
    # transform gene_to_ignore to a DataFrame
    # genes_to_ignore_df = pd.DataFrame(genes_to_ignore.items(),
    #                           columns=["contrast_field", "genes_to_ignore"])
    genes_to_ignore_df = pd.DataFrame(
        {"contrast_field": contrast_fields, "genes_to_ignore": genes_to_ignore}
    )
    return genes_to_ignore_df


# write results for loading into DESeq2
def write_DESeq2_inputs(
    pdata,
    layer=None,
    output_dir="",
    factor_fields=None,
    min_counts_per_sample_marking=20,
):
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
        obs_for_deseq[factor_fields].to_csv(
            col_metadata_file, sep="\t", index=True
        )
    else:
        obs_for_deseq.to_csv(col_metadata_file, sep="\t", index=True)
    # write var to a gene_metadata file
    pdata.var.to_csv(f"{output_dir}gene_metadata.tsv", sep="\t", index=True)
    # write the counts matrix of a specified layer to file
    if layer is None:
        # write the X numpy matrix transposed to file
        df = pd.DataFrame(
            pdata.X.T, index=pdata.var.index, columns=obs_for_deseq.index
        )
    else:
        df = pd.DataFrame(
            pdata.layers[layer].T,
            index=pdata.var.index,
            columns=obs_for_deseq.index,
        )
    df.to_csv(f"{output_dir}counts_matrix.tsv", sep="\t", index_label="")

    if factor_fields:
        df_genes_ignore = genes_to_ignore_per_contrast_field(
            count_matrix_df=df,
            samples_metadata=obs_for_deseq,
            sample_metadata_col_contrasts=factor_fields[0],
            min_counts_per_sample=min_counts_per_sample_marking,
        )
        df_genes_ignore.to_csv(
            f"{output_dir}genes_to_ignore_per_contrast_field.tsv", sep="\t"
        )


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
    >>> plot_pseudobulk_samples(pseudobulk,
    ...                         groupby=["bulk_labels", "louvain"],
    ...                         figsize=(10, 10))
    """
    fig = decoupler.plot_psbulk_samples(
        pseudobulk_data, groupby=groupby, figsize=figsize, return_fig=True
    )
    if save_path:
        fig.savefig(f"{save_path}/pseudobulk_samples.png")
    else:
        fig.show()


def plot_filter_by_expr(
    pseudobulk_data,
    group,
    min_count=None,
    min_total_count=None,
    save_path=None,
):
    """
    >>> import scanpy as sc
    >>> adata = sc.datasets.pbmc68k_reduced()
    >>> adata.X = abs(adata.X).astype(int)
    >>> pseudobulk = get_pseudobulk(adata, "bulk_labels", "louvain")
    >>> plot_filter_by_expr(pseudobulk, group="bulk_labels",
    ...                     min_count=10, min_total_count=200)
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
    >>> pdata_filt = filter_by_expr(pseudobulk,
    ...                             min_count=10, min_total_count=200)
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

    # Fields that will be created during the pseudobulking process
    pseudobulk_generated_fields = ['psbulk_n_cells', 'psbulk_counts']

    # Filter out the pseudobulk-generated fields from checking
    fields_to_check = [field for field in fields
                       if field not in pseudobulk_generated_fields]

    legend = ""
    if context:
        legend = f", passed in {context},"
    if obs:
        if not set(fields_to_check).issubset(set(adata.obs.columns)):
            raise ValueError(
                f"Some of the following fields {legend} are not present \
                    in adata.obs: {fields_to_check}. \
                        Possible fields are: {list(set(adata.obs.columns))}"
            )
    else:
        if not set(fields_to_check).issubset(set(adata.var.columns)):
            raise ValueError(
                f"Some of the following fields {legend} are not present \
                    in adata.var: {fields_to_check}. \
                        Possible fields are: {list(set(adata.var.columns))}"
            )


def main(args):
    # Load AnnData object from file
    adata = anndata.read_h5ad(args.adata_file)

    # Merge adata.obs fields specified in args.adata_obs_fields_to_merge
    if args.adata_obs_fields_to_merge:
        # first split potential groups by ":" and iterate over them
        for group in args.adata_obs_fields_to_merge.split(":"):
            fields = group.split(",")
            check_fields(fields, adata)
            merge_adata_obs_fields(fields, adata)

    check_fields([args.groupby, args.sample_key], adata)

    factor_fields = None
    if args.factor_fields:
        factor_fields = args.factor_fields.split(",")
        check_fields(factor_fields, adata)

    # Create pseudo replicates if specified
    if args.num_pseudo_replicates:
        adata = create_pseudo_replicates(
            adata, args.sample_key, args.num_pseudo_replicates, seed=args.seed
        )
        args.sample_key = f"{args.sample_key}_pseudo"

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

    print("Created pseudo-bulk AnnData, checking if fields still make sense.")
    print(
        "If this fails this check, it might mean that you asked for factors \
        that are not compatible with you sample identifiers (ie. asked for \
        phase in the factors, but each sample contains more than one phase,\
        try joining fields)."
    )
    if factor_fields:
        check_fields(
            factor_fields,
            pseudobulk_data,
            context=" after creation of pseudo-bulk AnnData",
        )
    print("Factors requested are adequate for the pseudo-bulked AnnData!")

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
        pseudobulk_data.write_h5ad(
            args.anndata_output_path, compression="gzip"
        )

    write_DESeq2_inputs(
        pseudobulk_data,
        output_dir=args.deseq2_output_path,
        factor_fields=factor_fields,
        min_counts_per_sample_marking=args.min_counts_per_sample_marking,
    )

    # if contrasts file is provided, produce a file with genes that should be
    # filtered for each contrasts
    if args.contrasts_file:
        contrast_genes_df = identify_genes_to_filter_per_contrast(
            contrast_file=args.contrasts_file,
            min_perc_cells_expression=args.min_gene_exp_perc_per_cell,
            adata=adata,
            obs_field=args.groupby
        )
        contrast_genes_df.to_csv(
            f"{args.save_path}/genes_to_filter_by_contrast.tsv",
            sep="\t",
            index=False,
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
    >>> merge_adata_obs_fields(["bulk_labels","louvain"], ad)
    >>> ad.obs.columns
    Index(['bulk_labels', 'n_genes', 'percent_mito', 'n_counts', 'S_score',
           'G2M_score', 'phase', 'louvain', 'bulk_labels_louvain'],
          dtype='object')
    """
    field_name = "_".join(obs_fields_to_merge)
    for field in obs_fields_to_merge:
        if field not in adata.obs.columns:
            raise ValueError(
                f"The '{field}' column is not present in adata.obs."
            )
        if field_name not in adata.obs.columns:
            adata.obs[field_name] = adata.obs[field].astype(str)
        else:
            adata.obs[field_name] = (
                adata.obs[field_name] + "_" + adata.obs[field].astype(str)
            )


def identify_genes_to_filter_per_contrast(
    contrast_file, min_perc_cells_expression, adata, obs_field
):
    """
    Identify genes to filter per contrast based on expression percentage.
    We need those genes to be under the threshold for all conditions
    in a contrast to be identified for further filtering. If
    one condition has the gene expressed above the threshold, the gene
    becomes of interest (it can be highly up or down regulated).

    Parameters
    ----------
    contrast_file : str
        Path to the contrasts file.
    min_perc_cells_expression : float
        Minimum percentage of cells that should express a gene.
    adata: adata
        Original AnnData file
    obs_field: str
        Field in the AnnData observations where the contrasts are defined.

    Returns
    -------
    None

    Examples
    --------
    >>> import anndata
    >>> import pandas as pd
    >>> import numpy as np
    >>> import os
    >>> from io import StringIO
    >>> contrast_file = StringIO(f"contrast{os.linesep}condition1-\
condition2{os.linesep}\
2*(condition1)-condition2{os.linesep}")
    >>> min_perc_cells_expression = 30.0
    >>> data = {
    ...     'obs': pd.DataFrame({'condition': ['condition1', 'condition1',
    ...                          'condition2', 'condition2']}),
    ...     'X': np.array([[1, 0, 0, 0, 0], [0, 0, 2, 2, 0],
    ...     [0, 0, 1, 1, 0], [0, 0, 0, 2, 0]]),
    ... }
    >>> adata = anndata.AnnData(X=data['X'], obs=data['obs'])
    >>> df = identify_genes_to_filter_per_contrast(
    ...     contrast_file, min_perc_cells_expression, adata, 'condition'
    ... ) # doctest:+ELLIPSIS
    Identifying genes to filter using ...
    >>> df.head() # doctest:+ELLIPSIS
                        contrast gene
    0      condition1-condition2...
    1      condition1-condition2...
    2  2*(condition1)-condition2...
    3  2*(condition1)-condition2...
    """
    import re

    # Implement the logic to identify genes to filter per contrast
    # This is a placeholder implementation
    print(
        f"Identifying genes to filter using {contrast_file} "
        f"with min expression {min_perc_cells_expression}%"
    )
    sides_regex = re.compile(r"[\+\-\*\/\(\)\^]+")

    contrasts = pd.read_csv(contrast_file, sep="\t")
    # Iterate over each line in the contrast file
    genes_filter_for_contrast = dict()
    for contrast in contrasts.iloc[:, 0]:
        conditions = set(sides_regex.split(contrast))

        selected_conditions = []
        failed_conditions = []
        for condition in conditions:
            # remove any starting or trailing whitespaces from condition
            condition = condition.strip()
            if len(condition) == 0:
                continue
            # check if the condition is simply a number, then skip it
            if condition.isnumeric():
                continue
            if condition not in adata.obs[obs_field].unique():
                # add condition to failed_conditions
                failed_conditions.append(condition)
                continue
            # append to selected_conditions
            selected_conditions.append(condition)

        if len(failed_conditions) > 0:
            raise ValueError(
                f"Condition(s) '{failed_conditions}' "
                f"from contrast {contrast} "
                f"is/are not present in the "
                f"obs_field '{obs_field}' from the AnnData object."
                f"Possible values are: "
                f"{', '.join(adata.obs[obs_field].unique())}.")
        # we want to find the genes that are below the threshold
        # of % of cells expressed for ALL the conditions in the
        # contrast. It is enough for one of the conditions
        # of the contrast to have the genes expressed above
        # the threshold of % of cells to be of interest.
        for condition in selected_conditions:
            # check the percentage of cells that express each gene
            # Filter the AnnData object based on the obs_field value
            adata_filtered = adata[adata.obs[obs_field] == condition]
            # Calculate the percentage of cells expressing each gene
            gene_expression = (adata_filtered.X > 0).mean(axis=0) * 100
            genes_to_filter = set(adata_filtered.var[
                gene_expression.transpose() < min_perc_cells_expression
            ].index.tolist())
            # Update the genes_filter_for_contrast dictionary
            if contrast in genes_filter_for_contrast.keys():
                genes_filter_for_contrast[contrast].intersection_update(
                    genes_to_filter
                )
            else:
                genes_filter_for_contrast[contrast] = genes_to_filter

    # write the genes_filter_for_contrast to pandas dataframe of two columns:
    # contrast and gene

    # Initialize an empty list to store the expanded pairs
    expanded_pairs = []

    # Iterate over the dictionary
    for contrast, genes in genes_filter_for_contrast.items():
        for gene in genes:
            expanded_pairs.append((contrast, gene))

    # Create the DataFrame
    contrast_genes_df = pd.DataFrame(
        expanded_pairs, columns=["contrast", "gene"]
    )

    return contrast_genes_df


if __name__ == "__main__":
    # Create argument parser
    parser = argparse.ArgumentParser(
        description="Perform pseudobulk analysis on an AnnData object"
    )

    # Add arguments
    parser.add_argument(
        "adata_file", type=str, help="Path to the AnnData file"
    )
    parser.add_argument(
        "-m",
        "--adata_obs_fields_to_merge",
        type=str,
        help="Fields in adata.obs to merge, comma separated. \
            You can have more than one set of fields, \
                separated by semi-colon ;",
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
    # add argument for min percentage of cells that should express a gene
    parser.add_argument(
        "--min_gene_exp_perc_per_cell",
        type=float,
        default=50,
        help="If all the conditions of one side of a contrast express a \
            gene in less than this percentage of cells, then the genes \
            will be added to a list of genes to ignore for that contrast.\
            Requires the contrast file to be provided.",
    )
    parser.add_argument(
        "--contrasts_file",
        type=str,
        required=False,
        help="Contrasts file, a one column tsv with a header, each line \
            represents a contrast as a combination of conditions at each \
            side of a substraction.",
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
        "--min_counts_per_sample_marking",
        type=int,
        default=20,
        help="Minimum count threshold per sample for \
            marking genes to be ignored after DE",
    )
    parser.add_argument(
        "--min_total_counts",
        type=int,
        help="Minimum total count threshold for filtering by expression",
    )
    parser.add_argument(
        "--num_pseudo_replicates",
        type=int,
        choices=range(3, 1000),
        help="Number of pseudo replicates to create per sample (at least 3)",
        required=False
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Random seed for pseudo replicate sampling",
    )
    parser.add_argument(
        "--anndata_output_path",
        type=str,
        help="Path to save the filtered AnnData object or pseudobulk data",
    )
    parser.add_argument(
        "--filter_expr",
        action="store_true",
        help="Enable filtering by expression",
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
    parser.add_argument(
        "--plot_filtering_figsize", type=int, default=[10, 10], nargs=2
    )

    # Parse the command line arguments
    args = parser.parse_args()

    # Call the main function
    main(args)
