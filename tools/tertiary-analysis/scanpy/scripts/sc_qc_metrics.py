import argparse

import matplotlib.pyplot as plt
import scanpy as sc

# import seaborn as sns


def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(
        description="Generate quality control metrics for single-cell RNA-seq."
    )
    parser.add_argument("adata_file", type=str, help="Path to AnnData object file")
    parser.add_argument(
        "--sample_field",
        type=str,
        default="Sample_ID",
        help="Field in the obs for the sample identifier"
    )
    parser.add_argument(
        "--output_format",
        type=str,
        default="pdf",
        choices=["pdf", "png"],
        help="Output format of the plots (default: pdf)",
    )
    parser.add_argument(
        "--plot_size",
        type=float,
        nargs=2,
        metavar=("width", "height"),
        help="Size of the plots",
    )
    # add an argument for general plot title font size
    parser.add_argument(
        "--title_font_size",
        type=int,
        default=12,
        help="General plot title font size",
    )
    # add an argument for general plot label font size
    parser.add_argument(
        "--label_font_size",
        type=int,
        default=8,
        help="General plot label font size",
    )
    # add an argument for general plot legend font size
    parser.add_argument(
        "--legend_font_size",
        type=int,
        default=10,
        help="General plot legend font size",
    )
    # add an argument for the gene symbols field
    parser.add_argument(
        "--gene_symbols_field",
        type=str,
        default="gene_symbols",
        help="Field in the var for the gene symbols",
    )
    parser.add_argument(
        "--percent_mito_field",
        type=str,
        default="pct_counts_mito",
        help="Field in the obs for the percentage of mitochondrial genes",
    )
    parser.add_argument(
        "--percent_ribo_field",
        type=str,
        default="pct_counts_ribo",
        help="Field in the obs for the percentage of ribosomal genes",
    )
    parser.add_argument(
        "--ribo_field",
        type=str,
        default="ribo",
        help="Field in the var for marking ribosomal genes",
    )
    parser.add_argument(
        "--mito_field",
        type=str,
        default="mito",
        help="Field in the var for marking mitochondrial genes",
    )
    parser.add_argument(
        "--doublet_score_field",
        type=str,
        default="doublet_score",
        help="Field in the obs for the doublet score",
    )
    # add an argument for an embedding to plot the cells
    parser.add_argument(
        "--embedding",
        type=str,
        default=None,
        help="Embedding to plot the cells",
    )
    args = parser.parse_args()

    # Load AnnData object
    adata = sc.read(args.adata_file)

    # Set plot size if provided
    if args.plot_size:
        sc.settings.figsize = tuple(args.plot_size)

    # set scanpy general plot font size and output format
    sc.settings.set_figure_params(scanpy=True,
                                  fontsize=args.label_font_size,
                                  format=args.output_format)
    # disable FutureWarning
    import warnings

    warnings.simplefilter(action="ignore", category=FutureWarning)

    run_quality_control = False
    if "n_genes_by_counts" not in adata.obs.columns:
        run_quality_control = True
    if "n_counts" not in adata.obs.columns:
        run_quality_control = True

    qc_vars = []
    fields = [
        "n_genes_by_counts",
        "total_counts",
        args.percent_mito_field,
        args.percent_ribo_field
    ]
    # calculate mitochondrial genes if not provided
    if args.percent_mito_field not in adata.obs.columns:
        qc_vars.append(args.mito_field)
    # calculate ribo metrics if not provided
    if args.ribo_field not in adata.var.columns:
        # create a new column with the name args.ribo_field where genes that
        # have in the gene symbols field the pattern ^RP[SL] are
        # marked as true
        print(f"Creating {args.ribo_field} column")
        adata.var[args.ribo_field] = adata.var[args.gene_symbols_field].str.contains(
                "^RP[SL]"
            )
        print(f"Number of ribosomal genes: {adata.var[args.ribo_field].sum()}")
    if args.percent_ribo_field not in adata.obs.columns:
        qc_vars.append(args.ribo_field)

    print(f"Calculating QC metrics for {len(qc_vars)} variables")

    if len(qc_vars) > 0 or run_quality_control:
        sc.pp.calculate_qc_metrics(
            adata,
            qc_vars=qc_vars,
            log1p=True,
            inplace=True,
        )
        adata.obs["n_counts"] = adata.obs["total_counts"]
        adata.obs["n_genes"] = adata.obs["n_genes_by_counts"]
        adata.var["n_counts"] = adata.var["total_counts"]
        adata.var["n_cells"] = adata.var["n_cells_by_counts"]

    # Define thresholds
    high_umi_threshold = adata.obs['n_counts'].quantile(0.95)  # Top 5% most UMI counts
    low_umi_threshold = adata.obs['n_counts'].quantile(0.05)   # Bottom 5% least UMI counts
    high_mito_threshold = adata.obs[args.percent_mito_field].quantile(0.90) # Top 10% pct mitochondrial genes

    from sklearn.linear_model import LinearRegression
    from sklearn.preprocessing import PolynomialFeatures

    # Polynomial regression to account for curvature in the n_counts vs. n_genes relationship
    poly = PolynomialFeatures(degree=2)
    X_poly = poly.fit_transform(adata.obs[['n_counts']])
    model = LinearRegression()
    model.fit(X_poly, adata.obs['n_genes'])
    predicted_counts = model.predict(X_poly)

    # Calculate residuals
    residuals = adata.obs['n_genes'] - predicted_counts
    outlier_threshold = residuals.abs().quantile(0.95)  # Top 5% residuals as outliers

    # Initialize diagnosis column
    adata.obs['auto_diagnosis'] = 'Healthy'

    # Identify outliers
    outliers = residuals.abs() > outlier_threshold
    adata.obs.loc[outliers, 'auto_diagnosis'] = 'Outlier'


    # Identify stressed/dying/apoptotic cells
    stressed_cells = (adata.obs['n_counts'] > high_umi_threshold) & (adata.obs[args.percent_mito_field] > high_mito_threshold)
    adata.obs.loc[stressed_cells, 'auto_diagnosis'] = 'Stressed/Dying/Apoptotic'

    # Identify poor-quality cells
    poor_quality_cells = (adata.obs['n_counts'] < low_umi_threshold) & (adata.obs[args.percent_mito_field] > high_mito_threshold)
    adata.obs.loc[poor_quality_cells, 'auto_diagnosis'] = 'Poor-Quality'

    # Print diagnosis summary
    print(adata.obs['auto_diagnosis'].value_counts())
    # make a barplot of the auto_diagnosis, omitting the healthy cells from the plot
    # but writing the number of healthy cells in the title. Plot per sample
    healthy_cells = adata.obs['auto_diagnosis'] == 'Healthy'
    healthy_count = healthy_cells.sum()

    # General quality for whole dataset
    plt.figure()
    ax = sc.pl.violin(
        adata,
        fields,
        jitter=False,
        multi_panel=True,
        show=False,
    )
    # ax.set_title("General QC")
    plt.savefig(f"general.{args.output_format}", bbox_inches="tight")
    plt.close()

    # Generate quality control plots
    generate_violin_plots(
        adata, args.sample_field, args.percent_mito_field,
        args.percent_ribo_field, format=args.output_format
    )
    generate_scatter_plot(
        adata,
        args.sample_field,
        percent_mito_field=args.percent_mito_field,
    )
    generate_scatter_plot(
        adata,
        args.sample_field,
        y='log1p_n_genes_by_counts',
        percent_mito_field=args.percent_mito_field,
    )
    if args.doublet_score_field in adata.obs.columns:
        generate_doublet_plot(
            adata,
            args.sample_field,
            double_score_field=args.doublet_score_field,
            format=args.output_format,
        )
    else:
        print(
            "Doublet score field provided not in adata.obs.columns, " + "skipping plot."
        )
    generate_complexity_plot(adata, args.sample_field, format=args.output_format)

    if args.embedding:
        generate_embedding_plot(
            adata,
            fields=fields + [args.sample_field, 'auto_diagnosis'],
            embedding=args.embedding,
            format=args.output_format,
        )

    generate_barplot(adata[~healthy_cells],
                     groups_field=args.sample_field,
                     props_field='auto_diagnosis',
                     figure_path='diagnosis_barplot.pdf',
                     topic_for_title=f"(Total Healthy/Unhealthy cells: {healthy_count}/{adata.n_obs - healthy_count})")
    # generate_scatter_by_sample(
    #     adata,
    #     sample_field=args.sample_field,
    #     format=args.output_format,
    #     percent_mito_field=args.percent_mito_field,
    # )


def generate_barplot(
    adata, groups_field, props_field, figure_path=None, topic_for_title=None
):
    """
    Generate a proportional bar plot from an AnnData object.

    Parameters:
    adata (AnnData): The input AnnData object containing the data to plot.
    groups_field (str): The column in adata.obs to group the data by.
    props_field (str): The column in adata.obs to plot as proportions.
    figure_path (str, optional): The path to save the generated figure. If not provided, the figure is not saved.
    topic_for_title (str, optional): The topic to be used in the figure title, goes after {props_field} proportion of {topic_for_title} per {groups_field}.

    Returns:
    matplotlib.figure.Figure: The generated bar plot.
    """
    props_plot_data = adata.obs[[groups_field, props_field]]
    # props_plot_data[groups_field] = props_plot_data[groups_field].cat.reorder_categories(['control', '2 days', '7 days', '10 days', '14 days'])
    # make a 100% stacked bar plot of props_plot_data, plotting phase counts grouped by cell_line_persister

    grouped = props_plot_data.groupby([groups_field, props_field]).size().unstack()
    # proportions = grouped.div(grouped.sum(axis=1), axis=0)
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    plt.gca().set_prop_cycle(color=colors[2:5])
    grouped.plot(kind="bar", stacked=False, figsize=(8, 6))
    if topic_for_title is not None:
        plt.title(f"{props_field} cells {topic_for_title}\nper {groups_field}")
    else:
        plt.title(f"{props_field} cells\nper {groups_field}")
    plt.xlabel(groups_field)
    # plt.xticks(rotation=45)
    plt.ylabel("Number of cells")

    # save plot to PDF file
    if figure_path is not None:
        plt.savefig(figure_path, bbox_inches="tight")
    return plt.figure()


def generate_embedding_plot(
            adata,
            fields,
            embedding,
            format="pdf"
        ):
    # Embedding plot
    plt.figure()
    sc.pl.embedding(
        adata,
        basis=embedding,
        color=fields,
        show=False,
        ncols=1
    )
    plt.savefig(f"embedding_plots.{format}", bbox_inches="tight")
    plt.close()


def generate_violin_plots(
    adata,
    sample_field,
    percent_mito_field="percent_mito",
    percent_ribo_field="percent_ribo",
    format="pdf",
    gene_symbols_field="gene_symbols",
):
    # Number of counts per cell
    plt.figure()
    ax = plt.gca()
    sc.pl.violin(
        adata,
        "n_counts",
        ax=ax,
        groupby=sample_field,
        title="Number of Counts per Cell (Separated by Sample)",
        show=False,
    )
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
    plt.savefig(f"n_counts_per_cell.{format}", bbox_inches="tight")
    plt.close()

    # Number of genes per cell
    plt.figure()
    ax = plt.gca()
    sc.pl.violin(
        adata,
        "n_genes",
        groupby=sample_field,
        ax=ax,
        title="Number of Genes per Cell (Separated by Sample)",
        show=False,
        # show=True,
        # save="_n_genes_per_cell",
    )
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
    plt.savefig(f"n_genes_per_cell.{format}", bbox_inches="tight")
    plt.close()

    # Percentage of mitochondrial genes per cell
    plt.figure()
    ax = plt.gca()
    sc.pl.violin(
        adata,
        percent_mito_field,
        groupby=sample_field,
        ax=ax,
        title="Percentage of Mitochondrial " + "Genes per Cell (Separated by Sample)",
        show=False,
    )
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
    plt.savefig(f"percent_mito_per_cell.{format}", bbox_inches="tight")
    plt.close()

    # Percentage of ribosomal genes per cell
    plt.figure()
    ax = plt.gca()
    sc.pl.violin(
        adata,
        percent_ribo_field,
        groupby=sample_field,
        ax=ax,
        title="Percentage of Ribosomal " + "Genes per Cell (Separated by Sample)",
        show=False,
    )
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
    plt.savefig(f"percent_ribo_per_cell.{format}", bbox_inches="tight")
    plt.close()

    # highest expressed genes per cell
    plt.figure()
    ax = sc.pl.highest_expr_genes(
        adata, n_top=30, gene_symbols=gene_symbols_field, show=False
    )
    # set title of ax
    ax.set_title("Highest expressed genes per cell\nby Sample")
    plt.savefig(f"highest_expr_genes.{format}", bbox_inches="tight")
    plt.close()

    for sample in adata.obs[sample_field].unique():
        plt.figure()
        ax = sc.pl.highest_expr_genes(
            adata[adata.obs[sample_field] == sample],
            n_top=30,
            gene_symbols=gene_symbols_field,
            show=False,
        )
        ax.set_title(f"Highest expressed genes {sample}")
        # sanitise sample for filename
        sample_fn = sample.replace(" ", "_")
        # generate filename based on sample for plot
        plt.savefig(f"highest_expr_genes_{sample_fn}.{format}", bbox_inches="tight")
        plt.close()


def generate_scatter_plot(
    adata,
    sample_field,
    y="n_genes",
    percent_mito_field="percent_mito",
):
    # Scatter plot of UMIs vs genes detected
    plt.figure()
    sc.pl.scatter(
        adata,
        x="n_counts",
        y=y,
        color=sample_field,
        title="UMIs vs Genes Detected (by Sample)",
        save=f"_umi_vs_{y}_detected",
        show=False,
    )
    plt.close()

    # UMIs vs genes detected scatterplot, colored by mitochondrial gene ratio
    plt.figure()
    sc.pl.scatter(
        adata,
        x="n_counts",
        y=y,
        color=percent_mito_field,
        title="UMIs vs Genes Detected (by Mitochondrial Gene Ratio)",
        save=f"_umi_vs_{y}_detected_colored_by_mito",
        show=False,
    )
    plt.close()

    plt.figure()
    sc.pl.scatter(
        adata,
        x='n_counts',
        y=y,
        color='auto_diagnosis',
        title="UMIs vs Genes Detected (by Mitochondrial Gene Ratio)",
        save=f"_umi_vs_{y}_detected_colored_by_auto_diagnosis",
        show=False
        )
    plt.close()


def generate_scatter_by_sample(
    adata, sample_field, percent_mito_field="percent_mito", format="pdf"
):
    sample_ids = adata.obs[sample_field].unique()
    num_samples = len(sample_ids)

    plt.figure(figsize=(10, 6 * num_samples))

    for idx, sample_id in enumerate(sample_ids, 1):
        plt.subplot(num_samples, 1, idx)
        adata_sample = adata[adata.obs[sample_field] == sample_id]
        sc.pl.scatter(
            adata_sample,
            x="n_counts",
            y="n_genes",
            title=f"Sample {sample_id}: UMI vs Genes Detected",
            color=percent_mito_field,
            show=False,
        )
        plt.title(f"Sample {sample_id}: UMI vs Genes Detected")
        plt.xlabel("UMIs")
        plt.ylabel("Genes Detected")
    plt.savefig(f"n_counts_n_genes_by_sample.{format}")
    plt.close()


def generate_doublet_plot(
    adata,
    sample_field,
    double_score_field="doublet_score",
    format="pdf",
):
    # Ratio of doublets per cell
    plt.figure()
    ax = plt.gca()
    sc.pl.violin(
        adata,
        double_score_field,
        groupby=sample_field,
        ax=ax,
        title="Doublet score distribution (Separated by Sample)",
        # save="_doublet_ratio",
        show=False,
    )
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
    plt.savefig(f"doublet_ratio.{format}", bbox_inches="tight")
    plt.close()


def generate_complexity_plot(adata, sample_field, format="pdf"):
    # Complexity distribution (log10 Genes per UMI)
    plt.figure()
    ax = plt.gca()
    sc.pl.violin(
        adata,
        "log1p_n_genes_by_counts",
        groupby=sample_field,
        ax=ax,
        title="Complexity Distribution (Log10 Genes per UMI)",
        # save="_complexity_distribution",
        show=False,
    )
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
    plt.savefig(f"complexity_distribution.{format}", bbox_inches="tight")
    plt.close()


if __name__ == "__main__":
    main()
