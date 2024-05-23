import argparse

import matplotlib.pyplot as plt
import scanpy as sc
import seaborn as sns


def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(
        description="Generate quality control metrics for single-cell RNA-seq."
    )
    parser.add_argument("adata_file", type=str, help="Path to AnnData object file")
    parser.add_argument(
        "sample_field", type=str, help="Field in the obs for the sample identifier"
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
    args = parser.parse_args()

    # Load AnnData object
    adata = sc.read(args.adata_file)

    # Set plot size if provided
    if args.plot_size:
        sc.settings.figsize = tuple(args.plot_size)

    # Set output format
    sc.settings.set_figure_params(format=args.output_format)

    run_quality_control = False
    if "n_genes_by_counts" not in adata.obs.columns:
        run_quality_control = True
    if "n_counts" not in adata.obs.columns:
        run_quality_control = True

    qc_vars = []
    # calculate mitochondrial genes if not provided
    if args.percent_mito_field not in adata.obs.columns:
        qc_vars.append(args.mito_field)
    # calculate ribo metrics if not provided
    if args.percent_ribo_field not in adata.obs.columns:
        qc_vars.append(args.ribo_field)

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
    # General quality for whole dataset
    plt.figure()
    ax = sc.pl.violin(
        adata,
        [
            "n_genes_by_counts",
            "total_counts",
            args.percent_mito_field,
            args.percent_ribo_field,
        ],
        jitter=False,
        multi_panel=True,
        show=False,
    )
    # ax.set_title("General QC")
    plt.savefig(f"general.{args.output_format}", bbox_inches="tight")
    plt.close()

    # Generate quality control plots
    generate_violin_plots(
        adata, args.sample_field, args.percent_mito_field, format=args.output_format
    )
    generate_scatter_plot(
        adata,
        args.sample_field,
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
        print("Doublet score field provided not in adata.obs.columns, skipping plot.")
    generate_complexity_plot(adata, args.sample_field, format=args.output_format)
    # generate_scatter_by_sample(
    #     adata,
    #     sample_field=args.sample_field,
    #     format=args.output_format,
    #     percent_mito_field=args.percent_mito_field,
    # )


def generate_violin_plots(
    adata,
    sample_field,
    percent_mito_field="percent_mito",
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
        show=False
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
        title="Percentage of Mitochondrial Genes per Cell (Separated by Sample)",
        show=False,
    )
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
    plt.savefig(f"percent_mito_per_cell.{format}", bbox_inches="tight")
    plt.close()

    # highest expressed genes per cell
    plt.figure()
    ax = sc.pl.highest_expr_genes(
        adata, n_top=30, gene_symbols=gene_symbols_field, show=False
    )
    # set title of ax
    ax.set_title(f"Highest expressed genes per cell (Separated by Sample)")
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
    percent_mito_field="percent_mito",
):
    # Scatter plot of UMIs vs genes detected
    plt.figure()
    sc.pl.scatter(
        adata,
        x="n_counts",
        y="n_genes",
        color=sample_field,
        title="UMIs vs Genes Detected (Separated by Sample)",
        save="_umi_vs_genes_detected",
        show=False,
    )
    plt.close()

    # UMIs vs genes detected scatterplot, colored by mitochondrial gene ratio
    plt.figure()
    sc.pl.scatter(
        adata,
        x="n_counts",
        y="n_genes",
        color=percent_mito_field,
        title="UMIs vs Genes Detected (Colored by Mitochondrial Gene Ratio)",
        save="_umi_vs_genes_detected_colored_by_mito",
        show=False,
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
