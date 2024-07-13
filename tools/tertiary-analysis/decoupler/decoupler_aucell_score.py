import argparse

import anndata
import decoupler as dc
import numba as nb
import pandas as pd


def read_gmt_long(gmt_file):
    """
    Reads a GMT file and produce a Pandas DataFrame in long format, ready to
    be passed to the AUCell method.

    Parameters
    ----------
    gmt_file : str
        Path to the GMT file.

    Returns
    -------
    pd.DataFrame
        A DataFrame with the gene sets. Each row represents a gene set to gene
        assignment, and the columns are "gene_set_name" and "genes".
    >>> line = "HALLMARK_NOTCH_SIGNALING\
    ...            \\thttp://www.gsea-msigdb.org/\
    ...            gsea/msigdb/human/geneset/HALLMARK_NOTCH_SIGNALING\
    ...            \\tJAG1\\tNOTCH3\\tNOTCH2\\tAPH1A\\tHES1\\tCCND1\
    ...            \\tFZD1\\tPSEN2\\tFZD7\\tDTX1\\tDLL1\\tFZD5\\tMAML2\
    ...            \\tNOTCH1\\tPSENEN\\tWNT5A\\tCUL1\\tWNT2\\tDTX4\
    ...            \\tSAP30\\tPPARD\\tKAT2A\\tHEYL\\tSKP1\\tRBX1\\tTCF7L2\
    ...            \\tARRB1\\tLFNG\\tPRKCA\\tDTX2\\tST3GAL6\\tFBXW11"
    >>> line2 = "HALLMARK_APICAL_SURFACE\
    ...            \\thttp://www.gsea-msigdb.org/\
    ...            gsea/msigdb/human/geneset/HALLMARK_APICAL_SURFACE\
    ...            \\tB4GALT1\\tRHCG\\tMAL\\tLYPD3\\tPKHD1\\tATP6V0A4\
    ...            \\tCRYBG1\\tSHROOM2\\tSRPX\\tMDGA1\\tTMEM8B\\tTHY1\
    ...            \\tPCSK9\\tEPHB4\\tDCBLD2\\tGHRL\\tLYN\\tGAS1\\tFLOT2\
    ...            \\tPLAUR\\tAKAP7\\tATP8B1\\tEFNA5\\tSLC34A3\\tAPP\
    ...            \\tGSTM3\\tHSPB1\\tSLC2A4\\tIL2RB\\tRTN4RL1\\tNCOA6\
    ...            \\tSULF2\\tADAM10\\tBRCA1\\tGATA3\\tAFAP1L2\\tIL2RG\
    ...            \\tCD160\\tADIPOR2\\tSLC22A12\\tNTNG1\\tSCUBE1\\tCX3CL1\
    ...            \\tCROCC\\n"
    >>> temp_dir = tempfile.gettempdir()
    >>> temp_gmt = os.path.join(temp_dir, "temp_file.gmt")
    >>> with open(temp_gmt, "w") as f:
    ...   f.write(line)
    ...   f.write(line2)
    288
    380
    >>> df = read_gmt_long(temp_gmt)
    >>> df.shape[0]
    76
    >>> len(df.loc[df["gene_set"] == "HALLMARK_APICAL_SURFACE"].gene.tolist())
    44
    """
    # Create a list of dictionaries, where each dictionary represents a
    # gene set
    gene_sets = {}

    # Read the GMT file into a list of lines
    with open(gmt_file, "r") as f:
        while True:
            line = f.readline()
            if not line:
                break
            fields = line.strip().split("\t")
            gene_sets[fields[0]] = fields[2:]

    return pd.concat(
        pd.DataFrame({"gene_set": k, "gene": v}) for k, v in gene_sets.items()
    )


def score_genes_aucell_mt(
    adata: anndata.AnnData,
    gene_set_gene: pd.DataFrame,
    use_raw=False,
    min_n_genes=5,
    var_gene_symbols_field=None,
):
    """Score genes using Aucell.

    Parameters
    ----------
    adata : anndata.AnnData
    gene_set_gene: pd.DataFrame with columns gene_set and gene
    use_raw : bool, optional, False by default.
    min_n_genes : int, optional, 5 by default.
    var_gene_symbols_field : str, optional, None by default. The field in var
    where gene symbols are stored

    >>> import scanpy as sc
    >>> import decoupler as dc
    >>> adata = sc.datasets.pbmc68k_reduced()
    >>> r_gene_list = adata.var[
    ...                  adata.var.index.str.startswith("RP")].index.tolist()
    >>> m_gene_list = adata.var[
    ...                  adata.var.index.str.startswith("M")].index.tolist()
    >>> gene_set = {}
    >>> gene_set["m"] = m_gene_list
    >>> gene_set["r"] = r_gene_list
    >>> gene_set_df = pd.concat(
    ...                  pd.DataFrame(
    ...                     {'gene_set':k, 'gene':v}
    ...     ) for k, v in gene_set.items())
    >>> score_genes_aucell_mt(adata, gene_set_df, use_raw=False)
    >>> "AUCell_m" in adata.obs.columns
    True
    >>> "AUCell_r" in adata.obs.columns
    True
    """

    # if var_gene_symbols_fiels is provided, transform gene_set_gene df so
    #  that gene contains gene ids instead of gene symbols
    if var_gene_symbols_field:
        # merge the index of var to gene_set_gene df based on
        # var_gene_symbols_field
        var_id_symbols = adata.var[[var_gene_symbols_field]]
        var_id_symbols["gene_id"] = var_id_symbols.index

        gene_set_gene = gene_set_gene.merge(
            var_id_symbols,
            left_on="gene",
            right_on=var_gene_symbols_field,
            how="left",
        )
        # this will still produce some empty gene_ids (genes in the
        # gene_set_gene df that are not in the var df), fill those
        # with the original gene symbol from the gene_set to avoid
        # deforming the AUCell calculation
        gene_set_gene["gene_id"] = gene_set_gene["gene_id"].fillna(
            gene_set_gene["gene"]
        )
        gene_set_gene["gene"] = gene_set_gene["gene_id"]

    # run decoupler's run_aucell
    dc.run_aucell(
        adata,
        net=gene_set_gene,
        source="gene_set",
        target="gene",
        use_raw=use_raw,
        min_n=min_n_genes,
    )
    for gs in gene_set_gene.gene_set.unique():
        if gs in adata.obsm["aucell_estimate"].keys():
            adata.obs[f"AUCell_{gs}"] = adata.obsm["aucell_estimate"][gs]


def run_for_genelists(
    adata,
    gene_lists,
    score_names,
    use_raw=False,
    gene_symbols_field=None,
    min_n_genes=5,
):
    if len(gene_lists) == len(score_names):
        for gene_list, score_names in zip(gene_lists, score_names):
            genes = gene_list.split(",")
            gene_sets = {}
            gene_sets[score_names] = genes
            gene_set_gene_df = pd.concat(
                pd.DataFrame({"gene_set": k, "gene": v})
                for k, v in gene_sets.items()
            )

            score_genes_aucell_mt(
                adata,
                gene_set_gene_df,
                use_raw,
                min_n_genes,
                var_gene_symbols_field=gene_symbols_field,
            )
    else:
        raise ValueError(
            "The number of gene lists (separated by :) and score names \
                (separated by :) must be the same"
        )


if __name__ == "__main__":
    # Create command-line arguments parser
    parser = argparse.ArgumentParser(description="Score genes using Aucell")
    parser.add_argument(
        "--input_file",
        type=str,
        help="Path to input AnnData file",
        required=True,
    )
    parser.add_argument(
        "--output_file", type=str, help="Path to output file", required=True
    )
    parser.add_argument(
        "--gmt_file", type=str, help="Path to GMT file", required=False
    )
    # add argument for gene sets to score
    parser.add_argument(
        "--gene_sets_to_score",
        type=str,
        required=False,
        help="Optional comma separated list of gene sets to score \
            (the need to be in the gmt file)",
    )
    # add argument for gene list (comma separated) to score
    parser.add_argument(
        "--gene_lists_to_score",
        type=str,
        required=False,
        help="Comma separated list of genes to score. You can have more \
            than one set of genes, separated by colon :",
    )
    # argument for the score name when using the gene list
    parser.add_argument(
        "--score_names",
        type=str,
        required=False,
        help="Name of the score column when using the gene list. You can \
            have more than one set of score names, separated by colon :. \
                It should be the same length as the number of gene lists.",
    )
    parser.add_argument(
        "--gene_symbols_field",
        type=str,
        help="Name of the gene symbols field in the AnnData object",
        required=True,
    )
    # argument for min_n Minimum of targets per source. If less, sources
    # are removed.
    parser.add_argument(
        "--min_n",
        type=int,
        required=False,
        default=5,
        help="Minimum of targets per source. If less, sources are removed.",
    )
    parser.add_argument("--use_raw", action="store_true", help="Use raw data")
    parser.add_argument(
        "--write_anndata",
        action="store_true",
        help="Write the modified AnnData object",
    )
    # argument for number of max concurrent processes
    parser.add_argument(
        "--max_threads",
        type=int,
        required=False,
        default=1,
        help="Number of max concurrent threads",
    )

    # Parse command-line arguments
    args = parser.parse_args()

    nb.set_num_threads(n=args.max_threads)

    # Load input AnnData object
    adata = anndata.read_h5ad(args.input_file)

    if args.gmt_file is not None:
        # Load MSigDB file in GMT format
        # msigdb = read_gmt(args.gmt_file)
        msigdb = read_gmt_long(args.gmt_file)

        gene_sets_to_score = (
            args.gene_sets_to_score.split(",")
            if args.gene_sets_to_score
            else []
        )
        if gene_sets_to_score:
            # we limit the GMT file read to the genesets specified in the
            # gene_sets_to_score argument
            msigdb = msigdb[msigdb["gene_set"].isin(gene_sets_to_score)]

        score_genes_aucell_mt(
            adata,
            msigdb,
            args.use_raw,
            args.min_n,
            var_gene_symbols_field=args.gene_symbols_field,
        )
    elif args.gene_lists_to_score is not None and args.score_names is not None:
        gene_lists = args.gene_lists_to_score.split(":")
        score_names = args.score_names.split(",")
        run_for_genelists(
            adata,
            gene_lists,
            score_names,
            args.use_raw,
            args.gene_symbols_field,
            args.min_n,
        )

    # Save the modified AnnData object or generate a file with cells as rows
    # and the new score_names columns
    if args.write_anndata:
        adata.write_h5ad(args.output_file)
    else:
        new_columns = [
            col for col in adata.obs.columns if col.startswith("AUCell_")
        ]
        adata.obs[new_columns].to_csv(args.output_file, sep="\t", index=True)
