import argparse
import os
import tempfile

import anndata
import decoupler as dc
import pandas as pd


def read_gmt(gmt_file):
    """
    Reads a GMT file into a Pandas DataFrame.

    Parameters
    ----------
    gmt_file : str
        Path to the GMT file.

    Returns
    -------
    pd.DataFrame
        A DataFrame with the gene sets. Each row represents a gene set, and the columns are "gene_set_name", "gene_set_url", and "genes".
    >>> line = "HALLMARK_NOTCH_SIGNALING\\thttp://www.gsea-msigdb.org/gsea/msigdb/human/geneset/HALLMARK_NOTCH_SIGNALING\\tJAG1\\tNOTCH3\\tNOTCH2\\tAPH1A\\tHES1\\tCCND1\\tFZD1\\tPSEN2\\tFZD7\\tDTX1\\tDLL1\\tFZD5\\tMAML2\\tNOTCH1\\tPSENEN\\tWNT5A\\tCUL1\\tWNT2\\tDTX4\\tSAP30\\tPPARD\\tKAT2A\\tHEYL\\tSKP1\\tRBX1\\tTCF7L2\\tARRB1\\tLFNG\\tPRKCA\\tDTX2\\tST3GAL6\\tFBXW11\\n"
    >>> line2 = "HALLMARK_APICAL_SURFACE\\thttp://www.gsea-msigdb.org/gsea/msigdb/human/geneset/HALLMARK_APICAL_SURFACE\\tB4GALT1\\tRHCG\\tMAL\\tLYPD3\\tPKHD1\\tATP6V0A4\\tCRYBG1\\tSHROOM2\\tSRPX\\tMDGA1\\tTMEM8B\\tTHY1\\tPCSK9\\tEPHB4\\tDCBLD2\\tGHRL\\tLYN\\tGAS1\\tFLOT2\\tPLAUR\\tAKAP7\\tATP8B1\\tEFNA5\\tSLC34A3\\tAPP\\tGSTM3\\tHSPB1\\tSLC2A4\\tIL2RB\\tRTN4RL1\\tNCOA6\\tSULF2\\tADAM10\\tBRCA1\\tGATA3\\tAFAP1L2\\tIL2RG\\tCD160\\tADIPOR2\\tSLC22A12\\tNTNG1\\tSCUBE1\\tCX3CL1\\tCROCC\\n"
    >>> temp_dir = tempfile.gettempdir()
    >>> temp_gmt = os.path.join(temp_dir, "temp_file.gmt")
    >>> with open(temp_gmt, "w") as f:
    ...   f.write(line)
    ...   f.write(line2)
    288
    380
    >>> df = read_gmt(temp_gmt)
    >>> df.shape[0]
    2
    >>> df.columns == ["gene_set_name", "genes"]
    array([ True,  True])
    >>> df.loc[df["gene_set_name"] == "HALLMARK_APICAL_SURFACE"].genes.tolist()[0].startswith("B4GALT1")
    True
    """
    # Read the GMT file into a list of lines
    with open(gmt_file, "r") as f:
        lines = f.readlines()

    # Create a list of dictionaries, where each dictionary represents a gene set
    gene_sets = []
    for line in lines:
        fields = line.strip().split("\t")
        gene_set = {"gene_set_name": fields[0], "genes": ",".join(fields[2:])}
        gene_sets.append(gene_set)

    # Convert the list of dictionaries to a DataFrame
    return pd.DataFrame(gene_sets)


def score_genes_aucell(
    adata: anndata.AnnData, gene_list: list, score_name: str, use_raw=False
):
    """Score genes using Aucell.

    Parameters
    ----------
    adata : anndata.AnnData
    gene_list : list
    score_names : str
    use_raw : bool, optional

    >>> import scanpy as sc
    >>> import decoupler as dc
    >>> adata = sc.datasets.pbmc68k_reduced()
    >>> gene_list = adata.var[adata.var.index.str.startswith("RP")].index.tolist()
    >>> score_genes_aucell(adata, gene_list, "ribosomal_aucell", use_raw=False)
    >>> "ribosomal_aucell" in adata.obs.columns
    True
    """
    # make a data.frame with two columns, geneset and gene_id, geneset filled with score_names and gene_id with gene_list, one row per element
    geneset_df = pd.DataFrame(
        {
            "gene_id": gene_list,
            "geneset": score_name,
        }
    )
    # run decoupler's run_aucell
    dc.run_aucell(
        adata, net=geneset_df, source="geneset", target="gene_id", use_raw=use_raw
    )
    # copy .obsm['aucell_estimate'] matrix columns to adata.obs using the column names
    adata.obs[score_name] = adata.obsm["aucell_estimate"][score_name]


def run_for_genelists(
    adata, gene_lists, score_names, use_raw=False, gene_symbols_field="gene_symbols"
):
    if len(gene_lists) == len(score_names):
        for gene_list, score_names in zip(gene_lists, score_names):
            genes = gene_list.split(",")
            ens_gene_ids = adata.var[adata.var[gene_symbols_field].isin(genes)].index
            score_genes_aucell(
                adata,
                ens_gene_ids,
                f"AUCell_{score_names}",
                use_raw,
            )
    else:
        raise ValueError(
            "The number of gene lists (separated by :) and score names (separated by :) must be the same"
        )


if __name__ == "__main__":
    # Create command-line arguments parser
    parser = argparse.ArgumentParser(description="Score genes using Aucell")
    parser.add_argument("--input_file", type=str, help="Path to input AnnData file")
    parser.add_argument("--output_file", type=str, help="Path to output file")
    parser.add_argument("--gmt_file", type=str, help="Path to GMT file", required=False)
    # add argument for gene sets to score
    parser.add_argument(
        "--gene_sets_to_score",
        type=str,
        required=False,
        help="Comma separated list of gene sets to score (the need to be in the gmt file)",
    )
    # add argument for gene list (comma separated) to score
    parser.add_argument(
        "--gene_lists_to_score",
        type=str,
        required=False,
        help="Comma separated list of genes to score. You can have more than one set of genes, separated by colon :",
    )
    # argument for the score name when using the gene list
    parser.add_argument(
        "--score_names",
        type=str,
        required=False,
        help="Name of the score column when using the gene list. You can have more than one set of score names, separated by colon :. It should be the same length as the number of gene lists.",
    )
    parser.add_argument(
        "--gene_symbols_field",
        type=str,
        help="Name of the gene symbols field in the AnnData object",
    )
    parser.add_argument("--use_raw", action="store_true", help="Use raw data")
    parser.add_argument(
        "--write_anndata", action="store_true", help="Write the modified AnnData object"
    )

    # Parse command-line arguments
    args = parser.parse_args()

    # Load input AnnData object
    adata = anndata.read_h5ad(args.input_file)

    if args.gene_sets_to_score is not None and args.gmt_file is not None:
        # Load MSigDB file in GMT format
        msigdb = read_gmt(args.gmt_file)

        gene_sets_to_score = args.gene_sets_to_score.split(",")
        # Score genes by their ensembl ids using the score_genes_aucell function
        for _, row in msigdb.iterrows():
            gene_set_name = row["gene_set_name"]
            if gene_set_name in gene_sets_to_score:
                genes = row["genes"].split(",")
                # Convert gene symbols to ensembl ids by using the columns gene_symbols and index in adata.var specific to the gene set
                ens_gene_ids = adata.var[
                    adata.var[args.gene_symbols_field].isin(genes)
                ].index
                score_genes_aucell(
                    adata, ens_gene_ids, f"AUCell_{gene_set_name}", args.use_raw
                )
    elif args.gene_lists_to_score is not None and args.score_names is not None:
        gene_lists = args.gene_lists_to_score.split(":")
        score_names = args.score_names.split(",")
        run_for_genelists(
            adata, gene_lists, score_names, args.use_raw, args.gene_symbols_field
        )

    # Save the modified AnnData object or generate a file with cells as rows and the new score_names columns
    if args.write_anndata:
        adata.write_h5ad(args.output_file)
    else:
        new_columns = [col for col in adata.obs.columns if col.startswith("AUCell_")]
        adata.obs[new_columns].to_csv(args.output_file, sep="\t", index=True)
