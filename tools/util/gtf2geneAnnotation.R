#!/usr/bin/env Rscript

# This script parses the GTF file to create a gene-wise annotation file with
# mitochondrial genes flagged, to assist in annotation and QC of single-cell
# expression data analysis.

suppressPackageStartupMessages(require(rtracklayer))

ucfirst <- function (str) {
  paste(toupper(substring(str, 1, 1)), tolower(substring(str, 2)), sep = "")
}

cl <- commandArgs(trailingOnly = TRUE)

gtf_file <- cl[1]
output <- cl[2]

gtf <- import(gtf_file, feature.type = 'gene')
anno <- cbind(chromosome = seqnames(gtf), as.data.frame(ranges(gtf)), elementMetadata(gtf))
anno$mito <- ucfirst(as.character(tolower(anno$gene_biotype) %in% c('mt_trna', 'mt_rrna', 'mt_trna_pseudogene') | tolower(anno$chromosome) %in% c('mt', 'mitochondrion_genome', 'mito')))

# Put the gene identifier first

anno <- anno[,c('gene_id', colnames(anno)[colnames(anno) != 'gene_id'])]

write.table(anno, file = output, sep = "\t", quote=FALSE, row.names = FALSE)