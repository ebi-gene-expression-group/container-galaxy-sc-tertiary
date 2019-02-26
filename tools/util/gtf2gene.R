#!/usr/bin/env Rscript

suppressPackageStartupMessages(require(rtracklayer))
args <- commandArgs(TRUE)

annotation <- elementMetadata(import( args[1] ))
genes <- unique(annotation[['gene_id']])
writeLines(genes[ ! is.na(genes)], con = 'genes.txt')
