#!/usr/bin/env Rscript

# This script parses the GTF file to create a feature-wise annotation file with
# mitochondrial features flagged, to assist in annotation and QC of single-cell
# expression data analysis.

suppressPackageStartupMessages(require(rtracklayer))
suppressPackageStartupMessages(require(optparse))

ucfirst <- function (str) {
  paste(toupper(substring(str, 1, 1)), tolower(substring(str, 2)), sep = "")
}

die <- function(message){
  write(message, stderr())
  q(status = 1)
}

cleanlist <- function(str){
  tolower(unlist(strsplit(str, ',')))
}

cl <- commandArgs(trailingOnly = TRUE)

option_list = list(
  make_option(
    c("-g", "--gtf-file"),
    action = "store",
    default = NA,
    type = 'character',
    help = "Path to a valid GTF file"
  ),
  make_option(
    c("-t", "--feature-type"),
    action = "store",
    default = 'gene',
    type = 'character',
    help = 'Feature type to use (default: gene)'
  ),
  make_option(
    c("-f", "--first-field"),
    action = "store",
    default = 'gene_id',
    type = 'character',
    help = 'Field to place first in output table (default: gene_id)'
  ),
  make_option(
    c("-r", "--no-header"),
    action = "store_false",
    default = TRUE,
    type = 'logical',
    help = 'Suppress header on output'
  ),
  make_option(
    c("-l", "--fields"),
    action = "store",
    default = NULL,
    type = 'character',
    help = 'Comma-separated list of output fields to retain (default: all)'
  ),
  make_option(
    c("-m", "--mito"),
    action = "store_true",
    default = FALSE,
    type = 'character',
    help = 'Mark mitochondrial elements with reference to chromsomes and biotypes'
  ),
  make_option(
    c("-n", "--mito-chr"),
    action = "store",
    default = 'mt,mitochondrion_genome,mito,m,chrM,chrMt',
    type = 'character',
    help = 'If specified, marks in a column called "mito" features on the specified chromosomes (case insensitive)'
  ),
  make_option(
    c("-p", "--mito-biotypes"),
    action = "store",
    default = 'mt_trna,mt_rrna,mt_trna_pseudogene',
    type = 'character',
    help = 'If specified,  marks in a column called "mito" features with the specified biotypes (case insensitve)'
  ),
  make_option(
    c("-c", "--filter-cdnas"),
    action = "store",
    default = NULL,
    type = 'character',
    help = 'If specified, sequences in the provided FASTA-format cDNAs file will be filtered to remove entries not present in the annotation'
  ),
  make_option(
    c("-d", "--filter-cdnas-field"),
    action = "store",
    default = 'transcript_id',
    type = 'character',
    help = 'Where --filter-cdnas is specified, what field should be used to compare to identfiers from the FASTA?'
  ),
  make_option(
    c("-e", "--filter-cdnas-output"),
    action = "store",
    default = 'filtered.fa.gz',
    type = 'character',
    help = 'Where --filter-cdnas is specified, what file should the filtered sequences be output to?'
  ),
  make_option(
    c("-u", "--version-transcripts"),
    action = "store_true",
    default = FALSE,
    type = 'logical',
    help = 'Where the GTF contains transcript versions, should these be appended to transcript identifiers? Useful when generating transcript/gene mappings for use with transcriptomes.'
  ),
  make_option(
    c("-o", "--output-file"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Output file path'
  )
)

opt <- parse_args(OptionParser(option_list = option_list), convert_hyphens_to_underscores = TRUE)

if (is.na(opt$gtf_file)){
  die('ERROR: No input GTF file specified')
}

if (is.na(opt$output_file)){
  die('ERROR: No output file specified')
}

# Import the GTF

print(paste('Reading', opt$gtf_file, 'elements of type', opt$feature_type))
gtf <- import(opt$gtf_file, feature.type = opt$feature_type )

# Combine basic info (chromosomes, coordinates) with annotation found in GTF attributes

anno <- cbind(chromosome = seqnames(gtf), as.data.frame(ranges(gtf)), elementMetadata(gtf))
print(paste('Found', nrow(anno), 'features'))

# Mark mitochondrial features

if (opt$mito){
  anno$mito <- ucfirst(as.character(tolower(anno$gene_biotype) %in% cleanlist(opt$mito_biotypes) | tolower(anno$chromosome) %in% cleanlist(opt$mito_chr)))
}

# If specified, put the desired field first

if (! is.na(opt$first_field)){
  if (! opt$first_field %in% colnames(anno)){
    die(paste(first_field, 'is not a valid field'))
  }
  anno <- anno[,c(opt$first_field, colnames(anno)[colnames(anno) != opt$first_field])]
}

# Version transcripts

if ( opt$feature_type == 'transcript' && opt$version_transcripts && all(c('transcript_id', 'transcript_version') %in% colnames(anno) )){
  has_transcript_version <- ! is.na(anno$transcript_version)
  anno$transcript_id[has_transcript_version] <- paste(anno$transcript_id[has_transcript_version], anno$transcript_version[has_transcript_version], sep='.')
}

# If specified, filter down a provided cDNA FASTA file

if (! is.null(opt$filter_cdnas)){
  
  print(paste("Filtering", opt$filter_cdnas, "to match the GTF"))
  
  suppressPackageStartupMessages(require(Biostrings))
  
  cdna <- readDNAStringSet(opt$filter_cdnas)
  cdna_transcript_names <- unlist(lapply(names(cdna), function(x) unlist(strsplit(x, ' '))[1]  ))
  
  # Filter out cDNAs without matching transcript entries in the GTF
  
  if (! any(cdna_transcript_names %in% anno[[opt$filter_cdnas_field]])){
    die(paste("ERROR: None of the input sequences have matching", opt$filter_cdnas_field, 'values in the GTF file'))
  }
  
  cdna <- cdna[which(cdna_transcript_names %in% anno[[opt$filter_cdnas_field]])]
  
  print(paste('Storing filtered seqeunces to', opt$filter_cdnas_output))
  writeXStringSet(x = cdna, filepath = opt$filter_cdnas_output, compress = 'gzip')
}

# If specified, subset to desired fields

if (! is.null(opt$fields) && opt$fields != ''){
  fields <- unlist(strsplit(opt$fields, ','))
  if (any(! fields %in% colnames(anno))){
    die(paste('ERROR:', fields, 'contains invalid field(s)'))
  }
  anno <- anno[,fields, drop = FALSE]
  anno <- anno[apply(anno, 1, function(x) all(! is.na(x))), ]
}

print(paste('Storing output to', opt$output_file))
write.table(anno, file = opt$output_file, sep = "\t", quote=FALSE, row.names = FALSE, col.names = opt$no_header)
