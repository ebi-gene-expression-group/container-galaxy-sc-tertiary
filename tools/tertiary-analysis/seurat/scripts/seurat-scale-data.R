#!/usr/bin/env Rscript 

# Load optparse we need to check inputs

suppressPackageStartupMessages(require(optparse))

# Load common functions

suppressPackageStartupMessages(require(workflowscriptscommon))

# parse options

option_list = list(
  make_option(
    c("-i", "--input-object-file"),
    action = "store",
    default = NA,
    type = 'character',
    help = "File name in which a serialized R matrix object may be found."
  ),
  make_option(
    c("--input-format"),
    action = "store",
    default = "seurat",
    type = 'character',
    help = "Either loom, seurat, anndata or singlecellexperiment for the input format to read."
  ),
  make_option(
    c("--output-format"),
    action = "store",
    default = "seurat",
    type = 'character',
    help = "Either loom, seurat, anndata or singlecellexperiment for the output format."
  ),
  make_option(
    c("-e", "--genes-use"),
    action = "store",
    default = NULL,
    type = 'character',
    help = "File with gene names to scale/center (one gene per line). Default is all genes in object@data."
  ),
  make_option(
    c("-v", "--vars-to-regress"),
    action = "store",
    default = NULL,
    type = 'character',
    help = "Comma-separated list of variables to regress out (previously latent.vars in RegressOut). For example, nUMI, or percent.mito."
  ),
  make_option(
    c("-m", "--model-use"),
    action = "store",
    default = 'linear',
    type = 'character',
    help = "Use a linear model or generalized linear model (poisson, negative binomial) for the regression. Options are 'linear' (default), 'poisson', and 'negbinom'."
  ),
  make_option(
    c("-u", "--use-umi"),
    action = "store",
    default = FALSE,
    type = 'logical',
    help = "Regress on UMI count data. Default is FALSE for linear modeling, but automatically set to TRUE if model.use is 'negbinom' or 'poisson'."
  ),
  make_option(
    c("-s", "--do-not-scale"),
    action = "store_true",
    default = FALSE,
    type = 'logical',
    help = "Skip the data scale."
  ),
  make_option(
    c("-c", "--do-not-center"),
    action = "store_true",
    default = FALSE,
    type = 'logical',
    help = "Skip data centering."
  ),  
  make_option(
    c("-x", "--scale-max"),
    action = "store",
    default = 10,
    type = 'double',
    help = "Max value to return for scaled data. The default is 10. Setting this can help reduce the effects of genes that are only expressed in a very small number of cells. If regressing out latent variables and using a non-linear model, the default is 50."
  ),
  make_option(
    c("-b", "--block-size"),
    action = "store",
    default = 1000,
    type = 'integer',
    help = "Default size for number of genes to scale at in a single computation. Increasing block.size may speed up calculations but at an additional memory cost."
  ),
  make_option(
    c("-d", "--min-cells-to-block"),
    action = "store",
    default = 1000,
    type = 'integer',
    help = "If object contains fewer than this number of cells, don't block for scaling calculations."
  ),
  make_option(
    c("-n", "--check-for-norm"),
    action = "store",
    default = TRUE,
    type = 'logical',
    help = "Check to see if data has been normalized, if not, output a warning (TRUE by default)."
  ),
  make_option(
    c("-o", "--output-object-file"),
    action = "store",
    default = NA,
    type = 'character',
    help = "File name in which to store serialized R object of type 'Seurat'.'"
  )
)

opt <- wsc_parse_args(option_list, mandatory = c('input_object_file', 'output_object_file'))

# Check parameter values

if ( ! file.exists(opt$input_object_file)){
  stop((paste('File', opt$input_object_file, 'does not exist')))
}

if (! is.null(opt$genes_use)){
  if (! file.exists(opt$genes_use)){
    stop((paste('Supplied genes file', opt$genes_use, 'does not exist')))
  }else{
    genes_use <- readLines(opt$genes_use)
  }
}else{
  genes_use <- NULL
}

# break up opt$vars_to_regress into a list if it has commas
opt$vars_to_regress <- unlist(strsplit(opt$vars_to_regress, ","))

# Now we're hapy with the arguments, load Seurat and do the work

suppressPackageStartupMessages(require(Seurat))
if(opt$input_format == "loom" | opt$output_format == "loom") {
  suppressPackageStartupMessages(require(SeuratDisk))
} else if(opt$input_format == "singlecellexperiment" | opt$output_format == "singlecellexperiment") {
  suppressPackageStartupMessages(require(scater))
}

# Input from serialized R object

seurat_object <- read_seurat4_object(input_path = opt$input_object_file, format = opt$input_format)
# https://stackoverflow.com/questions/9129673/passing-list-of-named-parameters-to-function
# might be useful
scaled_seurat_object <- ScaleData(seurat_object, 
                                  features = genes_use, 
                                  vars.to.regress = opt$vars_to_regress, 
                                  model.use = opt$model_use, 
                                  use.umi = opt$use_umi, 
                                  do.scale = !opt$do_not_scale, 
                                  do.center = !opt$do_not_center, 
                                  scale.max = opt$scale_max, 
                                  block.size = opt$block_size, 
                                  min.cells.to.block = opt$min_cells_to_block, 
                                  verbose = FALSE)


# Output to a serialized R object
write_seurat4_object(seurat_object = scaled_seurat_object, 
                     output_path = opt$output_object_file,
                     format = opt$output_format)
