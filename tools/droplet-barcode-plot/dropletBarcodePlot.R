#!/usr/bin/env Rscript

# This script parses the GTF file to create a feature-wise annotation file with
# mitochondrial features flagged, to assist in annotation and QC of single-cell
# expression data analysis.

suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(gridExtra))
suppressPackageStartupMessages(require(DropletUtils))
suppressPackageStartupMessages(require(Matrix))

die <- function(message){
  write(message, stderr())
  q(status = 1)
}

option_list = list(
  make_option(
    c("-b", "--barcode-frequencies"),
    action = "store",
    default = NA,
    type = 'character',
    help = "Path to a two-column tab-delimited file, with barcodes in the first column and frequencies in the second (ignored if --mtx-matrix supplied)"
  ),
  make_option(
    c("-m", "--mtx-matrix"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Matrix-market format matrix file, with cells by column (overrides --barcode-frequencies if supplied)'
  ),  
  make_option(
    c("-r", "--cells-by-row"),
    action = "store_true",
    default = FALSE,
    type = 'logical',
    help = 'For use with --mtx-matrix: force interpretation of matrix to assume cells are by row, rather than by column (default)'
  ),
  make_option(
    c("-l", "--label"),
    action = "store",
    default = '',
    type = 'character',
    help = 'Label to use in plot'
  ),
  make_option(
    c("-d", "--density-bins"),
    action = "store",
    default = 50,
    type = 'numeric',
    help = "Number of bins used to calculate density plot"
  ),
  make_option(
    c("-y", "--roryk-multiplier"),
    action = "store",
    default = 1.5,
    type = 'numeric',
    help = "Above-baseline multiplier to calculate roryk threshold"
  ),
  make_option(
    c("-o", "--output-plot"),
    action = "store",
    default = 'barcode_plot.png',
    type = 'character',
    help = "File path for output plot"
  ),
  make_option(
    c("-t", "--output-thresholds"),
    action = "store",
    default = 'barcode_thresholds.txt',
    type = 'character',
    help = "File path for output file containing calculted thresholds"
  )
)

opt <- parse_args(OptionParser(option_list = option_list), convert_hyphens_to_underscores = TRUE)

# Process inputs dependent on what has been provided

if (is.na(opt$mtx_matrix)){
  if (is.na(opt$barcode_frequencies)){
    die('ERROR: must supply --mtx-matrix or --barcode-frequencies')
  }else if (! file.exists(opt$barcode_frequencies)){
    die(paste('ERROR: barcode frequencies file', opt$barcode_frequencies, 'does not exist'))
  }else{
    barcode_counts <- read.delim(opt$barcode_frequencies, header = FALSE)
  }
}else if (! file.exists(opt$mtx_matrix)){
  die(paste('ERROR: MTX matrix file', opt$mtx_matrix, 'does not exist'))
}else{
  result_matrix <- Matrix::readMM(opt$mtx_matrix)
  if (opt$cells_by_row){
    barcode_counts <- data.frame(V1 = 1:nrow(result_matrix), V2=Matrix::rowSums(result_matrix))
  }else{
    barcode_counts <- data.frame(V1 = 1:ncol(result_matrix), V2=Matrix::colSums(result_matrix))
  }
}

# Pick a cutoff on count as per https://github.com/COMBINE-lab/salmon/issues/362#issuecomment-490160480

pick_roryk_cutoff = function(bcs, above_baseline_multiplier = 1.5){
  bcs_hist = hist(log10(bcs), plot=FALSE, n=opt$density_bins)
  mids = bcs_hist$mids
  vals = bcs_hist$count
  wdensity = vals * (10^mids) / sum(vals * (10^mids))
  baseline <- median(wdensity)
  
  # Find highest density in upper half of barcode distribution
  
  peak <- which(wdensity == max(wdensity[((length(wdensity)+1)/2):length(wdensity)]))
  
  # Cutoff is the point before the peak at which density falls below the multiplier of baseline
  
  10^mids[max(which(wdensity[1:peak] < (above_baseline_multiplier*baseline)))]
}

# Plot densities 

barcode_density_plot = function(bcs, roryk_cutoff, knee, inflection, name = '   ') {
  bcs_hist = hist(log10(bcs), plot=FALSE, n=opt$density_bins)
  counts = bcs_hist$count
  mids = bcs_hist$mids
  y = counts * (10^mids) / sum(counts * (10^mids))
  qplot(y, 10^mids) + geom_point() + theme_bw() + ggtitle(name) + ylab('Count') + xlab ('Density') +
    geom_hline(aes(yintercept = roryk_cutoff, color = paste('roryk_cutoff =', length(which(bcs > roryk_cutoff)), 'cells'))) + 
    geom_hline(aes(yintercept = inflection, color = paste('dropletutils_inflection =', length(which(bcs > inflection)), 'cells'))) +
    geom_hline(aes(yintercept = knee, color = paste('dropletutils_knee =', length(which(bcs > knee)), 'cells'))) +
    scale_y_continuous(trans='log10') + theme(axis.title.y=element_blank()) + labs(color='Thresholds')
}  

# Plot a more standard barcode rank plot

barcode_rank_plot <- function(br.out, roryk_total_cutoff, knee, inflection, name='no name'){
  ggplot(data.frame(br.out), aes(x=rank, y=total)) + geom_line() + scale_x_continuous(trans='log10') + scale_y_continuous(trans='log10') + theme_bw() + 
    geom_hline(aes(yintercept = knee, color = 'dropletutils_knee')) + 
    geom_hline(aes(yintercept = inflection, color = 'dropletutils_inflection')) +
    geom_hline(aes(yintercept = roryk_total_cutoff, color = 'roryk_cutoff')) +
    ggtitle(name) + ylab('Count') + xlab('Rank') + theme(legend.position = "none")
}

# Sort barcodes by descending frequency

barcode_counts <- barcode_counts[order(barcode_counts$V2, decreasing = TRUE), ]

roryk_count_cutoff <- pick_roryk_cutoff(barcode_counts$V2, opt$roryk_multiplier)
  
# Run dropletUtils' barcodeRanks to get knee etc
br.out <- barcodeRanks(t(barcode_counts[,2,drop=FALSE]))
  
dropletutils_knee <- metadata(br.out)$knee
dropletutils_inflection <- metadata(br.out)$inflection

plot_label <- paste(format(nrow(barcode_counts), big.mark = ','), 'cell barcodes')
if ((! is.na(opt$label)) && opt$label != ''){
  plot_label <- paste0(opt$label, ': ', plot_label)
}
  
plots <- list(
  dropletutils = barcode_rank_plot(br.out, roryk_count_cutoff, dropletutils_knee, dropletutils_inflection, name = plot_label),
  roryk = barcode_density_plot(barcode_counts$V2, roryk_count_cutoff, dropletutils_knee, dropletutils_inflection, name = '   ')
)

# Create output plot
png(width = 1000, height = 600, file=opt$output_plot)
grid.arrange(plots$dropletutils, plots$roryk, nrow=1)
dev.off()

# Return calculated thresholds
write.table(data.frame(dropletutils_knee = dropletutils_knee, dropletutils_inflection = dropletutils_inflection, roryk=roryk_count_cutoff), file = opt$output_thresholds, row.names = FALSE, quote = FALSE)
