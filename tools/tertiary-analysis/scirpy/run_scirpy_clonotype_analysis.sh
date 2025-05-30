#!/bin/bash

# Run scirpy_clonotype_analysis.py with GEX and VDJ data
# This script calls the scirpy_clonotype_analysis.py script with sample data

# Set variables for input files
GEX_FILE="test-data/10k_5p_Human_diseased_PBMC_ALL_Fresh_count_filtered_feature_bc_matrix_cell_annot.h5ad"
VDJ_FILE="test-data/10k_5p_Human_diseased_PBMC_ALL_Fresh_vdj_t_filtered_contig_annotations.csv"

# Set output directory
OUTPUT_DIR="./scirpy_output"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Run the analysis
python scirpy_clonotype_analysis.py \
    --gex "$GEX_FILE" \
    --vdj "$VDJ_FILE" \
    --vdj-format "10x" \
    --sample-name "PBMC_ALL" \
    --output-dir "$OUTPUT_DIR" \
    --min-cells 2 \
    --tcrdist-cutoff 15 \
    --plot-format "png" \
    --plot-dpi 120 \
    --output-tsv \
    --output-adata \
    --no-output-mudata \
    --metadata-obskey "gex:predicted_labels" \
    --perform-epitope-analysis

echo "Analysis complete. Results saved to $OUTPUT_DIR"
