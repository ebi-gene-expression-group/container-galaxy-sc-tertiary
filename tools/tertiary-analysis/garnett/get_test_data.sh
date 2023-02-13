#!/usr/bin/env bash

BASENAME_FILE='E-GEOD-83139'

MTX_LINK='EOD-83139/E-GEOD-83139.aggregated_filtered_normalised_counts.mtx'
BARCODE_LINK='http://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/sc_experiments/E-GEOD-83139/E-GEOD-83139.aggregated_filtered_normalised_counts.mtx_cols'
GENE_LINK='http://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/sc_experiments/E-GEOD-83139/E-GEOD-83139.aggregated_filtered_normalised_counts.mtx_rows'

function get_data {
  local link=$1
  local fname=$2

  if [ ! -f $fname ]; then
    echo "$fname not available locally, downloading.."
    wget -O $fname --retry-connrefused --waitretry=1 --read-timeout=20 --timeout=15 -t 3 $link
  fi
}

# get matrix data
pushd test-data
get_data $MTX_LINK $BASENAME_FILE".aggregated_filtered_normalised_counts.mtx"
get_data $BARCODE_LINK $BASENAME_FILE".aggregated_filtered_normalised_counts.mtx_cols"
get_data $GENE_LINK $BASENAME_FILE".aggregated_filtered_normalised_counts.mtx_rows"
