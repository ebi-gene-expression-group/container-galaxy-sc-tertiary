#!/usr/bin/env bash

BASENAME_FILE='mito_counted_anndata.h5ad'

MTX_LINK='https://zenodo.org/record/7053673/files/Mito-counted_AnnData'

# convenience for getting data
function get_data {
  local link=$1
  local fname=$2

  if [ ! -f $fname ]; then
    echo "$fname not available locally, downloading.."
    wget -O $fname --retry-connrefused --waitretry=1 --read-timeout=20 --timeout=15 -t 3 $link
  fi
}

# get matrix data
mkdir -p test-data
pushd test-data
get_data $MTX_LINK $BASENAME_FILE


BASENAME_FILE='pbmc3k_processed.h5ad'

MTX_LINK='https://zenodo.org/records/3752813/files/pbmc3k_processed.h5ad'

get_data $MTX_LINK $BASENAME_FILE