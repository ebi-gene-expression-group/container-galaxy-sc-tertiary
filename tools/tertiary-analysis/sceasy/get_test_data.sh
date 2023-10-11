#!/usr/bin/env bash

mkdir -p test-data

BASENAME_FILE='E-MTAB-6077-3k_features_90_cells'

LOOM_LINK='https://drive.google.com/uc?export=download&id=1qNk5cg8hJG3Nv1ljTKmUEnxTOf11EEZX'
SCE_LINK='https://drive.google.com/uc?export=download&id=1UKdyf3M01uAt7oBg93JfmRvNVB_jlUKe'
H5AD_LINK='https://drive.google.com/uc?export=download&id=1YpE0H_t_dkh17P-WBhPijKvRiGP0BlBz'
RDS_LINK='https://drive.google.com/uc?export=download&id=1KW_GX6xznSUpWRWUykpNaSbAhyClf7_n'

function get_data {
  local link=$1
  local fname=$2

  if [ ! -f $fname ]; then
    echo "$fname not available locally, downloading.."
    wget -O $fname --retry-connrefused --waitretry=1 --read-timeout=20 --timeout=15 -t 3 $link
  fi
}

pushd test-data
get_data $LOOM_LINK $BASENAME_FILE"_loom.h5"
get_data $RDS_LINK $BASENAME_FILE".rds"
get_data $SCE_LINK $BASENAME_FILE"_sce.rds"
get_data $H5AD_LINK $BASENAME_FILE".h5ad"