#!/usr/bin/env bash

BASENAME_FILE='E-MTAB-6077-3k_features_90_cells'

MTX_LINK='https://drive.google.com/uc?export=download&id=1-1ejn7scP80xsbrG0FtWzsozjg0hhc23'
RDS_LINK='https://drive.google.com/uc?export=download&id=1KW_GX6xznSUpWRWUykpNaSbAhyClf7_n'
NORM_LINK='https://drive.google.com/uc?export=download&id=1mvo3ENkBvEAOyWG6ejApzQTPDLX5yBKU'
FVG_LINK='https://drive.google.com/uc?export=download&id=13Fhruuj-vEEo1WM138ahtAYqfHc7LsaZ'
SCALED_LINK='https://drive.google.com/uc?export=download&id=18TK8us235LWNajarWDBAtASUXMYAxvw0'
PCA_LINK='https://drive.google.com/uc?export=download&id=1gf3BTB4dygDsom1TzjsBfgZnZepcoG5c'
NEIGHBOURS_LINK='https://drive.google.com/uc?export=download&id=1N2lHoKRBZ7pmAYGfghLWB9KUrLA5WoNX'
CLUSTERS_LINK='https://drive.google.com/uc?export=download&id=1HWxZWHbNUNo4z__9PhhL_CJOLzec_ETa'
TSNE_LINK='https://drive.google.com/uc?export=download&id=1qsvMr_GkCSp1dyTJt1BZ6cElJwFFX2zO'
MARKERS_LINK='https://drive.google.com/uc?export=download&id=18OmWNc7mF-4pzH6DQkPp1eKunN4BfvxD'

LOOM_LINK='https://drive.google.com/uc?export=download&id=1qNk5cg8hJG3Nv1ljTKmUEnxTOf11EEZX'
H5AD_LINK='https://drive.google.com/uc?export=download&id=1YpE0H_t_dkh17P-WBhPijKvRiGP0BlBz'
H5AD_SC182_LINK='https://drive.google.com/uc?export=download&id=16PUJ2KAkXT8F1UkfqU-9LWoOJUkUG1rp'
SCE_LINK='https://drive.google.com/uc?export=download&id=1UKdyf3M01uAt7oBg93JfmRvNVB_jlUKe'

# Seurat v4 exclusives
IFNB_BASE_FILE='ifnb_'

IFNB_CTRL_INT_LINK='https://drive.google.com/uc?export=download&id=15E_MLz-UclJYInNaA7YKLhLo5W-qlykL'
IFNB_STIM_INT_LINK='https://drive.google.com/uc?export=download&id=14iKgCJGPk16dEmpJJF-Gp_lBDcOdo-54'

## Classify and UMAP mapping
CLASSIFY_QUERY_LINK='https://drive.google.com/uc?export=download&id=1RFsHa_1EFD_n-19JH_cHGqxwO66QdmXN'
CLASSIFY_RESULTS_ANCHORS_OBJECT_LINK='https://drive.google.com/uc?export=download&id=1Xtv4K_CxIU1cJ8RjJ7NTvzLQkLvc8a3i'
UMAP_RESULT_OBJECT_LINK='https://oc.ebi.ac.uk/s/k4MdM07y9DAnurp/download'


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
get_data $MTX_LINK mtx.zip
unzip mtx.zip
rm -f mtx.zip

get_data $RDS_LINK $BASENAME_FILE".rds"
get_data $NORM_LINK $BASENAME_FILE"-normalised.rds"
get_data $FVG_LINK $BASENAME_FILE"-fvg.rds"
get_data $SCALED_LINK $BASENAME_FILE"-scaled.rds"
get_data $PCA_LINK $BASENAME_FILE"-pca.rds"
get_data $NEIGHBOURS_LINK $BASENAME_FILE"-neighbours.rds"
get_data $CLUSTERS_LINK $BASENAME_FILE"-clusters.rds"
get_data $TSNE_LINK $BASENAME_FILE"-tsne.rds"
get_data $MARKERS_LINK $BASENAME_FILE"-markers.csv.zip"

unzip $BASENAME_FILE"-markers.csv.zip"
rm -f $BASENAME_FILE"-markers.csv.zip"

get_data $LOOM_LINK $BASENAME_FILE"_loom.h5"
get_data $SCE_LINK $BASENAME_FILE"_sce.rds"
get_data $H5AD_LINK $BASENAME_FILE".h5ad"
get_data $H5AD_SC182_LINK $BASENAME_FILE"_sc182.h5ad"

get_data $IFNB_CTRL_INT_LINK $IFNB_BASE_FILE"ctrl_norm_fvg.rds"
get_data $IFNB_STIM_INT_LINK $IFNB_BASE_FILE"stim_norm_fvg.rds"

get_data $CLASSIFY_QUERY_LINK "Classify_query.rds"
get_data $CLASSIFY_RESULTS_ANCHORS_OBJECT_LINK "Classify_anchors.rds"
get_data $UMAP_RESULT_OBJECT_LINK "UMAP_result_integrated.rds"
