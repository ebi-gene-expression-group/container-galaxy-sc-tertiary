#!/usr/bin/env bash
TF_DATA_LINK='https://raw.githubusercontent.com/aertslab/scenic-nf/master/example/allTFs_hg38.txt'
MOTIF2TF_LINK='https://raw.githubusercontent.com/aertslab/scenic-nf/master/example/motifs.tbl'
RANKING_LINK='https://zenodo.org/records/13328724/files/genome-ranking_v2.feather'
LOOM_INPUT_LINK='https://raw.githubusercontent.com/aertslab/scenic-nf/master/example/expr_mat.loom'

REGULONS_LINK='https://zenodo.org/records/13328724/files/regulons.tsv'
TF2TARGETS_LINK='https://zenodo.org/records/13328724/files/tf2targets.tsv'

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
get_data $TF_DATA_LINK "allTFs_hg38.txt"
get_data $MOTIF2TF_LINK "motifs.tbl"
get_data $RANKING_LINK "genome-ranking_v2.feather"
get_data $LOOM_INPUT_LINK "expr_mat.loom"
get_data $REGULONS_LINK regulons.tsv
get_data $TF2TARGETS_LINK tf2targets.tsv