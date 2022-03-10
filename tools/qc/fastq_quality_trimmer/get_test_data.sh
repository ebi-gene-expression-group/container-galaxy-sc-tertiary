#!/usr/bin/env bash

BASE_LINK="https://raw.githubusercontent.com/agordon/fastx_toolkit/master/galaxy/test-data"

FQ_FILE="fastq_quality_trimmer.fastq"

FQ_LINK=$BASE_LINK"/"$FQ_FILE

function get_data {
  local link=$1
  local fname=$2

  if [ ! -f $fname ]; then
    echo "$fname not available locally, downloading.."
    wget -O $fname --retry-connrefused --waitretry=1 --read-timeout=20 --timeout=15 -t 3 $link
  fi
}

# Get test data
pushd test-data

get_data $FQ_LINK $FQ_FILE
