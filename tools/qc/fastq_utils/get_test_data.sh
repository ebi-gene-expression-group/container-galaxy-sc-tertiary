#!/usr/bin/env bash

BASE_LINK="https://raw.githubusercontent.com/nunofonseca/fastq_utils/master/tests"

BAR11_FILE="barcode_test_1.fastq.gz"
BAR12_FILE="barcode_test_2.fastq.gz"
BAR21_FILE="barcode_test2_1.fastq.gz"
BAR22_FILE="barcode_test2_2.fastq.gz"
INTER_FILE="inter.fastq.gz"

BAR11_LINK=$BASE_LINK"/"$BAR11_FILE
BAR12_LINK=$BASE_LINK"/"$BAR12_FILE
BAR21_LINK=$BASE_LINK"/"$BAR21_FILE
BAR22_LINK=$BASE_LINK"/"$BAR22_FILE
INTER_LINK=$BASE_LINK"/"$INTER_FILE

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

get_data $BAR11_LINK $BAR11_FILE
get_data $BAR12_LINK $BAR12_FILE
get_data $BAR21_LINK $BAR12_FILE
get_data $BAR22_LINK $BAR22_FILE
get_data $INTER_LINK $INTER_FILE
