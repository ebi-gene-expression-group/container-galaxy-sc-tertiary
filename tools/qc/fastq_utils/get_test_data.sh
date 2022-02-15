#!/usr/bin/env bash

BASE_LINK="https://raw.githubusercontent.com/nunofonseca/fastq_utils/master/tests"

BAR11_FILE="barcode_test_1.fastq.gz"
BAR12_FILE="barcode_test_2.fastq.gz"
BAR21_FILE="barcode_test2_1.fastq.gz"
BAR22_FILE="barcode_test2_2.fastq.gz"
INTER_FILE="inter.fastq.gz"
A1_FILE="a_1.fastq.gz"
POLYAT_FILE="poly_at.fastq.gz"
POLYAT3_FILE="poly_at_len3.fastq.gz"
TEST212_FILE="test_21_2.fastq.gz"

BAR11_LINK=$BASE_LINK"/"$BAR11_FILE
BAR12_LINK=$BASE_LINK"/"$BAR12_FILE
BAR21_LINK=$BASE_LINK"/"$BAR21_FILE
BAR22_LINK=$BASE_LINK"/"$BAR22_FILE
INTER_LINK=$BASE_LINK"/"$INTER_FILE
A1_LINK=$BASE_LINK"/"$A1_FILE
POLYAT_LINK=$BASE_LINK"/"$POLYAT_FILE
POLYAT3_LINK=$BASE_LINK"/"$POLYAT3_FILE
TEST212_LINK=$BASE_LINK"/"$TEST212_FILE

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
get_data $BAR21_LINK $BAR21_FILE
get_data $BAR22_LINK $BAR22_FILE
get_data $INTER_LINK $INTER_FILE
get_data $A1_LINK $A1_FILE
get_data $POLYAT_LINK $POLYAT_FILE
get_data $TEST212_LINK $TEST212_FILE
get_data $POLYAT3_LINK $POLYAT3_FILE
get_data $TEST212_LINK $TEST212_FILE
get_data $POLYAT3_LINK $POLYAT3_FILE
