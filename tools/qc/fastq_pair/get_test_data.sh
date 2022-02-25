#!/usr/bin/env bash

BASE_LINK="https://raw.githubusercontent.com/linsalrob/fastq-pair/master/test"

LEFT_FILE="left.fastq"
RIGHT_FILE="right.fastq"
LEFTGZ_FILE="left.fastq.gz"
RIGHTGZ_FILE="right.fastq.gz"

LEFT_LINK=$BASE_LINK"/"$LEFT_FILE
RIGHT_LINK=$BASE_LINK"/"$RIGHT_FILE
LEFTGZ_LINK=$BASE_LINK"/"$LEFTGZ_FILE
RIGHTGZ_LINK=$BASE_LINK"/"$RIGHTGZ_FILE

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

get_data $LEFT_LINK $LEFT_FILE
get_data $RIGHT_LINK $RIGHT_FILE
get_data $LEFTGZ_LINK $LEFTGZ_FILE
get_data $RIGHTGZ_LINK $RIGHTGZ_FILE
