#!/usr/bin/env bash 

wget http://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/github_test_data/container-galaxy-sc-tertiary/cell-types-analysis.tar.gz -P test-data
pushd test-data
tar -zxvf cell-types-analysis.tar.gz 
mv cell-types-analysis/* ./ 
rm -r cell-types-analysis
popd 
