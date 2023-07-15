#!/usr/bin/env bash
wget --retry-connrefused --waitretry=1 --read-timeout=20 --timeout=15 -t 3 http://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/sc_experiments/E-MTAB-6077/E-MTAB-6077.project.h5ad
mkdir -p test-data
mv E-MTAB-6077.project.h5ad test-data/
